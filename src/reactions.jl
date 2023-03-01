"""
Convert SBML.Reaction to MTK.Reaction
"""
function get_reactions(model::SBML.Model)
    subsdict = get_substitutions(model)  # Todo: replace with SUBSDICT
    rxs = Reaction[]
    for reaction in values(model.reactions)
        extensive_math = SBML.extensive_kinetic_math(model, reaction.kinetic_math)
        symbolic_math = interpret_as_num(extensive_math)
        reactant_references = reaction.reactants
        product_references = reaction.products
        if reaction.reversible
            symbolic_math = get_unidirectional_components(symbolic_math)
            kl_fw, kl_rv = [substitute(x, subsdict) for x in symbolic_math]
            enforce_rate = isequal(kl_rv, 0)
            add_reaction!(rxs, kl_fw, reactant_references, product_references, model;
                          enforce_rate = enforce_rate)
            add_reaction!(rxs, kl_rv, product_references, reactant_references, model;
                          enforce_rate = enforce_rate)
        else
            kl = substitute(symbolic_math, subsdict)  # Todo: use SUBSDICT
            add_reaction!(rxs, kl, reactant_references, product_references, model)
        end
    end
    rxs
end

"""
Infer forward and reverse components of bidirectional kineticLaw
"""
function get_unidirectional_components(bidirectional_math)
    bm = ModelingToolkit.value(bidirectional_math)  #  Symbolics.tosymbol(bidirectional_math)
    bm = simplify(bm; expand = true)
    if !SymbolicUtils.isadd(bm)
        @warn "Cannot separate bidirectional kineticLaw `$bidirectional_math` to forward and reverse part. Setting forward to `$bidirectional_math` and reverse to `0`. Stochastic simulations will be inexact."
        return (bidirectional_math, Num(0))
    end
    terms = SymbolicUtils.arguments(ModelingToolkit.value(bm))
    fw_terms = []
    rv_terms = []
    for term in terms
        if SymbolicUtils.ismul(ModelingToolkit.value(term)) && (term.coeff < 0)
            push!(rv_terms, Num(-term))  # PL: @Anand: Perhaps we should to create_var(term) or so?
        else
            push!(fw_terms, Num(term))  # PL: @Anand: Perhaps we should to create_var(term) or so?
        end
    end
    if (length(fw_terms) != 1) || (length(rv_terms) != 1)
        @warn "Cannot separate bidirectional kineticLaw `$bidirectional_math` to forward and reverse part. Setting forward to `$bidirectional_math` and reverse to `0`. Stochastic simulations will be inexact."
        return (bidirectional_math, Num(0))
    end
    return (fw_terms[1], rv_terms[1])
end

function add_reaction!(rxs::Vector{Reaction},
                       kl::Num,
                       reactant_references::Vector{SBML.SpeciesReference},
                       product_references::Vector{SBML.SpeciesReference},
                       model::SBML.Model;
                       enforce_rate = false)
    reactants, products, rstoichvals, pstoichvals = get_reagents(reactant_references,
                                                                 product_references, model)
    isnothing(reactants) && isnothing(products) && return
    rstoichvals = stoich_convert_to_ints(rstoichvals)
    pstoichvals = stoich_convert_to_ints(pstoichvals)
    kl, our = use_rate(kl, reactants, rstoichvals)
    our = enforce_rate ? true : our
    push!(rxs,
          Catalyst.Reaction(kl, reactants, products, rstoichvals, pstoichvals;
                            only_use_rate = our))
end

function stoich_convert_to_ints(xs)
    (xs !== nothing && all(isinteger(x) for x in xs)) ? Int.(xs) : xs
end

"""
Get reagents
"""
function get_reagents(reactant_references::Vector{SBML.SpeciesReference},
                      product_references::Vector{SBML.SpeciesReference},
                      model::SBML.Model)
    reactants = Num[]
    products = Num[]
    rstoich = Float64[]
    pstoich = Float64[]

    for rr in reactant_references
        sn = rr.species
        stoich = rr.stoichiometry
        if isnothing(stoich)
            @warn "Stoichiometries of SpeciesReferences are not defined. Setting to 1." maxlog=1
            stoich = 1.0
        end
        iszero(stoich) && @error("Stoichiometry of $sn must be non-zero")
        push!(reactants, create_var(sn, IV))
        push!(rstoich, stoich)
        if model.species[sn].boundary_condition == true
            push!(products, create_var(sn, IV))
            push!(pstoich, stoich)
        end
    end
    for pr in product_references
        sn = pr.species
        stoich = pr.stoichiometry
        if isnothing(stoich)
            @warn "Stoichiometries of SpeciesReferences are not defined. Setting to 1." maxlog=1
            stoich = 1.0
        end
        iszero(stoich) && @error("Stoichiometry of $sn must be non-zero")
        if model.species[sn].boundary_condition != true
            push!(products, create_var(sn, IV))
            push!(pstoich, stoich)
        end
    end

    if (length(reactants) == 0)
        reactants = nothing
        rstoich = nothing
    end
    if (length(products) == 0)
        products = nothing
        pstoich = nothing
    end
    (reactants, products, rstoich, pstoich)
end

"""
Get kineticLaw for use in MTK.Reaction
"""
function use_rate(kl::Num, react::Union{Vector{Num}, Nothing},
                  stoich::Union{Vector{<:Real}, Nothing})
    rate_const = get_massaction(kl, react, stoich)
    if !isnan(rate_const)
        kl = rate_const
        our = false
    else
        our = true
    end
    return (kl, our)
end

"""
Get rate constant of mass action kineticLaws
"""
function get_massaction(kl::Num, reactants::Union{Vector{Num}, Nothing},
                        stoich::Union{Vector{<:Real}, Nothing})
    function check_args(x::SymbolicUtils.BasicSymbolic{Real})
        check_args(Val(SymbolicUtils.istree(x)), x)
    end
    function check_args(::Val{true}, x::SymbolicUtils.BasicSymbolic{Real})
        for arg in SymbolicUtils.arguments(x)
            if isnan(check_args(arg)) || isequal(arg, Catalyst.DEFAULT_IV)
                return NaN
            end
        end
        return 0
    end
    check_args(::Val{false}, x::SymbolicUtils.BasicSymbolic{Real}) = isspecies(x) ? NaN : 0  # Species or Parameter leaf node
    check_args(::Real) = 0  # Real leaf node
    check_args(x) = throw(ErrorException("Cannot handle $(typeof(x)) types."))  # Unknow leaf node
    if isnothing(reactants) && isnothing(stoich)
        rate_const = kl
    elseif isnothing(reactants) | isnothing(stoich)
        throw(ErrorException("`reactants` and `stoich` are incosistent: `reactants` are $(reactants) and `stoich` is $(stoich)."))
    else
        rate_const = SymbolicUtils.simplify_fractions(kl / *((.^(reactants, stoich))...))
    end
    isnan(check_args(ModelingToolkit.value(rate_const))) ? NaN : rate_const
end
