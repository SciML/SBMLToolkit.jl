""" ReactionSystem constructor """
function ModelingToolkit.ReactionSystem(model::Model; kwargs...)  # Todo: requires unique parameters (i.e. SBML must have been imported with localParameter promotion in libSBML)
    checksupport(model)
    model = make_extensive(model)
    rxs = mtk_reactions(model)
    species = []
    for k in keys(model.species)
        push!(species, create_var(k,Catalyst.DEFAULT_IV))
    end
    params = vcat([create_param(k) for k in keys(model.parameters)], [create_param(k) for (k,v) in model.compartments if !isnothing(v.size)])
    ReactionSystem(rxs,Catalyst.DEFAULT_IV,species,params; kwargs...)
end

""" ReactionSystem constructor """
function ModelingToolkit.ReactionSystem(sbmlfile::String; kwargs...)
    model = readSBML(sbmlfile,SBML.convert_simplify_math)
    ReactionSystem(model; kwargs...)
end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(model::Model; kwargs...)
    rs = ReactionSystem(model; kwargs...)
    model = make_extensive(model)  # PL: consider making `make_extensive!` to avoid duplicate calling in ReactionSystem and here
    u0map = get_u0(model)
    parammap = get_paramap(model)
    defaults = Dict(vcat(u0map, parammap))
    convert(ODESystem, rs, defaults=defaults)
end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(sbmlfile::String; kwargs...)
    model = readSBML(sbmlfile)
    ODESystem(model; kwargs...)
end

""" ODEProblem constructor """
function ModelingToolkit.ODEProblem(model::Model,tspan;kwargs...)  # PL: Todo: add u0 and parameters argument
    odesys = ODESystem(model)
    ODEProblem(odesys, [], tspan; kwargs...)
end

""" ODEProblem constructor """
function ModelingToolkit.ODEProblem(sbmlfile::String,tspan;kwargs...)  # PL: Todo: add u0 and parameters argument
    odesys = ODESystem(sbmlfile)
    ODEProblem(odesys, [], tspan; kwargs...)
end

""" Check if conversion to ReactionSystem is possible """
function checksupport(model::Model)
    for s in values(model.species)
        if s.boundary_condition
            @warn "Species $(s.name) has `boundaryCondition` or is `constant`. This will lead to wrong results when simulating the `ReactionSystem`."
        end
    end
end

""" Convert intensive to extensive expressions """
function make_extensive(model)
    model = to_initial_amounts(model)
    model = to_extensive_math!(model)
    model
end

""" Convert initial_concentration to initial_amount """
function to_initial_amounts(model::Model)
    model = deepcopy(model)
    for specie in values(model.species)
        if isnothing(specie.initial_amount)
            compartment = model.compartments[specie.compartment]
            if !isnothing(compartment.size)
                specie.initial_amount = (specie.initial_concentration[1] * compartment.size, "")
            else
                @warn "Compartment $(compartment.name) has no `size`. Cannot calculate `initial_amount` of Species $(specie.name). Setting `initial_amount` to `initial_concentration`."
                specie.initial_amount = (specie.initial_concentration[1], "")
            end
            specie.initial_concentration = nothing
        end
    end
    model
end

""" Convert intensive to extensive mathematical expression """
function to_extensive_math!(model::SBML.Model)
    function conv(x::SBML.MathApply)
        SBML.MathApply(x.fn, SBML.Math[conv(x) for x in x.args])
    end
    function conv(x::SBML.MathIdent)
        x_new = x
        if x.id in keys(model.species)
            specie = model.species[x.id]
            if !specie.only_substance_units
                compartment = model.compartments[specie.compartment]
                if isnothing(compartment.size)
                    @warn "Specie $(x.id) hasOnlySubstanceUnits but its compartment $(compartment.name) has no size. Cannot auto-correct the rate laws $(x.id) is involved in. Please check manually."
                else
                    x_new = SBML.MathApply("/", SBML.Math[x,
                                SBML.MathVal(compartment.size)])
                    specie.only_substance_units = true
                end
            end
        end
        x_new
    end
    conv(x::SBML.MathVal) = x
    conv(x::SBML.MathLambda) =
        throw(DomainError(x, "can't translate lambdas to extensive units"))
    for reaction in values(model.reactions)
        reaction.kinetic_math = conv(reaction.kinetic_math)
    end
    model
end

""" Get dictonary to change types in kineticLaw """
function _get_substitutions(model)
    subsdict = Dict()
    for k in keys(model.species)
        push!(subsdict, Pair(create_var(k),create_var(k,Catalyst.DEFAULT_IV)))
    end
    for k in keys(model.parameters)
        push!(subsdict, Pair(create_var(k),create_param(k)))
    end
    for k in keys(model.compartments)
        push!(subsdict, Pair(create_var(k),create_param(k)))
    end
    subsdict
end

""" Convert SBML.Reaction to MTK.Reaction """
function mtk_reactions(model::Model)
    subsdict = _get_substitutions(model)
    rxs = []
    if length(model.reactions) == 0
        throw(ErrorException("Model contains no reactions."))
    end
    for reaction in values(model.reactions)
        reactants = Num[]
        rstoich = Num[]
        products = Num[]
        pstoich = Num[]
        for (k,v) in reaction.stoichiometry
            if v < 0
                push!(reactants, create_var(k,Catalyst.DEFAULT_IV))
                push!(rstoich, -v)
            elseif v > 0
                push!(products, create_var(k,Catalyst.DEFAULT_IV))
                push!(pstoich, v)
            else
                @error("Stoichiometry of $k must be non-zero")
            end
        end
        if (length(reactants)==0) reactants = nothing; rstoich = nothing end
        if (length(products)==0) products = nothing; pstoich = nothing end
        symbolic_math = convert(Num, reaction.kinetic_math)
        
        if reaction.reversible
            symbolic_math = getunidirectionalcomponents(symbolic_math)
            kl = [substitute(x, subsdict) for x in symbolic_math]
            push!(rxs, ModelingToolkit.Reaction(kl[1],reactants,products,rstoich,pstoich;only_use_rate=true))
            push!(rxs, ModelingToolkit.Reaction(kl[2],products,reactants,pstoich,rstoich;only_use_rate=true))
        else
            kl = substitute(symbolic_math, subsdict)
            push!(rxs, ModelingToolkit.Reaction(kl,reactants,products,rstoich,pstoich;only_use_rate=true))
        end
    end
    rxs
end

""" Infer forward and reverse components of bidirectional kineticLaw """
function getunidirectionalcomponents(bidirectional_math)
    err = "Cannot separate bidirectional kineticLaw `$bidirectional_math` to forward and reverse part. Please make reaction irreversible or rearrange kineticLaw to the form `term1 - term2`."
    bidirectional_math = Symbolics.tosymbol(bidirectional_math)
    bidirectional_math = simplify(bidirectional_math; expand=true)
    if SymbolicUtils.operation(bidirectional_math) != +
        throw(ErrorException(err))
    end
    terms = SymbolicUtils.arguments(bidirectional_math)
    fw_terms = []
    rv_terms = []
    for term in terms
        if (term isa SymbolicUtils.Mul) && (term.coeff < 0)
            push!(rv_terms, Num(-term))  # PL: @Anand: Perhaps we should to create_var(term) or so?
        else
            push!(fw_terms, Num(term))  # PL: @Anand: Perhaps we should to create_var(term) or so?
        end
    end
    if (length(fw_terms) != 1) || (length(rv_terms) != 1)
        throw(ErrorException(err))
    end
    return (fw_terms[1], rv_terms[1])
end

""" Extract u0map from Model """
function get_u0(model)
    u0map = []
    for (k,v) in model.species
        push!(u0map,Pair(create_var(k,Catalyst.DEFAULT_IV), v.initial_amount[1]))
    end
    u0map
end

""" Extract paramap from Model """
function get_paramap(model)
    paramap = Pair{Num, Float64}[]
    for (k,v) in model.parameters
        push!(paramap,Pair(create_param(k),v))
    end
    for (k,v) in model.compartments
        if !isnothing(v.size)
            push!(paramap,Pair(create_param(k),v.size))
        end
    end
    paramap
end

create_var(x, iv) = Num(Variable{Symbolics.FnType{Tuple{Any},Real}}(Symbol(x)))(iv).val
create_var(x) = Num(Variable(Symbol(x))).val
function create_param(x)
    p = Sym{Real}(Symbol(x))
    ModelingToolkit.toparam(p)
end
