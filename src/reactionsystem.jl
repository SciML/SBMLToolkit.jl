""" ReactionSystem constructor """
function Catalyst.ReactionSystem(model::SBML.Model; kwargs...)  # Todo: requires unique parameters (i.e. SBML must have been imported with localParameter promotion in libSBML)
    rxs = mtk_reactions(model)
    u0map = get_u0map(model)
    parammap = get_paramap(model)
    defs = ModelingToolkit._merge(Dict(u0map), Dict(parammap))

    algrules, obsrules, raterules = get_rules(model)
    for o in obsrules
        defs[o.lhs] = substitute(o.rhs, defs)
        # push!(defs, o.lhs => substitute(o.rhs, defs))
    end
    constraints_sys = ODESystem(vcat(algrules, raterules), Catalyst.DEFAULT_IV; name = gensym(:CONSTRAINTS))

    ReactionSystem(rxs, Catalyst.DEFAULT_IV, first.(u0map), first.(parammap); defaults = defs, name = gensym(:SBML),
        observed = obsrules, constraints = constraints_sys, kwargs...)
end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(model::SBML.Model; include_zero_odes = false, kwargs...)
    rs = ReactionSystem(model; kwargs...)
    sys = convert(ODESystem, rs; include_zero_odes = include_zero_odes)
    ODESystem(equations(sys);
        name = nameof(sys),
        defaults = ModelingToolkit.get_defaults(rs),
        continuous_events = get_events(model, rs))
end

""" Check if conversion to ReactionSystem is possible """
function checksupport(filename::String)
    not_implemented = ["listOfConstraints"]
    sbml = open(filename) do file
        read(file, String)
    end
    for item in not_implemented
        occursin(item, sbml) && throw(ErrorException("SBML models with $item are not yet implemented."))
    end
    occursin("<sbml xmlns:fbc=", sbml) && throw(ErrorException("This model was designed for constrained-based optimisation. Please use COBREXA.jl instead of SBMLToolkit."))
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
        push!(subsdict, Pair(create_var(k), create_var(k, Catalyst.DEFAULT_IV)))
    end
    for k in keys(model.parameters)
        push!(subsdict, Pair(create_var(k), create_param(k)))
    end
    for k in keys(model.compartments)
        push!(subsdict, Pair(create_var(k), create_param(k)))
    end
    subsdict
end

""" Convert SBML.Reaction to MTK.Reaction """
function mtk_reactions(model::SBML.Model)
    subsdict = _get_substitutions(model)
    rxs = []
    # if length(model.reactions) == 0
    #     throw(ErrorException("SBML.Model contains no reactions."))
    # end
    for reaction in values(model.reactions)
        extensive_math = SBML.extensive_kinetic_math(
            model, reaction.kinetic_math,
            handle_empty_compartment_size = _ -> 1.0)
        symbolic_math = Num(convert(Num, extensive_math,
            convert_time = (x::SBML.MathTime) -> Catalyst.DEFAULT_IV))

        rstoich = reaction.reactants
        pstoich = reaction.products
        if reaction.reversible
            symbolic_math = getunidirectionalcomponents(symbolic_math)
            kl_fw, kl_rv = [substitute(x, subsdict) for x in symbolic_math]
            reagents = getreagents(rstoich, pstoich, model)
            isnothing(reagents[1]) && isnothing(reagents[2]) && continue
            kl_fw, our = use_rate(kl_fw, reagents[1], reagents[3])
            kl_rv = from_noncombinatoric(kl_rv, reagents[3], our)
            push!(rxs, Catalyst.Reaction(kl_fw, reagents...; only_use_rate = our))

            reagents = getreagents(rstoich, pstoich, model; rev = true)
            kl_rv, our = use_rate(kl_rv, reagents[1], reagents[3])
            kl_rv = from_noncombinatoric(kl_rv, reagents[3], our)
            push!(rxs, Catalyst.Reaction(kl_rv, reagents...; only_use_rate = our))
        else
            kl = substitute(symbolic_math, subsdict)
            reagents = getreagents(rstoich, pstoich, model)
            isnothing(reagents[1]) && isnothing(reagents[2]) && continue
            kl, our = use_rate(kl, reagents[1], reagents[3])
            kl = from_noncombinatoric(kl, reagents[3], our)
            push!(rxs, Catalyst.Reaction(kl, reagents...; only_use_rate = our))
        end
    end
    rxs
end

function from_noncombinatoric(rl::Num, stoich::Union{Vector{<:Real},Nothing}, only_use_rate::Bool)
    if !isnothing(stoich) && !only_use_rate
        coef = 1
        for s in stoich
            isone(s) && continue
            coef *= factorial(s)
        end
        !isone(coef) && (rl *= coef)
    end
    rl
end

""" Get kineticLaw for use in MTK.Reaction """
function use_rate(kl::Num, react::Union{Vector{Num},Nothing}, stoich::Union{Vector{<:Real},Nothing})
    rate_const = getmassaction(kl, react, stoich)
    if !isnan(rate_const)
        kl = rate_const
        our = false
    else
        our = true
    end
    return (kl, our)
end

""" Get reagents """
function getreagents(rstoichdict::Dict{String,<:Real}, pstoichdict::Dict{String,<:Real}, model::SBML.Model; rev = false)
    reactants = Num[]
    products = Num[]
    rstoich = Float64[]
    pstoich = Float64[]

    if rev
        tmp = rstoichdict
        rstoichdict = pstoichdict
        pstoichdict = tmp
    end

    for (k, v) in rstoichdict
        iszero(v) && @error("Stoichiometry of $k must be non-zero")
        push!(reactants, create_var(k, Catalyst.DEFAULT_IV))
        push!(rstoich, v)
        if model.species[k].boundary_condition == true
            push!(products, create_var(k, Catalyst.DEFAULT_IV))
            push!(pstoich, v)
        end
    end

    for (k, v) in pstoichdict
        iszero(v) && @error("Stoichiometry of $k must be non-zero")
        if model.species[k].boundary_condition != true
            push!(products, create_var(k, Catalyst.DEFAULT_IV))
            push!(pstoich, v)
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

""" Infer forward and reverse components of bidirectional kineticLaw """
function getunidirectionalcomponents(bidirectional_math)
    err = "Cannot separate bidirectional kineticLaw `$bidirectional_math` to forward and reverse part. Please make reaction irreversible or rearrange kineticLaw to the form `term1 - term2`."
    bidirectional_math = Symbolics.tosymbol(bidirectional_math)
    bidirectional_math = simplify(bidirectional_math; expand = true)
    if (bidirectional_math isa Union{Real,Symbol}) || (SymbolicUtils.operation(bidirectional_math) != +)
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

""" Extract paramap from Model """
function get_paramap(model)
    paramap = Pair{Num,Float64}[]
    for (k, v) in model.parameters
        push!(paramap, Pair(create_param(k), v[1])) # [1] index is dropping unit (i think)
    end
    for (k, v) in model.compartments
        if !isnothing(v.size)
            push!(paramap, Pair(create_param(k), v.size))
        end
    end
    paramap
end

get_u0map(model) = [create_var(k, Catalyst.DEFAULT_IV) => v for (k, v) in SBML.initial_amounts(model, convert_concentrations = true)]
ModelingToolkit.defaults(model::SBML.Model) = ModelingToolkit._merge(get_u0map(model), get_paramap(model))

""" Get rate constant of mass action kineticLaws """
function getmassaction(kl::Num, reactants::Union{Vector{Num},Nothing}, stoich::Union{Vector{<:Real},Nothing})
    function check_args(x::SymbolicUtils.Symbolic{Real})
        for arg in SymbolicUtils.arguments(x)
            if isnan(check_args(arg)) || isequal(arg, Catalyst.DEFAULT_IV)
                return NaN
            end
        end
        return 0
    end
    check_args(x::Term{Real,Nothing}) = NaN  # Variable leaf node
    check_args(x::Sym{Real,Base.ImmutableDict{DataType,Any}}) = 0  # Parameter leaf node
    check_args(x::Real) = 0  # Real leaf node
    check_args(x) = throw(ErrorException("Cannot handle $(typeof(x)) types."))  # Unknow leaf node
    if isnothing(reactants) && isnothing(stoich)
        rate_const = kl
    elseif isnothing(reactants) | isnothing(stoich)
        throw(ErrorException("`reactants` and `stoich` are incosistent: `reactants` are $(reactants) and `stoich` is $(stoich)."))
    else
        rate_const = SymbolicUtils.simplify_fractions(kl / *((.^(reactants, stoich))...))
    end
    isnan(check_args(rate_const.val)) ? NaN : rate_const
end

function create_var(x)
    sym = Symbol(x)
    Symbolics.unwrap(first(@variables $sym))
end
function create_var(x, iv)
    sym = Symbol(x)
    Symbolics.unwrap(first(@variables $sym(iv)))
end
function create_param(x)
    sym = Symbol(x)
    Symbolics.unwrap(first(@parameters $sym))
end

function get_events(model, rs)
    evs = model.events
    mtk_evs = Pair{Vector{Equation},Vector{Equation}}[]
    for (_, e) in evs
        trig = [~(convert(Num, e.trigger).val.arguments...)] # need to convert S1 -> S1(t)
        mtk_evas = Equation[]
        for eva in e.event_assignments
            var = Symbol(eva.variable)
            # @info eva.math.val
            pair = ModelingToolkit.getvar(rs, var) ~ convert(Num, eva.math)
            push!(mtk_evas, pair)
        end
        push!(mtk_evs, trig => mtk_evas)
    end
    mtk_evs
end

function get_rules(model)
    subsdict = _get_substitutions(model)
    # goes into `observed`
    obseqs = Equation[]
    # these two go into `constraints` field of ReactionSystem
    algeqs = Equation[]
    raterules = Equation[]

    rules = model.rules
    for r in rules
        if r isa SBML.AlgebraicRule
            push!(algeqs, 0 ~ convert(Num, r.math))
        elseif r isa SBML.AssignmentRule
            @info r
            push!(obseqs, assignmentrule_to_obseq(model, r))
        elseif r isa SBML.RateRule
            push!(raterules, raterule_to_diffeq(model, r))
        else
            error()
        end
    end
    algeqs, obseqs, raterules = map(x -> substitute(x, subsdict), (algeqs, obseqs, raterules))
    algeqs, obseqs, raterules
end

function assignmentrule_to_obseq(model, rule)
    @info rule
    if haskey(model.species, rule.id)
        sym = Symbol(rule.id)
        var = Symbolics.unwrap(first(@variables $sym(Catalyst.DEFAULT_IV)))
        assignment = Num(convert(Num, rule.math, convert_time = (x::SBML.MathTime) -> Catalyst.DEFAULT_IV))
        return var ~ assignment
    elseif haskey(model.compartments, rule.id)
        error("not handling setting compartment volumes rn")
    elseif haskey(model.parameters, rule.id)
        error("not handling setting non-constant parameters rn")
    else
        error()
    end
end

function raterule_to_diffeq(model, rule)
    # the rule.id can be a species, speciesRef, compartment, or param. currently not doingthis
    # for now im just doing species
    D = Differential(Catalyst.DEFAULT_IV)
    if haskey(model.species, rule.id)
        sym = Symbol(rule.id)
        var = Symbolics.unwrap(first(@variables $sym(Catalyst.DEFAULT_IV)))
        assignment = Num(convert(Num, rule.math, convert_time = (x::SBML.MathTime) -> Catalyst.DEFAULT_IV))
        return D(var) ~ assignment
    elseif haskey(model.comnpartments, rule.id)
        error("not handling setting compartment volumes rn")
    elseif haskey(model.parameters, rule.id)
        error("not handling setting non-constant parameters rn")
    else
        error()
    end
end

