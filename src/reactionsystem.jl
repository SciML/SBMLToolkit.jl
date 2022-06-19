""" ReactionSystem constructor """
function Catalyst.ReactionSystem(model::SBML.Model; kwargs...)  # Todo: requires unique parameters (i.e. SBML must have been imported with localParameter promotion in libSBML)
    # length(model.events) > 0 ? error("Model contains events. Please import with `ODESystem(model)`") : nothing  @Anand: how to suppress this when called from ODESystem
    rxs = get_reactions(model)
    u0map, parammap = get_mappings(model)
    defs = ModelingToolkit._merge(Dict(u0map), Dict(parammap))

    algrules, obsrules, raterules = get_rules(model)
    for o in obsrules
        defs[o.lhs] = substitute(o.rhs, defs)
    end
    constraints_sys = ODESystem(vcat(algrules, raterules, obsrules),
                                IV; name = gensym(:CONSTRAINTS))
    ReactionSystem(rxs, IV, first.(u0map), first.(parammap);
                   defaults = defs, name = gensym(:SBML),
                   constraints = constraints_sys, kwargs...)
end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(model::SBML.Model; include_zero_odes = true, kwargs...)
    rs = ReactionSystem(model; kwargs...)
    convert(ODESystem, rs; include_zero_odes = include_zero_odes,
            continuous_events=get_events(model, rs), combinatoric_ratelaws=false)
end

function get_mappings(model::SBML.Model)
    inits = Dict(SBML.initial_amounts(model, convert_concentrations = true))
    u0map = Pair[]
    parammap = Pair[]
    for (k, v) in model.species
        if v.constant == true
            var = create_param(k; isconstantspecies=true)
            push!(parammap, var => inits[k])
        else
            var = create_var(k, IV;
                            isbcspecies=has_rule_type(k, model, SBML.RateRule) ||
                                has_rule_type(k, model, SBML.AssignmentRule) ||
                                (has_rule_type(k, model, SBML.AlgebraicRule) &&
                                (
                                    all([netstoich(k, r) == 0 for r in values(model.reactions)]) ||
                                    v.boundary_condition == true)
                                )
                            )  # To remove species that are otherwise defined
            push!(u0map, var => inits[k])
        end
    end
    for (k, v) in model.parameters
        if v.constant == false && SBML.seemsdefined(k, model)
            var = create_var(k, IV; isbcspecies=true)
            push!(u0map, var => v.value)
        else
            var = create_param(k)
            push!(parammap, var => v.value)
        end
    end
    for (k, v) in model.compartments
        if v.constant == false && SBML.seemsdefined(k, model)
            var = create_var(k, IV; isbcspecies=true)
            push!(u0map, var => v.size)
        else
            var = create_param(k)
            push!(parammap, var => v.size)
        end
    end
    u0map, parammap
end

function netstoich(id, reaction)
    netstoich = 0
    netstoich -= get(reaction.reactants, id, 0)
    netstoich += get(reaction.products, id, 0)
end

""" Check if conversion of file to ReactionSystem is possible """
function checksupport_file(filename::String)
    string = open(filename) do file
        read(file, String)
    end
    checksupport_string(string)
end

""" Check if conversion of xml-string to ReactionSystem is possible """
function checksupport_string(xml::String)
    not_implemented = ["listOfConstraints", "</delay>",
                       "<priority>", "spatialDimensions=\"0\"",
                       "factorial", "00387"]  # Case 00387 requires event directionality
    for item in not_implemented
        occursin(item, xml) && throw(ErrorException("SBML models with $item are not yet implemented."))
    end
    occursin("<sbml xmlns:fbc=", xml) && throw(ErrorException("This model was designed for constrained-based optimisation. Please use COBREXA.jl instead of SBMLToolkit."))
    !(occursin("<reaction", xml) || occursin("rateRule", xml)) && throw(ErrorException("Models that contain neither reactions or rateRules will fail in simulation."))
    true
end

### reactions.jl ###
""" Convert SBML.Reaction to MTK.Reaction """
function get_reactions(model::SBML.Model)
    subsdict = get_substitutions(model)  # Todo: replace with SUBSDICT
    rxs = Reaction[]
    for reaction in values(model.reactions)
        extensive_math = SBML.extensive_kinetic_math(
            model, reaction.kinetic_math)
        symbolic_math = interpret_as_num(extensive_math)
        rstoich = reaction.reactants
        pstoich = reaction.products
        if reaction.reversible
            symbolic_math = get_unidirectional_components(symbolic_math)
            kl_fw, kl_rv = [substitute(x, subsdict) for x in symbolic_math]
            add_reaction!(rxs, kl_fw, rstoich, pstoich, model)
            add_reaction!(rxs, kl_rv, pstoich, rstoich, model)
        else
            kl = substitute(symbolic_math, subsdict)  # Todo: use SUBSDICT
            add_reaction!(rxs, kl, rstoich, pstoich, model)
        end
    end
    rxs
end

""" Infer forward and reverse components of bidirectional kineticLaw """
function get_unidirectional_components(bidirectional_math)
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

function add_reaction!(rxs::Vector{Reaction},
                       kl::Num,
                       rstoich::Dict{String,Float64}, pstoich::Dict{String,Float64},
                       model::SBML.Model)  
    reactants, products, rstoichvals, pstoichvals = get_reagents(rstoich, pstoich, model)
    isnothing(reactants) && isnothing(products) && return
    rstoichvals = stoich_convert_to_ints(rstoichvals)
    pstoichvals = stoich_convert_to_ints(pstoichvals)
    kl, our = use_rate(kl, reactants, rstoichvals)
    push!(rxs, Catalyst.Reaction(kl, reactants, products, rstoichvals, pstoichvals; only_use_rate = our))
end

function stoich_convert_to_ints(xs)
    (xs !== nothing && all(isinteger(x) for x in xs)) ? Int.(xs) : xs
end

""" Get reagents """
function get_reagents(rstoichdict::Dict{String,<:Real},
                      pstoichdict::Dict{String,<:Real},
                      model::SBML.Model)
    reactants = Num[]
    products = Num[]
    rstoich = Float64[]
    pstoich = Float64[]

    for (k, v) in rstoichdict
        iszero(v) && @error("Stoichiometry of $k must be non-zero")
        push!(reactants, create_var(k, IV))
        push!(rstoich, v)
        if model.species[k].boundary_condition == true
            push!(products, create_var(k, IV))
            push!(pstoich, v)
        end
    end
    for (k, v) in pstoichdict
        iszero(v) && @error("Stoichiometry of $k must be non-zero")
        if model.species[k].boundary_condition != true
            push!(products, create_var(k, IV))
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

""" Get kineticLaw for use in MTK.Reaction """
function use_rate(kl::Num, react::Union{Vector{Num},Nothing}, stoich::Union{Vector{<:Real},Nothing})
    rate_const = get_massaction(kl, react, stoich)
    if !isnan(rate_const)
        kl = rate_const
        our = false
    else
        our = true
    end
    return (kl, our)
end

""" Get rate constant of mass action kineticLaws """
function get_massaction(kl::Num, reactants::Union{Vector{Num},Nothing}, stoich::Union{Vector{<:Real},Nothing})
    function check_args(x::SymbolicUtils.Symbolic{Real})
        for arg in SymbolicUtils.arguments(x)
            if isnan(check_args(arg)) || isequal(arg, Catalyst.DEFAULT_IV)
                return NaN
            end
        end
        return 0
    end
    check_args(_::Term{Real,Nothing}) = NaN  # Variable leaf node
    check_args(_::Sym{Real,Base.ImmutableDict{DataType,Any}}) = 0  # Parameter leaf node
    check_args(_::Real) = 0  # Real leaf node
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

### rules.jl ###
function get_rules(model)
    subsdict = get_substitutions(model)  # Todo: use SUBSDICT
    obseqs = Equation[]
    algeqs = Equation[]
    raterules = Equation[]
    for r in model.rules
        if r isa SBML.AlgebraicRule
            push!(algeqs, 0 ~ interpret_as_num(r.math))
        elseif r isa SBML.AssignmentRule
            var, ass = get_var_and_assignment(model, r)
            push!(obseqs, var ~ ass)
        elseif r isa SBML.RateRule
            var, ass = get_var_and_assignment(model, r)
            push!(raterules, D(var) ~ ass)
        else
            error("Rule must be of type SBML.AlgebraicRule, SBML.AssignmentRule, or SBML.RateRule.")
        end
    end
    algeqs, obseqs, raterules = map(x -> substitute(x, subsdict), (algeqs, obseqs, raterules))
    algeqs, obseqs, raterules
end

function get_var_and_assignment(model, rule)
    if !haskey(merge(model.species, model.compartments, model.parameters), rule.id)
        error("Cannot find target for rule with ID `$rule.id`")
    end
    var = create_var(rule.id, IV)
    math = SBML.extensive_kinetic_math(model, rule.math)
    vc = get_volume_correction(model, rule.id)
    if !isnothing(vc)
        math = SBML.MathApply("*", [SBML.MathIdent(vc), math])
    end
    assignment = interpret_as_num(math)
    var, assignment
end

function get_volume_correction(model, s_id)
    haskey(model.species, s_id) || return nothing
    sp = model.species[s_id]
    comp = model.compartments[sp.compartment]
    sp.only_substance_units == true && return nothing  
    isnothing(comp.size) && !SBML.seemsdefined(sp.compartment, model) &&
        comp.spatial_dimensions != 0 &&  # remove this line when SBML test suite is fixed
        throw(
            DomainError(
                sp.compartment,
                "compartment size is insufficiently defined",
            ),
        )  
    sp.compartment
end

### utils.jl ###
# Conversion to symbolics
symbolicsRateOf(x) = Differential(t)(x)

const IV = Catalyst.DEFAULT_IV
const symbolics_mapping = Dict(SBML.default_function_mapping..., "rateOf" => symbolicsRateOf)
const D = Differential(IV)
# const SUBSDICT = get_substitutions(model)

map_symbolics_ident(x) = begin
    sym = Symbol(x.id)
    first(@variables $sym)
end

interpret_as_num(x::SBML.Math) = SBML.interpret_math(
    x;
    map_apply = (x::SBML.MathApply, interpret::Function) ->
        Num(symbolics_mapping[x.fn](interpret.(x.args)...)),
    map_const = (x::SBML.MathConst) -> Num(SBML.default_constants[x.id]),
    map_ident = map_symbolics_ident,
    map_lambda = (_, _) ->
        throw(ErrorException("Symbolics.jl does not support lambda functions")),
    map_time = (x::SBML.MathTime) -> IV,
    map_value = (x::SBML.MathVal) -> Num(x.val),
)

""" Get dictonary to change types in kineticLaw """
function get_substitutions(model)
    u0map, parammap = get_mappings(model)
    subsdict = Dict()
    for item in first.(u0map)
        k = create_var(string(item.f.name))
        subsdict[k] = item
    end
    for item in first.(parammap)
        k = create_var(string(item.name))
        subsdict[k] = item
    end
    subsdict
end

function create_var(x; isbcspecies=false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@variables $sym [isbcspecies=isbcspecies]))
end
function create_var(x, iv; isbcspecies=false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@variables $sym(iv) [isbcspecies=isbcspecies]))
end
function create_param(x; isconstantspecies=false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@parameters $sym [isconstantspecies=isconstantspecies]))
end

function has_rule_type(id::String, m::SBML.Model, T::Type{<:SBML.Rule})
    T == SBML.AlgebraicRule && return any(SBML.isfreein(id, r.math) for r in m.rules if r isa SBML.AlgebraicRule)
    any(r.id == id for r in m.rules if r isa T)
end

### events.jl ###
"""
    Creates ContinuousVectorCallbacks


Note that one limitation of Event support is that ReactionSystems do not have a field for it yet.
So in order for the system to have events, you must call `ODESystem(m::SBML.Model)` rather than `convert(ODESystem, ReactionSystem(m::SBML.Model))`
"""
function get_events(model, rs)  # Todo: implement up or downpass and parameters
    subsdict = get_substitutions(model)  # Todo: use SUBSDICT
    evs = model.events
    mtk_evs = Pair{Vector{Equation},Vector{Equation}}[]
    for (_, e) in evs
        trigger = SBML.extensive_kinetic_math(model, e.trigger)
        trigger = Symbolics.unwrap(interpret_as_num(trigger))
        lhs, rhs = map(x -> substitute(x, subsdict), trigger.arguments)
        trig = [lhs ~ rhs]
        mtk_evas = Equation[]
        for eva in e.event_assignments
            math = eva.math
            if haskey(model.species, eva.variable)
                vc = get_volume_correction(model, eva.variable)
                if !isnothing(vc)
                    math = SBML.MathApply("*", [SBML.MathIdent(vc), math])
                end
            end
            var = create_var(eva.variable, IV)
            math = substitute(Symbolics.unwrap(interpret_as_num(math)),
                              subsdict)
            effect = var ~ math
            push!(mtk_evas, effect)
        end
        push!(mtk_evs, trig => mtk_evas)
    end
    mtk_evs
end
