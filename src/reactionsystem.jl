# Conversion to symbolics
symbolicsRateOf(x) = Symbolics.Differential(convert(Num, MathTime("t")))(x)

symbolics_mapping = Dict(SBML.default_function_mapping..., "rateOf" => symbolicsRateOf)

map_symbolics_time_ident(x) = begin
    sym = Symbol(x.id)
    first(@variables $sym)
end

interpret_as_num(x::SBML.Math) = SBML.interpret_math(
    x;
    map_apply = (x::SBML.MathApply, interpret::Function) ->
        Num(symbolics_mapping[x.fn](interpret.(x.args)...)),
    map_const = (x::SBML.MathConst) -> Num(SBML.default_constants[x.id]),
    map_ident = map_symbolics_time_ident,
    map_lambda = (_, _) ->
        throw(ErrorException("Symbolics.jl does not support lambda functions")),
    map_time = (x::SBML.MathTime) -> Catalyst.DEFAULT_IV,
    map_value = (x::SBML.MathVal) -> Num(x.val),
)

""" ReactionSystem constructor """
function Catalyst.ReactionSystem(model::SBML.Model; kwargs...)  # Todo: requires unique parameters (i.e. SBML must have been imported with localParameter promotion in libSBML)
    # length(model.events) > 0 ? error("Model contains events. Please import with `ODESystem(model)`") : nothing  @Anand: how to suppress this when called from ODESystem
    rxs = mtk_reactions(model)
    u0map = get_u0map(model)
    parammap = get_paramap(model)
    defs = ModelingToolkit._merge(Dict(u0map), Dict(parammap))

    algrules, obsrules, raterules = get_rules(model)
    # constant_eqs = fix_constants_at_init(model)  # Todo PL: take this out
    unassigned_pars = fix_unassigned_nonconstant_par_to_init(model)  # @Sam: if I do not do that, there are too few equations. Do you have better suggestions?
    for o in obsrules
        defs[o.lhs] = substitute(o.rhs, defs)
    end
    constraints_sys = ODESystem(vcat(algrules, raterules, obsrules, unassigned_pars),
                                Catalyst.DEFAULT_IV; name = gensym(:CONSTRAINTS))

    # Hacky way to keep zero_ode species with contant=boundaryCondition=false in system
    rs = ReactionSystem(rxs, Catalyst.DEFAULT_IV, first.(u0map), first.(parammap); defaults = defs, name = gensym(:SBML),
        constraints = constraints_sys, kwargs...)
    odessys = convert(ODESystem, rs; include_zero_odes=true, combinatoric_ratelaws=false)
    unchanged_eqs = fix_zero_odes_to_init(model, equations(odessys))  # Todo PL: Use input=true instead
    constraints_sys = ODESystem(vcat(algrules, raterules, obsrules, unassigned_pars, unchanged_eqs),
                                Catalyst.DEFAULT_IV; name = gensym(:CONSTRAINTS))
    ReactionSystem(rxs, Catalyst.DEFAULT_IV, first.(u0map), first.(parammap); defaults = defs, name = gensym(:SBML),
        constraints = constraints_sys, kwargs...)
end

function fix_zero_odes_to_init(model, equations)
    fix_to_init = Equation[]
    inits = Dict(SBML.initial_amounts(model, convert_concentrations = true))
    D = Differential(Catalyst.DEFAULT_IV)
    assigned_species = [r.id for r in values(model.rules) if r isa Union{SBML.AssignmentRule, SBML.RateRule}]
    algebraic_species = get_algebraic_species(model)
    for (k, v) in model.species
        var = create_var(k, Catalyst.DEFAULT_IV; constant=v.constant, boundary_condition=v.boundary_condition)
        if v.constant == false && !(k in vcat(assigned_species, algebraic_species))
            if v.boundary_condition == false
                eq = D(var) ~ 0
                if eq in equations
                    push!(fix_to_init, var ~ inits[k])
                end
            elseif v.boundary_condition == true
                push!(fix_to_init, var ~ inits[k])
            end
        end
    end
    fix_to_init
end

function fix_unassigned_nonconstant_par_to_init(model)
    fix_to_init = Equation[]
    assigned_pars = [r.id for r in values(model.rules) if r isa Union{SBML.AssignmentRule, SBML.RateRule}]
    algebraic_pars = get_algebraic_species(model; kind="parameter")
    for (k, v) in model.parameters
        var = create_var(k, Catalyst.DEFAULT_IV)
        if v.constant == false && !(k in vcat(assigned_pars, algebraic_pars))
            push!(fix_to_init, var ~ v.value)
        end
    end
    fix_to_init
end

function get_algebraic_species(model; kind="species")
    algebraic_rules = [r.math for r in model.rules if r isa SBML.AlgebraicRule]
    algebraic_species = String[]
    for math in algebraic_rules
        parse_math!(math, algebraic_species, model; kind=kind)
    end
    algebraic_species
end

function parse_math!(math, algebraic_species, model; kind="species")
    k = kind == "species" ? model.species : kind == "parameter" ? model.parameters : kind == "compartment" ? model.compartments : throw(ArgumentError("kind must be either 'species', 'parameter' or 'compartment'"))
    if math isa SBML.MathIdent && math.id in keys(k)
        push!(algebraic_species, math.id)
    elseif math isa Union{SBML.MathApply, SBML.MathLambda}
        [parse_math!(arg, algebraic_species, model; kind=kind) for arg in math.args]
    end
    nothing
end

# function fix_constants_at_init(model)
#     fixed_at_init = Equation[]
#     inits = Dict(SBML.initial_amounts(model, convert_concentrations = true))
#     for (k, v) in model.species
#         if v.constant==true
#             push!(fixed_at_init, create_var(k, Catalyst.DEFAULT_IV; constant=v.constant, boundary_condition=v.boundary_condition) ~ inits[k])
#         end
#     end
#     fixed_at_init
# end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(model::SBML.Model; include_zero_odes = false, kwargs...)
    rs = ReactionSystem(model; kwargs...)
    convert(ODESystem, rs; include_zero_odes = include_zero_odes,
            continuous_events=get_events(model, rs), combinatoric_ratelaws=false)
end

""" Check if conversion to ReactionSystem is possible """
function checksupport(filename::String)
    not_implemented = ["listOfConstraints", "</delay>", "<priority>"]
    sbml = open(filename) do file
        read(file, String)
    end
    for item in not_implemented
        occursin(item, sbml) && throw(ErrorException("SBML models with $item are not yet implemented."))
    end
    occursin("<sbml xmlns:fbc=", sbml) && throw(ErrorException("This model was designed for constrained-based optimisation. Please use COBREXA.jl instead of SBMLToolkit."))
end

""" Get dictonary to change types in kineticLaw """
function _get_substitutions(model)
    subsdict = Dict()
    for (k, v) in model.species
        push!(subsdict, Pair(create_var(k), create_var(k, Catalyst.DEFAULT_IV; constant=v.constant, boundary_condition=v.boundary_condition)))
    end
    for (k, v) in model.parameters
        if v.constant !== nothing && v.constant
            push!(subsdict, Pair(create_var(k), create_param(k)))
        else
            push!(subsdict, Pair(create_var(k), create_var(k, Catalyst.DEFAULT_IV)))
        end
    end
    for (k, v) in model.compartments
        if v.constant
            push!(subsdict, Pair(create_var(k), create_param(k)))
        else
            push!(subsdict, Pair(create_var(k), create_var(k, Catalyst.DEFAULT_IV)))
        end
    end
    subsdict
end

function stoich_convert_to_ints(xs)
    (xs !== nothing && all(isinteger(x) for x in xs)) ? Int.(xs) : xs
end

""" Convert SBML.Reaction to MTK.Reaction """
function mtk_reactions(model::SBML.Model)
    subsdict = _get_substitutions(model)
    rxs = []
    for reaction in values(model.reactions)
        extensive_math = SBML.extensive_kinetic_math(
            model, reaction.kinetic_math)
        symbolic_math = interpret_as_num(extensive_math)

        rstoich = reaction.reactants
        pstoich = reaction.products
        if reaction.reversible
            symbolic_math = getunidirectionalcomponents(symbolic_math)
            kl_fw, kl_rv = [substitute(x, subsdict) for x in symbolic_math]
            reactants, products, rstoichvals, pstoichvals = getreagents(rstoich, pstoich, model)
            isnothing(reactants) && isnothing(products) && continue
            rstoichvals = stoich_convert_to_ints(rstoichvals)
            pstoichvals = stoich_convert_to_ints(pstoichvals)
            kl_fw, our = use_rate(kl_fw, reactants, rstoichvals)
            push!(rxs, Catalyst.Reaction(kl_fw, reactants, products, rstoichvals, pstoichvals; only_use_rate = our))

            reagents = getreagents(rstoich, pstoich, model; rev = true)
            reactants_rev, products_rev, rstoichvals_rev, pstoichvals_rev = reagents
            rstoichvals_rev = stoich_convert_to_ints(rstoichvals_rev)
            pstoichvals_rev = stoich_convert_to_ints(pstoichvals_rev)
            kl_rv, our = use_rate(kl_rv, reactants_rev, rstoichvals_rev)
            push!(rxs, Catalyst.Reaction(kl_rv, reactants_rev, products_rev, rstoichvals_rev, pstoichvals_rev; only_use_rate = our))
        else
            kl = substitute(symbolic_math, subsdict)
            reactants, products, rstoichvals, pstoichvals = getreagents(rstoich, pstoich, model)
            isnothing(reactants) && isnothing(products) && continue
            rstoichvals = stoich_convert_to_ints(rstoichvals)
            pstoichvals = stoich_convert_to_ints(pstoichvals)
            kl, our = use_rate(kl, reactants, rstoichvals)
            push!(rxs, Catalyst.Reaction(kl, reactants, products, rstoichvals, pstoichvals; only_use_rate = our))
        end
    end
    rxs
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
function get_paramap(model)  # Todo PL: merge with get_u0map.
    assigned_pars = [r.id for r in values(model.rules) if r isa Union{SBML.AssignmentRule, SBML.RateRule}]
    algebraic_pars = get_algebraic_species(model; kind="parameter")
    paramap = Pair{Num,Float64}[]
    for (k, v) in model.parameters
        if v.constant || isnothing(v.value)  # !(k in vcat(assigned_pars, algebraic_pars))
            push!(paramap, Pair(create_param(k), v.value))
        end
    end
    algebraic_comps = get_algebraic_species(model; kind="compartment")
    for (k, v) in model.compartments
        # if !isnothing(v.size) && (v.constant || !(k in vcat(assigned_pars, algebraic_comps)))
        if isnothing(v.size) || v.constant
            push!(paramap, Pair(create_param(k), v.size))
        end
    end
    paramap
end

""" Extract paramap from Model """
function get_u0map(model)
    u0s = Pair[]
    inits = Dict(SBML.initial_amounts(model, convert_concentrations = true))

    for (k, v) in model.species
        p = create_var(k, Catalyst.DEFAULT_IV; constant=v.constant, boundary_condition=v.boundary_condition) => inits[k]
        push!(u0s, p)
    end

    for (k, v) in model.compartments
        if !isnothing(v.size) && !v.constant
            push!(u0s, Pair(create_var(k, Catalyst.DEFAULT_IV), v.size))
        end
    end

    for (k, v) in model.parameters
        if !isnothing(v.value) && !v.constant
            push!(u0s, Pair(create_var(k, Catalyst.DEFAULT_IV), v.value))
        end
    end
    u0s
end

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
function create_var(x, iv; constant=false, boundary_condition=false)
    sym = Symbol(x)
    # Symbolics.unwrap(first(@variables $sym(iv)))
    Symbolics.unwrap(first(@variables $sym(iv) [isconstant=constant, isbc=boundary_condition]))
end
function create_param(x)
    sym = Symbol(x)
    Symbolics.unwrap(first(@parameters $sym))
end
# function create_param(x, input)
#     sym = Symbol(x)
#     Symbolics.unwrap(first(@parameters $sym [input=input]))
# end

function get_rules(model)
    subsdict = _get_substitutions(model)
    # these three go into `constraints` field of ReactionSystem
    obseqs = Equation[]
    algeqs = Equation[]
    raterules = Equation[]

    rules = model.rules
    for r in rules
        if r isa SBML.AlgebraicRule
            push!(algeqs, 0 ~ interpret_as_num(r.math))
        elseif r isa SBML.AssignmentRule
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

function rule_to_var_and_eq(model, rule; volume_correction=nothing)
    sym = Symbol(rule.id)
    var = Symbolics.unwrap(first(@variables $sym(Catalyst.DEFAULT_IV)))
    math = SBML.extensive_kinetic_math(model, rule.math)
    if !isnothing(volume_correction)
        math = SBML.MathApply("*", [SBML.MathIdent(volume_correction), math])
    end
    assignment = interpret_as_num(math)
    var, assignment
end

handle_empty_compartment_size = (id::String) -> throw(
    DomainError(
        "Non-substance-only-unit reference to species `$id' in an unsized compartment."
    ))

function volume_correction(model, s_id)
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

function assignmentrule_to_obseq(model, rule)
    if haskey(model.species, rule.id)
        vc = volume_correction(model, rule.id)
        var, assignment = rule_to_var_and_eq(model, rule; volume_correction=vc)
        return var ~ assignment
    elseif haskey(model.compartments, rule.id)
        var, assignment = rule_to_var_and_eq(model, rule)
        return var ~ assignment
    elseif haskey(model.parameters, rule.id)
        var, assignment = rule_to_var_and_eq(model, rule)
        return var ~ assignment
    else
        error("invalid rule: $rule")
    end
end

function raterule_to_diffeq(model, rule)
    D = Differential(Catalyst.DEFAULT_IV)
    if haskey(model.species, rule.id)
        vc = volume_correction(model, rule.id)
        var, assignment = rule_to_var_and_eq(model, rule; volume_correction=vc)
        return D(var) ~ assignment
    elseif haskey(model.compartments, rule.id)
        var, assignment = rule_to_var_and_eq(model, rule)
        return D(var) ~ assignment
    elseif haskey(model.parameters, rule.id)
        var, assignment = rule_to_var_and_eq(model, rule)
        return D(var) ~ assignment
    else
        error()
    end
end

"""
    Creates ContinuousVectorCallbacks


Note that one limitation of Event support is that ReactionSystems do not have a field for it yet.
So in order for the system to have events, you must call `ODESystem(m::SBML.Model)` rather than `convert(ODESystem, ReactionSystem(m::SBML.Model))`
"""
function get_events(model, rs)
    subsdict = _get_substitutions(model)
    evs = model.events
    mtk_evs = Pair{Vector{Equation},Vector{Equation}}[]
    for (_, e) in evs
        trigger = SBML.extensive_kinetic_math(model, e.trigger)
        args = Symbolics.unwrap(interpret_as_num(trigger))
        lhs, rhs = map(x -> substitute(x, subsdict), args.arguments)
        trig = [lhs ~ rhs]
        mtk_evas = Equation[]
        for eva in e.event_assignments
            if haskey(model.species, eva.variable)
                vc = volume_correction(model, eva.variable)
                if !isnothing(vc)
                    math = SBML.MathApply("*", [SBML.MathIdent(vc), eva.math])
                end
            end
            var = Symbol(eva.variable)
            pair = ModelingToolkit.getvar(rs, var) ~ Symbolics.unwrap(interpret_as_num(math))
            push!(mtk_evas, pair)
        end
        push!(mtk_evs, trig => mtk_evas)
    end
    mtk_evs
end
