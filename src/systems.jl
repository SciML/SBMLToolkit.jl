"DefaultImporter to use in conjunction with `readSBML`"
struct DefaultImporter end

"ReactionSystemImporter to use in conjunction with `readSBML`"
struct ReactionSystemImporter end

"ODESystemImporter to use in conjunction with `readSBML`"
struct ODESystemImporter end

"""
    readSBML(sbmlfile::String, ::DefaultImporter)

Create a `SBML.Model` from an SBML file, using the default import settings for use as Catalyst and ModelingToolkit types.

See also [`Model`](@ref) and [`DefaultImporter`](@ref).
"""
function SBML.readSBML(sbmlfile::String, ::DefaultImporter)  # Returns an SBML.Model
    SBMLToolkit.checksupport_file(sbmlfile)
    return readSBML(sbmlfile, importdefaults)
end

"""
    readSBML(sbmlfile::String, ::ReactionSystemImporter)

Create a `Catalyst.ReactionSystem` from an SBML file, using the default import settings.

See also [`Model`](@ref) and [`ReactionSystemImporter`](@ref).
"""
function SBML.readSBML(sbmlfile::String, ::ReactionSystemImporter; kwargs...)  # Returns a Catalyst.ReactionSystem
    return ReactionSystem(readSBML(sbmlfile::String, DefaultImporter()), kwargs...)
end

"""
    readSBML(sbmlfile::String, ::ODESystemImporter)

Create a `ModelingToolkit.ODESystem` from an SBML file, using the default import settings.

See also [`Model`](@ref) and [`ODESystemImporter`](@ref).
"""
function SBML.readSBML(
        sbmlfile::String, ::ODESystemImporter;
        include_zero_odes::Bool = true, kwargs...
    )  # Returns an MTK.ODESystem
    odesys = convert(
        ODESystem, readSBML(sbmlfile, ReactionSystemImporter(), kwargs...),
        include_zero_odes = include_zero_odes
    )
    return complete(odesys)
end

"""
    ReactionSystem(model::SBML.Model; kwargs...)

Create a `ReactionSystem` from an `SBML.Model`.

See also [`ODESystem`](@ref).
"""
function Catalyst.ReactionSystem(model::SBML.Model; kwargs...)  # Todo: requires unique parameters (i.e. SBML must have been imported with localParameter promotion in libSBML)
    # length(model.events) > 0 ? error("Model contains events. Please import with `ODESystem(model)`") : nothing  @Anand: how to suppress this when called from ODESystem
    rxs = get_reactions(model)
    u0map, parammap, initial_assignment_map = get_mappings(model)
    defs = Dict{Num, Any}()
    for (k, v) in vcat(u0map, parammap, initial_assignment_map)  # initial_assignments override u0map and parammap
        defs[k] = v
    end
    # defs = ModelingToolkit._merge(Dict(u0map), Dict(parammap))
    algrules, obsrules, raterules = get_rules(model)
    obsrules_rearranged = Equation[]
    for o in obsrules
        rhs = o.rhs
        for r in raterules
            if isequal(rhs, r.lhs)
                rhs = r.rhs
            end
        end
        defs[o.lhs] = ModelingToolkit.fixpoint_sub(rhs, defs)
        #  ModelingToolkit._merge(defs,
        #                         Dict(Catalyst.DEFAULT_IV.val => 0)))
        push!(obsrules_rearranged, 0 ~ rhs - o.lhs)
    end
    raterules_subs = []
    for o in raterules
        rhs = o.rhs
        for r in raterules
            if isequal(rhs, r.lhs)
                rhs = r.rhs
            end
        end
        defs[o.lhs] = ModelingToolkit.fixpoint_sub(rhs, defs)
        #  ModelingToolkit._merge(defs,
        #                         Dict(Catalyst.DEFAULT_IV.val => 0)))
        push!(raterules_subs, rhs ~ o.lhs)
    end
    zero_rates = []
    for (k, v) in merge(model.parameters, model.compartments)
        if is_event_assignment(k, model)
            if v.constant
                ErrorException("$k appears in event assignment but should be constant.")
            end
            push!(zero_rates, D(create_var(k, IV; isbcspecies = true)) ~ 0)
        end
    end
    if haskey(kwargs, :defaults)
        defs = ModelingToolkit._merge(defs, kwargs[:defaults])
        kwargs = filter(x -> !isequal(first(x), :defaults), kwargs)
    end
    unknown_keys = first.(u0map)
    param_keys = first.(parammap)
    rs = _build_reaction_system(
        [rxs..., algrules..., raterules_subs..., obsrules_rearranged..., zero_rates...],
        IV, unknown_keys, param_keys, defs,
        get_events(model);
        kwargs...
    )
    return complete(rs)  # Todo: maybe add a `complete=True` kwarg
end

# Catalyst v16 replaced the single `defaults` kwarg of `ReactionSystem` with
# separate `bindings` (parameter values) and `initial_conditions` (species/unknown
# values). Older Catalyst (v14, v15) still uses `defaults`. This shim splits the
# merged dictionary and dispatches on the installed Catalyst version so the same
# SBMLToolkit code base supports both APIs.
@static if pkgversion(Catalyst) >= v"16"
    function _build_reaction_system(eqs, iv, unknowns, ps, defs, cevs; kwargs...)
        unknown_set = Set(SymbolicUtils.unwrap(u) for u in unknowns)
        param_set = Set(SymbolicUtils.unwrap(p) for p in ps)
        bindings = Dict{Any, Any}()
        initial_conditions = Dict{Any, Any}()
        for (k, v) in defs
            if SymbolicUtils.unwrap(k) in unknown_set
                initial_conditions[k] = v
            elseif _references_non_parameter(v, param_set)
                # SBML `initialAssignment`s on parameters may reference species
                # (e.g. `parameter S3 := k1*S2`); those cannot live in `bindings`
                # under Catalyst v16 (`check_bindings` rejects non-parameter symbols)
                # but are accepted as `initial_conditions`, which is evaluated at
                # t=0 with the same SBML semantics.
                initial_conditions[k] = v
            else
                bindings[k] = v
            end
        end
        return ReactionSystem(
            eqs, iv, unknowns, ps;
            bindings = bindings, initial_conditions = initial_conditions,
            name = gensym(:SBML), continuous_events = cevs,
            combinatoric_ratelaws = false, kwargs...
        )
    end

    function _references_non_parameter(v, param_set)
        # Note: Symbolics.Num <: Real <: Number, so a `v isa Number` early-return
        # would swallow symbolic values. Use `get_variables` directly — it returns
        # an empty iterator for plain numerics.
        for var in Symbolics.get_variables(v)
            SymbolicUtils.unwrap(var) in param_set || return true
        end
        return false
    end
else
    function _build_reaction_system(eqs, iv, unknowns, ps, defs, cevs; kwargs...)
        return ReactionSystem(
            eqs, iv, unknowns, ps;
            defaults = defs, name = gensym(:SBML), continuous_events = cevs,
            combinatoric_ratelaws = false, kwargs...
        )
    end
end

"""
    ODESystem(model::SBML.Model; include_zero_odes = true, kwargs...)

Create an `ODESystem` from an `SBML.Model`.

See also [`ReactionSystem`](@ref).
"""
function ModelingToolkit.ODESystem(
        model::SBML.Model; include_zero_odes::Bool = true,
        kwargs...
    )
    rs = ReactionSystem(model; kwargs...)
    odesys = _rs_to_odesys(rs; include_zero_odes = include_zero_odes)
    return complete(odesys)
end

# Catalyst v16 / MTK v11 removed the `convert(ODESystem, rs; ...)` path: ODESystem
# is now an `IntermediateDeprecationSystem` alias and the conversion is done via
# `Catalyst.ode_model`. Older Catalyst still uses the convert path.
@static if pkgversion(Catalyst) >= v"16"
    _rs_to_odesys(rs; kwargs...) = Catalyst.ode_model(rs; kwargs...)
else
    _rs_to_odesys(rs; kwargs...) = convert(ODESystem, rs; kwargs...)
end

function get_mappings(model::SBML.Model)
    inits = Dict(SBML.initial_amounts(model, convert_concentrations = true))
    u0map = Pair[]
    parammap = Pair[]
    initial_assignment_map = Pair[]

    for (k, v) in model.species
        var = create_symbol(k, model)
        if v.constant == true
            push!(parammap, var => inits[k])
        else
            push!(u0map, var => inits[k])
        end
    end

    for (k, v) in model.parameters
        var = create_symbol(k, model)
        if v.constant == false &&
                (SBML.seemsdefined(k, model) || is_event_assignment(k, model))
            push!(u0map, var => v.value)
        elseif v.constant == true && isnothing(v.value)  # Todo: maybe add this branch also to model.compartments
            val = model.initial_assignments[k]
            push!(parammap, var => interpret_as_num(val, model))
        else
            push!(parammap, var => v.value)
        end
    end

    for (k, v) in model.compartments
        var = create_symbol(k, model)
        if v.constant == false && SBML.seemsdefined(k, model)
            push!(u0map, var => v.size)
        else
            push!(parammap, var => v.size)
        end
    end

    for (k, v) in model.initial_assignments
        var = create_symbol(k, model)
        push!(initial_assignment_map, var => interpret_as_num(v, model))
    end
    return u0map, parammap, initial_assignment_map
end

function netstoich(id, reaction)
    netstoich = 0
    rdict = Dict(
        getproperty.(
            reaction.reactants, :species
        ) .=> getproperty.(reaction.reactants, :stoichiometry)
    )
    pdict = Dict(getproperty.(reaction.products, :species) .=> getproperty.(reaction.products, :stoichiometry))
    netstoich -= get(rdict, id, 0)
    return netstoich += get(pdict, id, 0)
end

"""
    checksupport_file(filename::String)

Check if SBML file is supported by SBMLToolkit.jl.
"""
function checksupport_file(filename::String)
    string = open(filename) do file
        read(file, String)
    end
    return checksupport_string(string)
end

"""
    checksupport_string(filename::String)

Check if SBML passed as string is supported by SBMLToolkit.jl.
"""
function checksupport_string(xml::String)
    not_implemented = [
        "listOfConstraints", "/delay",
        "<priority>",
        "factorial", "id=\"case00387\"",  # Case 00387 requires event directionality
        "id=\"case01071\"",  # require event directionality, I think
        "</eventAssignment>\n          <eventAssignment",
    ]
    for item in not_implemented
        occursin(item, xml) &&
            throw(ErrorException("SBML models with $item are not yet implemented."))
    end
    occursin("<sbml xmlns:fbc=", xml) &&
        throw(ErrorException("This model was designed for constrained-based optimisation. Please use COBREXA.jl instead of SBMLToolkit."))
    occursin("<sbml xmlns:comp=", xml) &&
        throw(ErrorException("This model uses the SBML \"comp\" package, which is not yet implemented."))
    !(occursin("<reaction", xml) || occursin("rateRule", xml)) &&
        throw(ErrorException("Models that contain neither reactions or rateRules will fail in simulation."))
    return true
end
