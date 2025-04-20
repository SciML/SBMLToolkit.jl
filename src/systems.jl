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
    readSBML(sbmlfile, importdefaults)
end

"""
    readSBML(sbmlfile::String, ::ReactionSystemImporter)

Create a `Catalyst.ReactionSystem` from an SBML file, using the default import settings.

See also [`Model`](@ref) and [`ReactionSystemImporter`](@ref).
"""
function SBML.readSBML(sbmlfile::String, ::ReactionSystemImporter; kwargs...)  # Returns a Catalyst.ReactionSystem
    ReactionSystem(readSBML(sbmlfile::String, DefaultImporter()), kwargs...)
end

"""
    readSBML(sbmlfile::String, ::ODESystemImporter)

Create a `ModelingToolkit.ODESystem` from an SBML file, using the default import settings.

See also [`Model`](@ref) and [`ODESystemImporter`](@ref).
"""
function SBML.readSBML(sbmlfile::String, ::ODESystemImporter;
        include_zero_odes::Bool = true, kwargs...)  # Returns an MTK.ODESystem
    odesys = convert(ODESystem, readSBML(sbmlfile, ReactionSystemImporter(), kwargs...),
        include_zero_odes = include_zero_odes)
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
    rs = ReactionSystem(
        [rxs..., algrules..., raterules_subs..., obsrules_rearranged..., zero_rates...],
        IV, first.(u0map), first.(parammap);
        defaults = defs, name = gensym(:SBML),
        continuous_events = get_events(model),
        combinatoric_ratelaws = false, kwargs...)
    return rs
end

"""
    ODESystem(model::SBML.Model; include_zero_odes = true, kwargs...)

Create an `ODESystem` from an `SBML.Model`.

See also [`ReactionSystem`](@ref).
"""
function ModelingToolkit.ODESystem(model::SBML.Model; include_zero_odes::Bool = true,
        kwargs...)
    rs = ReactionSystem(model; kwargs...)
    odesys = convert(ODESystem, rs; include_zero_odes = include_zero_odes)
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
    u0map, parammap, initial_assignment_map
end

function netstoich(id, reaction)
    netstoich = 0
    rdict = Dict(getproperty.(reaction.reactants, :species) .=>
        getproperty.(reaction.reactants, :stoichiometry))
    pdict = Dict(getproperty.(reaction.products, :species) .=>
        getproperty.(reaction.products, :stoichiometry))
    netstoich -= get(rdict, id, 0)
    netstoich += get(pdict, id, 0)
end

"""
    checksupport_file(filename::String)

Check if SBML file is supported by SBMLToolkit.jl.
"""
function checksupport_file(filename::String)
    string = open(filename) do file
        read(file, String)
    end
    checksupport_string(string)
end

"""
    checksupport_string(filename::String)

Check if SBML passed as string is supported by SBMLToolkit.jl.
"""
function checksupport_string(xml::String)
    not_implemented = ["listOfConstraints", "/delay",
        "<priority>",
        "factorial", "id=\"case00387\"",  # Case 00387 requires event directionality
        "id=\"case01071\"",  # require event directionality, I think
        "</eventAssignment>\n          <eventAssignment"]
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
    true
end
