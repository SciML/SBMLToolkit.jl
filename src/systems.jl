"""
    ReactionSystem(model::SBML.Model; kwargs...)

Create a `ReactionSystem` from an `SBML.Model`.

See also [`ODESystem`](@ref).
"""
function Catalyst.ReactionSystem(model::SBML.Model; kwargs...)  # Todo: requires unique parameters (i.e. SBML must have been imported with localParameter promotion in libSBML)
    # length(model.events) > 0 ? error("Model contains events. Please import with `ODESystem(model)`") : nothing  @Anand: how to suppress this when called from ODESystem
    rxs = get_reactions(model)
    u0map, parammap = get_mappings(model)
    defs = ModelingToolkit._merge(Dict(u0map), Dict(parammap))
    algrules, obsrules, raterules = get_rules(model)
    obsrules_rearranged = Equation[]
    for o in obsrules
        rhs = o.rhs
        for r in raterules
            if isequal(rhs, r.lhs)
                rhs = r.rhs
            end
        end
        defs[o.lhs] = substitute(rhs,
                                 ModelingToolkit._merge(defs,
                                                        Dict(Catalyst.DEFAULT_IV.val => 0)))
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
        defs[o.lhs] = substitute(rhs,
                                 ModelingToolkit._merge(defs,
                                                        Dict(Catalyst.DEFAULT_IV.val => 0)))
        push!(raterules_subs, rhs ~ o.lhs)
    end
    if haskey(kwargs, :defaults)
        defs = ModelingToolkit._merge(defs, kwargs[:defaults])
        kwargs = filter(x -> !isequal(first(x), :defaults), kwargs)
    end
    ReactionSystem([rxs..., algrules..., raterules_subs..., obsrules_rearranged...],
                   IV, first.(u0map), first.(parammap);
                   defaults = defs, name = gensym(:SBML),
                   continuous_events = get_events(model),
                   combinatoric_ratelaws = false, kwargs...)
end

"""
    ODESystem(model::SBML.Model; include_zero_odes = true, kwargs...)

Create an `ODESystem` from an `SBML.Model`.

See also [`ReactionSystem`](@ref).
"""
function ModelingToolkit.ODESystem(model::SBML.Model; include_zero_odes = true, kwargs...)
    rs = ReactionSystem(model; kwargs...)
    convert(ODESystem, rs; include_zero_odes = include_zero_odes)
end

function get_mappings(model::SBML.Model)
    inits = Dict(SBML.initial_amounts(model, convert_concentrations = true))
    u0map = Pair[]
    parammap = Pair[]
    for (k, v) in model.species
        if v.constant == true
            var = create_param(k; isconstantspecies = true)
            push!(parammap, var => inits[k])
        else
            var = create_var(k, IV;
                             isbcspecies = has_rule_type(k, model, SBML.RateRule) ||
                                           has_rule_type(k, model, SBML.AssignmentRule) ||
                                           (has_rule_type(k, model, SBML.AlgebraicRule) &&
                                            (all([netstoich(k, r) == 0
                                                  for r in values(model.reactions)]) ||
                                             v.boundary_condition == true)))  # To remove species that are otherwise defined
            push!(u0map, var => inits[k])
        end
    end
    for (k, v) in model.parameters
        if v.constant == false &&
           (SBML.seemsdefined(k, model) || is_event_assignment(k, model))
            var = create_var(k, IV; isbcspecies = true)
            push!(u0map, var => v.value)
        else
            var = create_param(k)
            push!(parammap, var => v.value)
        end
    end
    for (k, v) in model.compartments
        if v.constant == false && SBML.seemsdefined(k, model)
            var = create_var(k, IV; isbcspecies = true)
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
