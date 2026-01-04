"""
Creates ContinuousVectorCallbacks
"""
function get_events(model)  # Todo: implement up or downpass and parameters
    subsdict = get_substitutions(model)  # Todo: use SUBSDICT
    evs = model.events
    mtk_evs = Pair{Vector{Equation}, Vector{Equation}}[]
    for (_, e) in evs
        trigger = SBML.extensive_kinetic_math(model, e.trigger.math)
        trigger = Symbolics.unwrap(interpret_as_num(trigger, model))
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
            var = create_var(eva.variable, IV; irreducible = true)
            math = substitute(
                Symbolics.unwrap(interpret_as_num(math, model)),
                subsdict
            )
            effect = var ~ math
            push!(mtk_evas, effect)
        end
        push!(mtk_evs, trig => mtk_evas)
    end
    !isempty(evs) &&
        @warn "SBMLToolkit currently fires events regardless of uppass or downpass trigger."
    return isempty(mtk_evs) ? nothing : mtk_evs
end
