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
    algeqs, obseqs, raterules = map(x -> substitute(x, subsdict),
                                    (algeqs, obseqs, raterules))
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
        throw(DomainError(sp.compartment,
                          "compartment size is insufficiently defined"))
    sp.compartment
end
