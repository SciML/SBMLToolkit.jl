function get_rules(model)
    subsdict = get_substitutions(model)  # Todo: use SUBSDICT
    obseqs = Equation[]
    algeqs = Equation[]
    raterules = Equation[]
    for r in model.rules
        if r isa SBML.AlgebraicRule
            push!(algeqs, 0 ~ interpret_as_num(r.math, model))
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


extensive_kinetic_math(m::SBML.Model, formula::SBML.Math) = SBML.interpret_math(  # TODO: move this to SBML.jl
    formula,
    map_apply = (x, rec) -> SBML.MathApply(x.fn, rec.(x.args)),
    map_const = identity,
    map_ident = (x::SBML.MathIdent) -> begin
        haskey(m.reactions, x.id)  && return m.reactions[x.id].kinetic_math
        haskey(m.species, x.id) || return x
        sp = m.species[x.id]
        sp.only_substance_units && return x
        if isnothing(m.compartments[sp.compartment].size) &&
           !seemsdefined(sp.compartment, m)
            if m.compartments[sp.compartment].spatial_dimensions == 0
                # If the comparment ID doesn't seem directly defined anywhere
                # and it is a zero-dimensional unsized compartment, just avoid
                # any sizing questions.
                return x
            else
                # In case the compartment is expected to be defined, complain.
                throw(
                    DomainError(
                        sp.compartment,
                        "compartment size is insufficiently defined",
                    ),
                )
            end
        else
            # Now we are sure that the model either has the compartment with
            # constant size, or the definition is easily reachable. So just use
            # the compartment ID as a variable to compute the concentration (or
            # area-centration etc, with different dimensionalities) by dividing
            # it.
            return SBML.MathApply("/", [x, SBML.MathIdent(sp.compartment)])
        end
    end,
    map_lambda = (x, _) -> error(
        ErrorException("converting lambdas to extensive kinetic math is not supported"),
    ),
    map_time = identity,
    map_avogadro = identity,
    map_value = identity,
)


function get_var_and_assignment(model, rule)
    if !haskey(merge(model.species, model.compartments, model.parameters), rule.variable)
        error("Cannot find target for rule with ID `$(rule.variable)`")
    end
    var = create_var(rule.variable, IV)
    math = extensive_kinetic_math(model, rule.math)
    vc = get_volume_correction(model, rule.variable)
    if !isnothing(vc)
        math = SBML.MathApply("*", [SBML.MathIdent(vc), math])
    end
    assignment = interpret_as_num(math, model)
    if rule isa SBML.RateRule && haskey(model.species, rule.variable)
        sp = model.species[rule.variable]
        comp = model.compartments[sp.compartment]
        comp.constant == false && sp.only_substance_units == false &&
            begin
                c = create_var(sp.compartment, IV)
                assignment = c * assignment + var / c * D(c)
            end
    end
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
