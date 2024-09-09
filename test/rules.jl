using SBMLToolkit
using SBML, SBMLToolkitTestSuite
using Catalyst, ModelingToolkit, OrdinaryDiffEq
using Test

const IV = default_t()
@parameters k1, compartment
@species S1(IV), S2(IV)

function readmodel(sbml)
    SBMLToolkit.readSBMLFromString(sbml, doc -> begin
        set_level_and_version(3, 2)(doc)
        convert_promotelocals_expandfuns(doc)
    end)
end

# Test get_rules
sbml, _, _ = SBMLToolkitTestSuite.read_case("00029")  # assignmentRule
m = complete(readmodel(sbml))
a, o, r = SBMLToolkit.get_rules(m)
o_true = [S1 ~ 7 * compartment]
@test isequal(o, o_true)

sbml, _, _ = SBMLToolkitTestSuite.read_case("00031")  # rateRule
m = complete(readmodel(sbml))
a, o, r = SBMLToolkit.get_rules(m)
r_true = [default_time_deriv()(S1) ~ 7 * compartment]
@test isequal(r, r_true)

sbml, _, _ = SBMLToolkitTestSuite.read_case("00039")  # algebraicRule
m = complete(readmodel(sbml))
a, o, r = SBMLToolkit.get_rules(m)
a_true = [0 ~ S1 + S2 - k1]
@test isequal(a, a_true)

# Test get_var_and_assignment
sbml, _, _ = SBMLToolkitTestSuite.read_case("00031")
m = complete(readmodel(sbml))
var, assignment = SBMLToolkit.get_var_and_assignment(m, m.rules[1])
var_true = S1
assignment_true = 7 * compartment
@test isequal(var, var_true)
@test isequal(assignment, assignment_true)

r = SBML.AssignmentRule("S2", SBML.MathVal(1))
@test_throws ErrorException("Cannot find target for rule with ID `S2`") SBMLToolkit.get_var_and_assignment(
    m,
    r)

# Test get_volume_correction
vc = SBMLToolkit.get_volume_correction(m, "notaspecies")
@test isnothing(vc)

vc = SBMLToolkit.get_volume_correction(m, "S1")
@test isequal(vc, "compartment")

sbml, _, _ = SBMLToolkitTestSuite.read_case("00060")  # hOSU="true" species
m = complete(readmodel(sbml))
vc = SBMLToolkit.get_volume_correction(m, "S1")
@test isnothing(vc)

# tests that non-constant parameters become variables
sbml, _, _ = SBMLToolkitTestSuite.read_case("00033")
m = complete(readmodel(sbml))
@named sys = ODESystem(m)
@species k1(IV)
@test isequal(k1, unknowns(sys)[end])

# tests that non-constant compartments become variables
sbml, _, _ = SBMLToolkitTestSuite.read_case("00051")  # hOSU="true" species
m = complete(readmodel(sbml))
@named sys = ODESystem(m)
@species C(IV)
@test isequal(C, unknowns(sys)[end])
