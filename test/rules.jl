using SBMLToolkit
using SBMLToolkitTestSuite
using Catalyst, ModelingToolkit, OrdinaryDiffEq
using Test

const IV = Catalyst.DEFAULT_IV
@parameters k1, compartment
@variables S1(IV), S2(IV)

# myread(fn) = readSBML(fn, doc -> begin
#                           set_level_and_version(3, 2)(doc)
#                           convert_simplify_math(doc)
#                       end)
readmodel(sbml) = SBMLToolkit.readSBMLFromString(sbml, doc -> begin
                            set_level_and_version(3, 2)(doc)
                            convert_simplify_math(doc)
                        end)

# Test get_rules
sbml, _, _ = SBMLToolkitTestSuite.read_case("00029")  # assignmentRule
m = readmodel(sbml)
a, o, r = SBMLToolkit.get_rules(m)
o_true = [S1 ~ 7 * compartment]
@test isequal(o, o_true)

sbml, _, _ = SBMLToolkitTestSuite.read_case("00031")  # rateRule
m = readmodel(sbml)
a, o, r = SBMLToolkit.get_rules(m)
r_true = [Differential(IV)(S1) ~ 7 * compartment]
@test isequal(r, r_true)

sbml, _, _ = SBMLToolkitTestSuite.read_case("00039")  # algebraicRule
m = readmodel(sbml)
a, o, r = SBMLToolkit.get_rules(m)
a_true = [0 ~ S1 + S2 - k1]
@test isequal(a, a_true)

# Test get_var_and_assignment
sbml, _, _ = SBMLToolkitTestSuite.read_case("00031")
m = readmodel(sbml)
var, assignment = SBMLToolkit.get_var_and_assignment(m, m.rules[1])
var_true = S1
assignment_true = 7 * compartment
@test isequal(var, var_true)
@test isequal(assignment, assignment_true)

# Todo: add rule where the call to extensive_kinetic_math is needed

# Test get_volume_correction
vc = SBMLToolkit.get_volume_correction(m, "notaspecies")
@test isnothing(vc)

vc = SBMLToolkit.get_volume_correction(m, "S1")
@test isequal(vc, "compartment")

sbml, _, _ = SBMLToolkitTestSuite.read_case("00060")  # hOSU="true" species
m = readmodel(sbml)
vc = SBMLToolkit.get_volume_correction(m, "S1")
@test isnothing(vc)

# # assignment rule
# fn = "data/00038-sbml-l3v2.xml" # this case is for observable eqs
# m = myread(fn)
# @named rs = ReactionSystem(m)

# sys = convert(ODESystem, rs; include_zero_odes = true)
# @test length(equations(sys)) == 3
# ssys = structural_simplify(sys)
# @test length(observed(ssys)) == 1
# prob = ODEProblem(ssys, [], (0, 10.0))
# sol = solve(prob, Tsit5())

# @variables t S3(t)
# obsvar_sol = sol[S3]
# @test !all(isequal(first(obsvar_sol)), obsvar_sol) # check its nonconstant

# # rate rule
# fn = "data/00031-sbml-l3v2.xml" # this case is for observable eqs
# m = myread(fn)
# @named rs = ReactionSystem(m)
# sys = convert(ODESystem, rs; include_zero_odes = true)
# @variables t S1(t)
# @parameters compartment
# D = Differential(t)
# @test equations(sys) == [D(S1) ~ 7 * compartment]

# # algebraic rule
# fn = "data/00039-sbml-l3v2.xml" # this case is for observable eqs
# m = myread(fn)
# @named rs = ReactionSystem(m)
# sys = convert(ODESystem, rs; include_zero_odes = true)
# ssys = structural_simplify(sys)
# prob = ODEProblem(ssys, [], (0, 10.0))
# sol = solve(prob, QNDF())
# @variables t S2(t)
# obsvar_sol = sol[S2]

# tests that non-constant parameters become states
sbml, _, _ = SBMLToolkitTestSuite.read_case("00033")
m = readmodel(sbml)
@named sys = ODESystem(m)
@variables k1(IV)
@test isequal(k1, states(sys)[end])

# tests that non-constant compartments become states
# WARNING the abserr wrt the reference solution is 0.04690613469254479
# this is much higher than other tests, so simulations may be incorrect 
# however, a very similar case (00053) passed reference tests so it may be nothing
sbml, _, _ = SBMLToolkitTestSuite.read_case("00051")  # hOSU="true" species
m = readmodel(sbml)
@named sys = ODESystem(m)
@variables C(IV)
@test isequal(C, states(sys)[end])

# # tests that rules for non-hOSU species are multiplied with their compartment size
# fn = "data/00328-sbml-l3v2.xml"
# m = myread(fn)
# @named rs = ReactionSystem(m)
# sys = convert(ODESystem, rs; include_zero_odes = true, combinatoric_ratelaws = false)
# ssys = structural_simplify(sys)
# rhs = last(equations(ssys)).rhs
# @variables S4(t)
# @parameters k2, compartment
# rhs_true = -0.5k2 * compartment
# @test isequal(rhs, rhs_true)
