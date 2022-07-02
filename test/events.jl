using SBMLToolkit
using Catalyst, SBMLToolkitTestSuite
using Test

const IV = Catalyst.DEFAULT_IV
@parameters compartment
@variables S1(IV) S2(IV)

readmodel(sbml) = SBMLToolkit.readSBMLFromString(sbml, doc -> begin
                            set_level_and_version(3, 2)(doc)
                            convert_simplify_math(doc)
                        end)

# Test get_events
sbml, _, _ = SBMLToolkitTestSuite.read_case("00001")
m = readmodel(sbml)
@test isnothing(SBMLToolkit.get_events(m))


@variables
sbml, _, _ = SBMLToolkitTestSuite.read_case("00026")  # 1 single trigger, single affect
m = readmodel(sbml)
events = SBMLToolkit.get_events(m)
events_true = [[S1 / compartment ~ 0.1] => [S1 ~ compartment]]
@test isequal(events, events_true)

sbml, _, _ = SBMLToolkitTestSuite.read_case("00041")  # multiple events 
m = readmodel(sbml)
events = SBMLToolkit.get_events(m)
events_true =  [[S1 / compartment ~ 0.1] => [S1 ~ compartment],
                [S2 / compartment ~ 0.5] => [S2 ~ 0]]
@test isequal(events, events_true)

# # 1 single trigger, single affect
# fn = "data/00026-sbml-l3v2.xml"
# m = myread(fn)
# @named rs = ReactionSystem(m)
# sys = convert(ODESystem, rs)
# ssys = structural_simplify(sys)
# @test length(ModelingToolkit.get_continuous_events(ssys)) == 1
# prob = ODEProblem(ssys, [], (0, 5.0); saveat = 0:0.1:5.0)
# sol = solve(prob, Tsit5())
# @test sol.destats.ncondition > 0

# # multiple events 
# fn = "data/00041-sbml-l3v2.xml"
# m = myread(fn)
# @named sys = ODESystem(m)
# ssys = structural_simplify(sys)
# @test length(ModelingToolkit.get_continuous_events(ssys)) == 2
# prob = ODEProblem(ssys, [], (0, 5.0); saveat = 0:0.1:5.0)
# sol = solve(prob, Tsit5())
# @test sol.destats.ncondition > 0
