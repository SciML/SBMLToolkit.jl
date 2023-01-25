using SBMLToolkit
using Catalyst, SBMLToolkitTestSuite
using Test

const IV = Catalyst.DEFAULT_IV
@parameters compartment
@variables S1(IV) S2(IV)

function readmodel(sbml)
    SBMLToolkit.readSBMLFromString(sbml, doc -> begin
                                       set_level_and_version(3, 2)(doc)
                                       convert_promotelocals_expandfuns(doc)
                                   end)
end

# Test get_events
sbml, _, _ = SBMLToolkitTestSuite.read_case("00001")
m = readmodel(sbml)
@test isnothing(SBMLToolkit.get_events(m))

sbml, _, _ = SBMLToolkitTestSuite.read_case("00026")  # 1 single trigger, single affect
m = readmodel(sbml)
events = SBMLToolkit.get_events(m)
events_true = [[S1 / compartment ~ 0.1] => [S1 ~ compartment]]
@test isequal(events, events_true)

sbml, _, _ = SBMLToolkitTestSuite.read_case("00041")  # multiple events 
m = readmodel(sbml)
events = SBMLToolkit.get_events(m)
events_true = [[S1 / compartment ~ 0.1] => [S1 ~ compartment],
    [S2 / compartment ~ 0.5] => [S2 ~ 0]]
@test isequal(events, events_true)
