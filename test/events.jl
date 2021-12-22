myread(fn) = readSBML(fn, doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_simplify_math(doc)
end)

fn = "data/00041-sbml-l3v2.xml"
m = myread(fn)
sys = ODESystem(m) # internally calls structural_simplify
evs = ModelingToolkit.get_continuous_events(sys)
@test length(evs) == 2
prob = ODEProblem(sys, [], (0, 5.0); saveat = 0:0.1:5.0)
sol = solve(prob, Tsit5())
@test sol.destats.ncondition > 0
