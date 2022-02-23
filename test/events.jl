myread(fn) = readSBML(fn, doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_simplify_math(doc)
end)

# 1 single trigger, single affect
fn = "data/00026-sbml-l3v2.xml"
m = myread(fn)
@named sys = ODESystem(m)
ssys = structural_simplify(sys)
@test length(ModelingToolkit.get_continuous_events(ssys)) == 1
prob = ODEProblem(ssys, [], (0, 5.0); saveat = 0:0.1:5.0)
sol = solve(prob, Tsit5())
@test sol.destats.ncondition > 0

# multiple events 
fn = "data/00041-sbml-l3v2.xml"
m = myread(fn)
@named sys = ODESystem(m)
ssys = structural_simplify(sys)
@test length(ModelingToolkit.get_continuous_events(ssys)) == 2
prob = ODEProblem(ssys, [], (0, 5.0); saveat = 0:0.1:5.0)
sol = solve(prob, Tsit5())
@test sol.destats.ncondition > 0
