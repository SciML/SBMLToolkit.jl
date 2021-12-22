myread(fn) = readSBML(fn, doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_simplify_math(doc)
end)

# assignment rule
fn = "data/00038-sbml-l3v2.xml" # this case is for observable eqs
m = myread(fn)
rs = ReactionSystem(m)

sys = convert(ODESystem, rs; include_zero_odes = false)
@test length(equations(sys)) == 3
ssys = structural_simplify(sys)
@test length(observed(ssys)) == 1
prob = ODEProblem(ssys, [], (0, 10.0))
sol = solve(prob, Tsit5())

@variables t S3(t)
obsvar_sol = sol[S3]
@test !all(isequal(first(obsvar_sol)), obsvar_sol) # check its nonconstant

# rate rule
fn = "data/00031-sbml-l3v2.xml" # this case is for observable eqs
m = myread(fn)
rs = ReactionSystem(m)
sys = convert(ODESystem, rs; include_zero_odes = false)
@variables t S1(t)
D = Differential(t)
@test equations(sys) == [D(S1) ~ 7]

# algebraic rule
fn = "data/00039-sbml-l3v2.xml" # this case is for observable eqs
m = myread(fn)
rs = ReactionSystem(m)
sys = convert(ODESystem, rs; include_zero_odes = false)
ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [], (0, 10.0))
sol = solve(prob, Tsit5())
@variables t S2(t)
obsvar_sol = sol[S2]
