myread(fn) = readSBML(fn, doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_simplify_math(doc)
end)

# assignment rule
fn = "data/00038-sbml-l3v2.xml" # this case is for observable eqs
m = myread(fn)
@named rs = ReactionSystem(m)

sys = convert(ODESystem, rs; include_zero_odes = true)
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
@named rs = ReactionSystem(m)
sys = convert(ODESystem, rs; include_zero_odes = true)
@variables t S1(t)
@parameters compartment
D = Differential(t)
@test equations(sys) == [D(S1) ~ 7*compartment]

# algebraic rule
fn = "data/00039-sbml-l3v2.xml" # this case is for observable eqs
m = myread(fn)
get_mappings(m)
@named rs = ReactionSystem(m)
sys = convert(ODESystem, rs; include_zero_odes = true)
ssys = structural_simplify(sys)
prob = ODEProblem(ssys, [], (0, 10.0))
sol = solve(prob, QNDF())
@variables t S2(t)
obsvar_sol = sol[S2]

# tests that non-constant parameters become states
fn = "data/00033-sbml-l3v2.xml"
m = myread(fn)
@named rs = ReactionSystem(m)
sys = convert(ODESystem, rs)
@variables k1(t)
@test isequal(k1, states(sys)[end])

# tests that non-constant compartments become states
# WARNING the abserr wrt the reference solution is 0.04690613469254479
# this is much higher than other tests, so simulations may be incorrect 
# however, a very similar case (00053) passed reference tests so it may be nothing
fn = "data/00051-sbml-l3v2.xml"
m = myread(fn)
@named rs = ReactionSystem(m)
sys = convert(ODESystem, rs)
@variables C(t)
@test isequal(C, states(sys)[end])

# tests that rules for non-hOSU species are multiplied with their compartment size
fn = "data/00328-sbml-l3v2.xml"
m = myread(fn)
@named rs = ReactionSystem(m)
sys = convert(ODESystem, rs; include_zero_odes = true, combinatoric_ratelaws=false)
ssys = structural_simplify(sys)
rhs = last(equations(ssys)).rhs
@variables S4(t)
@parameters k2, compartment
rhs_true = -0.5k2 * compartment
@test isequal(rhs, rhs_true)
