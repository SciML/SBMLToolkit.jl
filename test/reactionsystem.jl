cd(@__DIR__)
sbmlfile = joinpath("data", "reactionsystem_01.xml")
@parameters t, k1, c1
@variables s1(t), s2(t), s1s2(t)

kinetic_params = Dict{String,Tuple{Float64,String}}()

COMP1 = SBML.Compartment("c1", true, 3, 2.0, "nl", nothing, nothing)
SPECIES1 = SBML.Species(name = "s1", compartment = "c1", initial_amount = 1.0, substance_units = "substance", only_substance_units = true)  # Todo: Maybe not support units in initial_concentration?
SPECIES2 = SBML.Species(name = "s2", compartment = "c1", initial_amount = 1.0, substance_units = "substance/nl", only_substance_units = false)
KINETICMATH1 = SBML.MathIdent("k1")
KINETICMATH2 = SBML.MathApply("*", SBML.Math[
    SBML.MathIdent("k1"), SBML.MathIdent("s2")])
KINETICMATH3 = SBML.MathApply("-", SBML.Math[KINETICMATH2, KINETICMATH1])
REACTION1 = SBML.Reaction(
    reactants = Dict(),
    products = Dict("s1" => 1),
    kinetic_parameters = kinetic_params,
    kinetic_math = KINETICMATH1,
    reversible = false)
REACTION2 = SBML.Reaction(
    reactants = Dict("s2" => 1),
    products = Dict(),
    kinetic_parameters = kinetic_params,
    kinetic_math = KINETICMATH2,
    reversible = false)
REACTION3 = SBML.Reaction(
    reactants = Dict("s2" => 1),
    products = Dict(),
    kinetic_parameters = kinetic_params,
    kinetic_math = KINETICMATH3,
    reversible = true)
PARAM1 = SBML.Parameter(name = "k1", value = 1.0, constant = true)
MODEL1 = SBML.Model(
    parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s1" => SPECIES1),
    reactions = Dict("r1" => REACTION1),
)  # PL: For instance in the compartments dict, we may want to enforce that key and compartment.name are identical
MODEL2 = SBML.Model(
    parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s2" => SPECIES2),
    reactions = Dict("r2" => REACTION2),
)
MODEL3 = SBML.Model(
    parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s2" => SPECIES2),
    reactions = Dict("r3" => REACTION3),
)

# Test ReactionSystem constructor
rs = ReactionSystem(MODEL1)
@test isequal(Catalyst.get_eqs(rs), Catalyst.Reaction[Catalyst.Reaction(k1, nothing, [s1], nothing, [1.0])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s1])
@test isequal(Catalyst.get_ps(rs), [k1, c1])
@named rs = ReactionSystem(MODEL1)
isequal(nameof(rs), :rs)

rs = ReactionSystem(readSBML(sbmlfile))
@test isequal(Catalyst.get_eqs(rs), Catalyst.Reaction[Catalyst.Reaction(0.25c1 * k1, [s1, s2], [s1s2], [1.0, 1.0], [1.0])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s1, s1s2, s2])
@test isequal(Catalyst.get_ps(rs), [k1, c1])
@named rs = ReactionSystem(MODEL1)
isequal(nameof(rs), :rs)

rs = ReactionSystem(MODEL3)  # Contains reversible reaction
@test isequal(Catalyst.get_eqs(rs),
    Catalyst.Reaction[
        Catalyst.Reaction(0.5k1, [s2], nothing,
            [1], nothing),
        Catalyst.Reaction(k1, nothing, [s2],
            nothing, [1])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s2])
@test isequal(Catalyst.get_ps(rs), [k1, c1])

@test_nowarn convert(ModelingToolkit.ODESystem, rs)
@test_nowarn structural_simplify(convert(ModelingToolkit.ODESystem, rs))

# Test ODESystem constructor
odesys = ODESystem(MODEL1)
trueeqs = Equation[Differential(t)(s1)~k1]
@test isequal(Catalyst.get_eqs(odesys), trueeqs)
@test isequal(Catalyst.get_iv(odesys), t)
@test isequal(Catalyst.get_states(odesys), [s1])
@test isequal(Catalyst.get_ps(odesys), [k1, c1])
u0 = [s1 => 1.0]
par = [k1 => 1.0, c1 => 2.0]
@test isequal(ModelingToolkit.defaults(odesys), ModelingToolkit._merge(u0, par))  # PL: @Anand: for some reason this does not work with `Catalyst.get_default()`
@named odesys = ODESystem(MODEL1)
isequal(nameof(odesys), :odesys)
@test_nowarn structural_simplify(odesys)

odesys = ODESystem(readSBML(sbmlfile))
m = readSBML(sbmlfile)
trueeqs = Equation[Differential(t)(s1)~-0.25c1*k1*s1*s2,
    Differential(t)(s1s2)~0.25c1*k1*s1*s2,
    Differential(t)(s2)~-0.25c1*k1*s1*s2]
@test isequal(Catalyst.get_eqs(odesys), trueeqs)
@test isequal(Catalyst.get_iv(odesys), t)
@test isequal(Catalyst.get_states(odesys), [s1, s1s2, s2])
@test isequal(Catalyst.get_ps(odesys), [k1, c1])
u0 = [s1 => 2 * 1.0, s2 => 2 * 1.0, s1s2 => 2 * 1.0]
par = [k1 => 1.0, c1 => 2.0]
@test isequal(ModelingToolkit.defaults(odesys), ModelingToolkit._merge(u0, par))
@named odesys = ODESystem(MODEL1)
isequal(nameof(odesys), :odesys)

@test_nowarn ODEProblem(odesys, [], [0.0, 1.0], [])

# Test ODEProblem
oprob = ODEProblem(ODESystem(MODEL1), [], [0.0, 1.0], [])
sol = solve(oprob, Tsit5())
@test isapprox(sol.u, [[1.0], [2.0]])

@test_nowarn ODEProblem(ODESystem(readSBML(sbmlfile)), [], [0.0, 1.0], [])

# Test checksupport
# @test_nowarn SBMLToolkit.checksupport(MODEL1)
# sbc = SBML.Species("sbc", "c1", true, nothing, nothing, (1., "substance"), nothing, true)
# mod = SBML.Model(Dict("k1" => 1.), Dict(), Dict("c1" => COMP1), Dict("sbc" => sbc), Dict("r1" => REACTION1), Dict(), Dict())
# @test_logs (:warn, "Species sbc has `boundaryCondition` or is `constant`. This will lead to wrong results when simulating the `ReactionSystem`.") SBMLToolkit.checksupport(mod)

# Test checksupport
@test_nowarn SBMLToolkit.checksupport(sbmlfile)
@test_throws ErrorException SBMLToolkit.checksupport(joinpath("data", "unsupported.sbml"))



# Test _get_substitutions
truesubs = Dict(Num(Symbolics.variable(:c1; T = Real)) => c1,
    Num(Symbolics.variable(:k1; T = Real)) => k1,
    Num(Symbolics.variable(:s1; T = Real)) => s1)
subs = SBMLToolkit._get_substitutions(MODEL1)
@test isequal(subs, truesubs)

# Test mtk_reactions
reaction = SBMLToolkit.mtk_reactions(MODEL1)[1]
truereaction = Catalyst.Reaction(k1, nothing, [s1], nothing, [1])  # Todo: implement Sam's suggestion on mass action kinetics
@test isequal(reaction, truereaction)

km = SBML.MathTime("x")
reac = SBML.Reaction(
    reactants = Dict("s1" => 1),
    products = Dict(),
    kinetic_parameters = kinetic_params,
    kinetic_math = km,
    reversible = false)
mod = SBML.Model(
    parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s1" => SPECIES1),
    reactions = Dict("r1" => reac),
)

@test isequal(Catalyst.DEFAULT_IV, SBMLToolkit.mtk_reactions(mod)[1].rate)


# test from_noncombinatoric
@test isequal(4k1, SBMLToolkit.from_noncombinatoric(k1, [2, 2], false))
@test isequal(k1, SBMLToolkit.from_noncombinatoric(k1, nothing, false))

# Test use_rate
@test isequal(SBMLToolkit.use_rate(k1 * s1, [s1], [1]), (k1, false))  # Case hOSU=true
@test isequal(SBMLToolkit.use_rate(k1 * s1 * s2 / (c1 + s2), [s1], [1]), (k1 * s1 * s2 / (c1 + s2), true))  # Case Michaelis-Menten kinetics

# Test getreagents
@test isequal((nothing, [s1], nothing, [1.0]), SBMLToolkit.getreagents(REACTION1.reactants, REACTION1.products, MODEL1))

constspec = SBML.Species(name = "constspec", compartment = "c1", boundary_condition = true, initial_amount = 1.0, substance_units = "substance", only_substance_units = true)
ncs = SBMLToolkit.create_var("constspec", Catalyst.DEFAULT_IV)
kineticmath = SBML.MathApply("-", SBML.Math[
    SBML.MathApply("*", SBML.Math[
        SBML.MathIdent("k1"),
        SBML.MathIdent("constspec")]),
    SBML.MathIdent("k1")])
constreac = SBML.Reaction(
    reactants = Dict("constspec" => 1),
    products = Dict(),
    kinetic_parameters = kinetic_params,
    kinetic_math = kineticmath,
    reversible = false)

constmod = SBML.Model(
    parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("constspec" => constspec),
    reactions = Dict("r1" => constreac),
)
@test isequal(([ncs], [ncs], [1.0], [1.0]), SBMLToolkit.getreagents(constreac.reactants, constreac.products, constmod))
@test isequal((nothing, nothing, nothing, nothing), SBMLToolkit.getreagents(constreac.reactants, constreac.products, constmod; rev = true))

# Test getunidirectionalcomponents
km = SBML.MathApply("-", SBML.Math[KINETICMATH1, SBML.MathIdent("c1")])
sm = convert(Num, km)
kl = SBMLToolkit.getunidirectionalcomponents(sm)
@test isequal(kl, (k1, c1))

km = SBML.MathApply("-", SBML.Math[KINETICMATH1, KINETICMATH2])
sm = convert(Num, km)
fw, rv = SBMLToolkit.getunidirectionalcomponents(sm)
rv = substitute(rv, Dict(SBMLToolkit.create_var("s2") => SBMLToolkit.create_var("s2", Catalyst.DEFAULT_IV)))
@test isequal((fw, rv), (k1, k1 * s2))

km = SBML.MathIdent("s1s2")
sm1 = convert(Num, km)
sm2 = sm - sm1
@test_throws ErrorException SBMLToolkit.getunidirectionalcomponents(sm2)
@test_throws ErrorException SBMLToolkit.getunidirectionalcomponents(:k1)

# Test get_paramap
trueparamap = [k1 => 1.0, c1 => 2.0]
paramap = SBMLToolkit.get_paramap(MODEL1)
@test isequal(paramap, trueparamap)

# Test getmassaction
@test isequal(SBMLToolkit.getmassaction(k1 * s1, [s1], [1]), k1)  # Case hOSU=true
@test isequal(SBMLToolkit.getmassaction(k1 * s1 / c1, [s1], [1]), k1 / c1)  # Case hOSU=false
@test isequal(SBMLToolkit.getmassaction(k1 + c1, nothing, nothing), k1 + c1)  # Case zero order kineticLaw
@test isnan(SBMLToolkit.getmassaction(k1 * s1 * s2 / (c1 + s2), [s1], [1]))  # Case Michaelis-Menten kinetics
@test isnan(SBMLToolkit.getmassaction(k1 * s1 * Catalyst.DEFAULT_IV, [s1], [1]))  # Case kineticLaw with time

@test isnan(SBMLToolkit.getmassaction(k1 * s1 * s2, [s1], [1]))
@test isnan(SBMLToolkit.getmassaction(k1 + c1, [s1], [1]))
@test_throws ErrorException SBMLToolkit.getmassaction(k1, nothing, [1])

# default test
@test ModelingToolkit.defaults(m) == ModelingToolkit.defaults(ReactionSystem(m))
