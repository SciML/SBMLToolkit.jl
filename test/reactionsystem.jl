sbmlfile = joinpath("data", "reactionsystem_01.xml")
@parameters t, k1, c1
@variables s1(t), s2(t), s1s2(t)

COMP1 = SBML.Compartment("c1", true, 3, 2., "nl") 
SPECIES1 = SBML.Species("s1", "c1", false, nothing, nothing, (1., "substance"), nothing, true)  # Todo: Maybe not support units in initial_concentration?
SPECIES2 = SBML.Species("s2", "c1", false, nothing, nothing, nothing, (1., "substance/nl"), false)
KINETICMATH1 = SBML.MathIdent("k1")
KINETICMATH2 = SBML.MathApply("*", SBML.Math[
SBML.MathIdent("k1"), SBML.MathIdent("s2")])
KINETICMATH3 = SBML.MathApply("-", SBML.Math[KINETICMATH2, KINETICMATH1])
REACTION1 = SBML.Reaction(Dict("s1" => 1), Dict(), (NaN, ""), (NaN, ""), NaN,
                        nothing, KINETICMATH1, false)
REACTION2 = SBML.Reaction(Dict(), Dict("s2" => 1), (NaN, ""), (NaN, ""), NaN,
                        nothing, KINETICMATH2, false)
REACTION3 = SBML.Reaction(Dict(), Dict("s2" => 1), (NaN, ""), (NaN, ""), NaN,
                        nothing, KINETICMATH3, true)
MODEL1 = SBML.Model(Dict("k1" => 1.), Dict(), Dict("c1" => COMP1), Dict("s1" => SPECIES1), Dict("r1" => REACTION1), Dict(), Dict())  # PL: For instance in the compartments dict, we may want to enforce that key and compartment.name are identical
MODEL2 = SBML.Model(Dict("k1" => 1.), Dict(), Dict("c1" => COMP1), Dict("s2" => SPECIES2), Dict("r2" => REACTION2), Dict(), Dict())
MODEL3 = SBML.Model(Dict("k1" => 1.), Dict(), Dict("c1" => COMP1), Dict("s2" => SPECIES2), Dict("r3" => REACTION3), Dict(), Dict())

# Test ReactionSystem constructor
rs = ReactionSystem(MODEL1)
@test isequal(Catalyst.get_eqs(rs), ModelingToolkit.Reaction[ModelingToolkit.Reaction(k1, nothing, [s1], nothing, [1.])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s1])
@test isequal(Catalyst.get_ps(rs), [k1,c1])
@named rs = ReactionSystem(MODEL1)
isequal(nameof(rs), :rs)

rs = ReactionSystem(readSBML(sbmlfile))
@test isequal(Catalyst.get_eqs(rs), ModelingToolkit.Reaction[ModelingToolkit.Reaction(0.25c1*k1, [s1, s2], [s1s2], [1., 1.], [1.])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s1, s1s2, s2])
@test isequal(Catalyst.get_ps(rs), [k1,c1])
@named rs = ReactionSystem(MODEL1)
isequal(nameof(rs), :rs)

rs = ReactionSystem(MODEL3)  # Contains reversible reaction
@test isequal(Catalyst.get_eqs(rs),
            ModelingToolkit.Reaction[
                ModelingToolkit.Reaction(0.5k1, [s2], nothing,
                    [1.], nothing),
                ModelingToolkit.Reaction(k1, nothing, [s2],
                    nothing, [1.])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s2])
@test isequal(Catalyst.get_ps(rs), [k1, c1])

@test_nowarn convert(ModelingToolkit.ODESystem, rs)
@test_nowarn structural_simplify(convert(ModelingToolkit.ODESystem, rs))

# Test ODESystem constructor
odesys = ODESystem(MODEL1)
trueeqs = Equation[Differential(t)(s1) ~ k1]
@test isequal(Catalyst.get_eqs(odesys), trueeqs)
@test isequal(Catalyst.get_iv(odesys), t)
@test isequal(Catalyst.get_states(odesys), [s1])
@test isequal(Catalyst.get_ps(odesys), [k1, c1])
u0 = [s1 => 1.]
par = [k1 => 1., c1 => 2.]
@test isequal(ModelingToolkit.defaults(odesys), ModelingToolkit._merge(u0, par))  # PL: @Anand: for some reason this does not work with `Catalyst.get_default()`
@named odesys = ODESystem(MODEL1)
isequal(nameof(odesys), :odesys)
@test_nowarn structural_simplify(odesys)

odesys = ODESystem(readSBML(sbmlfile))
m = readSBML(sbmlfile)
trueeqs = Equation[Differential(t)(s1) ~ -0.25c1 * k1 * s1 * s2,
                Differential(t)(s1s2) ~ 0.25c1 * k1 * s1 * s2,
                Differential(t)(s2) ~ -0.25c1 * k1 * s1 * s2]
@test isequal(Catalyst.get_eqs(odesys), trueeqs)
@test isequal(Catalyst.get_iv(odesys), t)
@test isequal(Catalyst.get_states(odesys), [s1, s1s2, s2])
@test isequal(Catalyst.get_ps(odesys), [k1, c1])
u0 = [s1 => 2*1., s2 => 2*1., s1s2 => 2*1.]
par = [k1 => 1., c1 => 2.]
@test isequal(ModelingToolkit.defaults(odesys), ModelingToolkit._merge(u0, par))
@named odesys = ODESystem(MODEL1)
isequal(nameof(odesys), :odesys)

@test_nowarn ODEProblem(odesys, [], [0., 1.], [])

# Test ODEProblem
oprob = ODEProblem(ODESystem(MODEL1), [], [0., 1.], [])
sol = solve(oprob, Tsit5())
@test isapprox(sol.u, [[1.], [2.]])

@test_nowarn ODEProblem(ODESystem(readSBML(sbmlfile)), [], [0., 1.], [])

# Test checksupport
# @test_nowarn SBMLToolkit.checksupport(MODEL1)
# sbc = SBML.Species("sbc", "c1", true, nothing, nothing, (1., "substance"), nothing, true)   
# mod = SBML.Model(Dict("k1" => 1.), Dict(), Dict("c1" => COMP1), Dict("sbc" => sbc), Dict("r1" => REACTION1), Dict(), Dict())
# @test_logs (:warn, "Species sbc has `boundaryCondition` or is `constant`. This will lead to wrong results when simulating the `ReactionSystem`.") SBMLToolkit.checksupport(mod)

# Test checksupport
@test_nowarn SBMLToolkit.checksupport(sbmlfile)
@test_throws ErrorException SBMLToolkit.checksupport(joinpath("data", "unsupported.sbml"))

# Test _get_substitutions
truesubs = Dict(Num(Variable(:c1)) => c1,
        Num(Variable(:k1)) => k1,
        Num(Variable(:s1)) => s1)
subs = SBMLToolkit._get_substitutions(MODEL1)
@test isequal(subs, truesubs)

# Test mtk_reactions
reaction = SBMLToolkit.mtk_reactions(MODEL1)[1]
truereaction = ModelingToolkit.Reaction(k1, nothing, [s1], nothing, [1]; only_use_rate=true)  # Todo: implement Sam's suggestion on mass action kinetics
@test isequal(reaction, truereaction)

km = SBML.MathTime("x")
reac = SBML.Reaction(Dict("s1" => 1), Dict(), (NaN, ""), (NaN, ""), NaN,
nothing, km, false)
mod = SBML.Model(Dict(), Dict(), Dict("c1" => COMP1), Dict("s1" => SPECIES1), Dict("r1" => reac), Dict(), Dict())
@test isequal(Catalyst.DEFAULT_IV, SBMLToolkit.mtk_reactions(mod)[1].rate)


# test from_noncombinatoric
@test isequal(4k1, SBMLToolkit.from_noncombinatoric(k1, [2,2], false))
@test isequal(k1, SBMLToolkit.from_noncombinatoric(k1, nothing, false))

# Test use_rate
@test isequal(SBMLToolkit.use_rate(k1*s1, [s1], [1]), (k1, false))  # Case hOSU=true
@test isequal(SBMLToolkit.use_rate(k1*s1*s2/(c1+s2), [s1], [1]), (k1*s1*s2/(c1+s2), true))  # Case Michaelis-Menten kinetics

# Test getreagents
@test isequal((nothing, [s1], nothing, [1.]), SBMLToolkit.getreagents(REACTION1.reactants, REACTION1.products, MODEL1))

constspec = SBML.Species("constspec", "c1", true, nothing, nothing, (1., "substance"), nothing, true)
ncs = SBMLToolkit.create_var("constspec",Catalyst.DEFAULT_IV)
kineticmath = SBML.MathApply("-", SBML.Math[
SBML.MathApply("*", SBML.Math[
    SBML.MathIdent("k1"),
    SBML.MathIdent("constspec")]),
SBML.MathIdent("k1")])
constreac = SBML.Reaction(Dict(), Dict("constspec" => 1), (NaN, ""), (NaN, ""), NaN,nothing, kineticmath, false)
constmod = SBML.Model(Dict("k1" => 1.), Dict(), Dict("c1" => COMP1), Dict("constspec" => constspec), Dict("r1" => constreac), Dict(), Dict())
@test isequal(([ncs], [ncs], [1.], [1.]), SBMLToolkit.getreagents(constreac.reactants, constreac.products, constmod))
@test isequal((nothing, nothing, nothing, nothing), SBMLToolkit.getreagents(constreac.reactants, constreac.products, constmod))

# Test getunidirectionalcomponents
km = SBML.MathApply("-", SBML.Math[KINETICMATH1, SBML.MathIdent("c1")])
sm = convert(Num, km)
kl = SBMLToolkit.getunidirectionalcomponents(sm)
@test isequal(kl, (k1, c1))

km = SBML.MathApply("-", SBML.Math[KINETICMATH1, KINETICMATH2])
sm = convert(Num, km)
fw, rv = SBMLToolkit.getunidirectionalcomponents(sm)
rv = substitute(rv, Dict(SBMLToolkit.create_var("s2") => SBMLToolkit.create_var("s2",Catalyst.DEFAULT_IV)))
@test isequal((fw, rv), (k1, k1*s2))

km = SBML.MathIdent("s1s2")
sm1 = convert(Num, km)
sm2 = sm - sm1
@test_throws ErrorException SBMLToolkit.getunidirectionalcomponents(sm2)
@test_throws ErrorException SBMLToolkit.getunidirectionalcomponents(:k1)

# Test get_paramap
trueparamap = [k1 => 1., c1 => 2.]
paramap = SBMLToolkit.get_paramap(MODEL1)
@test isequal(paramap, trueparamap)

# Test getmassaction
@test isequal(SBMLToolkit.getmassaction(k1*s1, [s1], [1]), k1)  # Case hOSU=true
@test isequal(SBMLToolkit.getmassaction(k1*s1/c1, [s1], [1]), k1/c1)  # Case hOSU=false
@test isequal(SBMLToolkit.getmassaction(k1+c1, nothing, nothing), k1+c1)  # Case zero order kineticLaw
@test isnan(SBMLToolkit.getmassaction(k1*s1*s2/(c1+s2), [s1], [1]))  # Case Michaelis-Menten kinetics
@test isnan(SBMLToolkit.getmassaction(k1*s1*Catalyst.DEFAULT_IV, [s1], [1]))  # Case kineticLaw with time

@test isnan(SBMLToolkit.getmassaction(k1*s1*s2, [s1], [1]))  
@test isnan(SBMLToolkit.getmassaction(k1+c1, [s1], [1]))
@test_throws ErrorException SBMLToolkit.getmassaction(k1, nothing, [1])

# default test
@test ModelingToolkit.defaults(m) == ModelingToolkit.defaults(ReactionSystem(m))
