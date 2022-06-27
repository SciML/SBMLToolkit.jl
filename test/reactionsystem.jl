using SBMLToolkit
using Catalyst, SBML, ModelingToolkit
using OrdinaryDiffEq
using Test

cd(@__DIR__)
sbmlfile = joinpath("data", "reactionsystem_01.xml")
@parameters t, k1, c1
@variables s1(t), s2(t), s1s2(t), s3(t)
const IV = Catalyst.DEFAULT_IV

kinetic_params = Dict{String,Tuple{Float64,String}}()

COMP1 = SBML.Compartment("c1", true, 3, 2.0, "nl", nothing, nothing)
SPECIES1 = SBML.Species(name = "s1", compartment = "c1", initial_amount = 1.0, substance_units = "substance", only_substance_units = true, boundary_condition = false, constant = false)  # Todo: Maybe not support units in initial_concentration?
SPECIES2 = SBML.Species(name = "s2", compartment = "c1", initial_amount = 1.0, substance_units = "substance/nl", only_substance_units = false)
SPECIES3 = SBML.Species(name = "s3", compartment = "c1", initial_amount = 1.0, substance_units = "substance", only_substance_units = false, constant=true)
KINETICMATH1 = SBML.MathIdent("k1")
KINETICMATH2 = SBML.MathApply("*", SBML.Math[
    SBML.MathIdent("k1"), SBML.MathIdent("s2")])
KINETICMATH3 = SBML.MathApply("-",
    SBML.Math[SBML.MathApply("*", SBML.Math[
        SBML.MathIdent("k1"), SBML.MathIdent("s1")]),
    KINETICMATH1])
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
    reactants = Dict("s1" => 1),
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
    species = Dict("s1" => SPECIES1),
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
@test isequal(Catalyst.get_eqs(rs), Catalyst.Reaction[Catalyst.Reaction(k1 / c1, [s1, s2], [s1s2], [1.0, 1.0], [1.0])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s1, s1s2, s2])
@test isequal(Catalyst.get_ps(rs), [k1, c1])
@named rs = ReactionSystem(MODEL1)
isequal(nameof(rs), :rs)

rs = ReactionSystem(MODEL3)  # Contains reversible reaction
@test isequal(Catalyst.get_eqs(rs),
    Catalyst.Reaction[
        Catalyst.Reaction(k1, [s1], nothing,
            [1], nothing),
        Catalyst.Reaction(k1, nothing, [s1],
            nothing, [1])])
@test isequal(Catalyst.get_iv(rs), t)
@test isequal(Catalyst.get_states(rs), [s1])
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
trueeqs = Equation[Differential(t)(s1)~-((k1*s1*s2) / c1),
    Differential(t)(s1s2)~(k1*s1*s2) / c1,
    Differential(t)(s2)~-((k1*s1*s2) / c1)]
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
@test_nowarn SBMLToolkit.checksupport_file(sbmlfile)
@test_throws ErrorException SBMLToolkit.checksupport_file(joinpath("data", "unsupported.sbml"))
@test_nowarn SBMLToolkit.checksupport_string("all good <reaction")
@test_throws ErrorException SBMLToolkit.checksupport_string("contains </delay>")

# Test get_substitutions
truesubs = Dict(Num(Symbolics.variable(:c1; T = Real)) => c1,
    Num(Symbolics.variable(:k1; T = Real)) => k1,
    Num(Symbolics.variable(:s1; T = Real)) => s1)
subs = SBMLToolkit.get_substitutions(MODEL1)
@test isequal(subs, truesubs)

# Test get_reactions
reaction = SBMLToolkit.get_reactions(MODEL1)[1]
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
@test isequal(IV, SBMLToolkit.get_reactions(mod)[1].rate)

# Test use_rate
@test isequal(SBMLToolkit.use_rate(k1 * s1, [s1], [1]), (k1, false))  # Case hOSU=true
@test isequal(SBMLToolkit.use_rate(k1 * s1 * s2 / (c1 + s2), [s1], [1]), (k1 * s1 * s2 / (c1 + s2), true))  # Case Michaelis-Menten kinetics

# Test get_reagents
@test isequal((nothing, [s1], nothing, [1.0]), SBMLToolkit.get_reagents(REACTION1.reactants, REACTION1.products, MODEL1))

s = SBML.Species(name = "s", compartment = "c1", boundary_condition = true,
                 initial_amount = 1.0, substance_units = "substance",
                 only_substance_units = true)
var = SBMLToolkit.create_var("s", IV)
r = SBML.Reaction(reactants = Dict("s" => 1), products = Dict(),
                  reversible=false)

m = SBML.Model(species = Dict("s" => s), reactions = Dict("r1" => r))
@test isequal(([var], [var], [1.0], [1.0]), SBMLToolkit.get_reagents(r.reactants, r.products, m))
@test isequal((nothing, nothing, nothing, nothing), SBMLToolkit.get_reagents(r.products, r.reactants, m))

# Test get_unidirectional_components
km = SBML.MathApply("-", SBML.Math[KINETICMATH1, SBML.MathIdent("c1")])
sm = SBMLToolkit.interpret_as_num(km)
kl = SBMLToolkit.get_unidirectional_components(sm)
@test isequal(kl, (k1, c1))

km = SBML.MathApply("-", SBML.Math[KINETICMATH1, KINETICMATH2])
sm = SBMLToolkit.interpret_as_num(km)
fw, rv = SBMLToolkit.get_unidirectional_components(sm)
rv = substitute(rv, Dict(SBMLToolkit.create_var("s2") => SBMLToolkit.create_var("s2", IV)))
@test isequal((fw, rv), (k1, k1 * s2))

km = SBML.MathIdent("s1s2")
sm1 = SBMLToolkit.interpret_as_num(km)
sm2 = sm - sm1
@test isequal(SBMLToolkit.get_unidirectional_components(sm2), (sm2, Num(0)))
@test isequal(SBMLToolkit.get_unidirectional_components(:k1), (:k1, Num(0)))

# Test get_mappings
u0map, parammap = SBMLToolkit.get_mappings(MODEL1)
u0map_true = [s1 => 1.0]
parammap_true = [k1 => 1.0, c1 => 2.0]
@test isequal(u0map, u0map_true)
@test isequal(parammap, parammap_true)

p = SBML.Parameter(name = "k2", value = 1.0, constant = false)
m = SBML.Model(
    parameters = Dict("k2" => p),
    rules = SBML.Rule[SBML.RateRule("k2", KINETICMATH2)])
u0map, parammap = SBMLToolkit.get_mappings(m)
u0map_true = [SBMLToolkit.create_var("k2", IV; isbcspecies=true) => 1.0]
@test isequal(u0map, u0map_true)
@test Catalyst.isbc(first(u0map[1]))

m = SBML.Model(
    species = Dict("s2" => SPECIES2),
    rules = SBML.Rule[SBML.AlgebraicRule(KINETICMATH2)])
u0map, parammap = SBMLToolkit.get_mappings(m)
Catalyst.isbc(first(u0map[1]))

# Test get_massaction
@test isequal(SBMLToolkit.get_massaction(k1 * s1, [s1], [1]), k1)  # Case hOSU=true
@test isequal(SBMLToolkit.get_massaction(k1 * s1 / c1, [s1], [1]), k1 / c1)  # Case hOSU=false
@test isequal(SBMLToolkit.get_massaction(k1 + c1, nothing, nothing), k1 + c1)  # Case zero order kineticLaw
@test isnan(SBMLToolkit.get_massaction(k1 * s1 * s2 / (c1 + s2), [s1], [1]))  # Case Michaelis-Menten kinetics
@test isnan(SBMLToolkit.get_massaction(k1 * s1 * IV, [s1], [1]))  # Case kineticLaw with time

@test isnan(SBMLToolkit.get_massaction(k1 * s1 * s2, [s1], [1]))
@test isnan(SBMLToolkit.get_massaction(k1 + c1, [s1], [1]))
@test_throws ErrorException SBMLToolkit.get_massaction(k1, nothing, [1])
