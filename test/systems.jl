using SBMLToolkit
using Catalyst, SBML
using Test

cd(@__DIR__)
sbmlfile = joinpath("data", "reactionsystem_01.xml")
const IV = default_t()
@parameters k1, c1
@species s1(IV), s2(IV), s1s2(IV)

COMP1 = SBML.Compartment("c1", true, 3, 2.0, "nl", nothing, nothing, nothing, nothing,
    SBML.CVTerm[])
SPECIES1 = SBML.Species(name = "s1", compartment = "c1", initial_amount = 1.0,
    substance_units = "substance", only_substance_units = true,
    boundary_condition = false, constant = false)  # Todo: Maybe not support units in initial_concentration?
SPECIES2 = SBML.Species(name = "s2", compartment = "c1", initial_amount = 1.0,
    substance_units = "substance/nl", only_substance_units = false)
KINETICMATH1 = SBML.MathIdent("k1")
KINETICMATH2 = SBML.MathApply("*", SBML.Math[SBML.MathIdent("k1"), SBML.MathIdent("s2")])
KINETICMATH3 = SBML.MathApply("-",
    SBML.Math[SBML.MathApply("*",
            SBML.Math[SBML.MathIdent("k1"),
                SBML.MathIdent("s1")]),
        KINETICMATH1])
REACTION1 = SBML.Reaction(
    products = [
        SBML.SpeciesReference(species = "s1", stoichiometry = 1.0)
    ],
    kinetic_math = KINETICMATH1,
    reversible = false)
REACTION2 = SBML.Reaction(
    reactants = [
        SBML.SpeciesReference(species = "s1", stoichiometry = 1.0)
    ],
    kinetic_math = KINETICMATH3,
    reversible = true)
PARAM1 = SBML.Parameter(name = "k1", value = 1.0, constant = true)
MODEL1 = SBML.Model(parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s1" => SPECIES1),
    reactions = Dict("r1" => REACTION1))  # PL: For instance in the compartments dict, we may want to enforce that key and compartment.name are identical
MODEL2 = SBML.Model(parameters = Dict("k1" => PARAM1),
    compartments = Dict("c1" => COMP1),
    species = Dict("s1" => SPECIES1),
    reactions = Dict("r3" => REACTION2))

# Test ReactionSystem constructor
rs = ReactionSystem(MODEL1)
@test isequal(Catalyst.get_eqs(rs),
    Catalyst.Reaction[Catalyst.Reaction(k1, nothing, [s1], nothing, [1.0])])
@test isequal(Catalyst.get_iv(rs), IV)
@test isequal(Catalyst.get_species(rs), [s1])
@test issetequal(Catalyst.get_ps(rs), [k1, c1])
@named rs = ReactionSystem(MODEL1)
isequal(nameof(rs), :rs)

rs = ReactionSystem(readSBML(sbmlfile))
@test isequal(Catalyst.get_eqs(rs),
    Catalyst.Reaction[Catalyst.Reaction(k1 / c1, [s1, s2], [s1s2], [1.0, 1.0],
        [1.0])])
@test isequal(Catalyst.get_iv(rs), IV)
@test isequal(Catalyst.get_species(rs), [s1, s1s2, s2])
@test issetequal(Catalyst.get_ps(rs), [k1, c1])
@named rs = ReactionSystem(MODEL1)
isequal(nameof(rs), :rs)

rs = ReactionSystem(MODEL2)  # Contains reversible reaction
@test isequal(Catalyst.get_eqs(rs),
    Catalyst.Reaction[Catalyst.Reaction(k1, [s1], nothing,
            [1], nothing),
        Catalyst.Reaction(k1, nothing, [s1],
            nothing, [1])])
@test isequal(Catalyst.get_iv(rs), IV)
@test isequal(Catalyst.get_species(rs), [s1])
@test issetequal(Catalyst.get_ps(rs), [k1, c1])

@test convert(ModelingToolkit.ODESystem, rs) isa ODESystem
@test structural_simplify(convert(ModelingToolkit.ODESystem, rs)) isa ODESystem

# Test ODESystem constructor
odesys = ODESystem(MODEL1)
trueeqs = Equation[default_time_deriv()(s1) ~ k1]
@test isequal(Catalyst.get_eqs(odesys), trueeqs)
@test isequal(Catalyst.get_iv(odesys), IV)
@test isequal(Catalyst.get_unknowns(odesys), [s1])
@test issetequal(Catalyst.get_ps(odesys), [k1, c1])
u0 = [s1 => 1.0]
par = [k1 => 1.0, c1 => 2.0]
@test isequal(ModelingToolkit.defaults(odesys), ModelingToolkit._merge(u0, par))  # PL: @Anand: for some reason this does not work with `Catalyst.get_default()`
@named odesys = ODESystem(MODEL1)
isequal(nameof(odesys), :odesys)
@test structural_simplify(odesys) isa ODESystem

odesys = ODESystem(readSBML(sbmlfile))
m = readSBML(sbmlfile)
trueeqs = Equation[default_time_deriv()(s1) ~ -((k1 * s1 * s2) / c1),
    default_time_deriv()(s1s2) ~ (k1 * s1 * s2) / c1,
    default_time_deriv()(s2) ~ -((k1 * s1 * s2) / c1)]
@test isequal(Catalyst.get_eqs(odesys), trueeqs)
@test isequal(Catalyst.get_iv(odesys), IV)
@test isequal(Catalyst.get_unknowns(odesys), [s1, s1s2, s2])
@test issetequal(Catalyst.get_ps(odesys), [k1, c1])
u0 = [s1 => 2 * 1.0, s2 => 2 * 1.0, s1s2 => 2 * 1.0]
par = [k1 => 1.0, c1 => 2.0]
@test isequal(ModelingToolkit.defaults(odesys), ModelingToolkit._merge(u0, par))
@named odesys = ODESystem(MODEL1)
isequal(nameof(odesys), :odesys)

@test ODEProblem(odesys, [], [0.0, 1.0], []) isa ODEProblem

# # Test ODEProblem
# oprob = ODEProblem(ODESystem(MODEL1), [], [0.0, 1.0], [])
# sol = solve(oprob, Tsit5())
# @test isapprox(sol.u, [[1.0], [2.0]])

# @test_nowarn ODEProblem(ODESystem(readSBML(sbmlfile)), [], [0.0, 1.0], [])

# Test get_mappings
u0map, parammap, initial_assignment_map = SBMLToolkit.get_mappings(MODEL1)
u0map_true = [s1 => 1.0]
parammap_true = [k1 => 1.0, c1 => 2.0]
initial_assignment_map_true = Pair[]
@test isequal(u0map, u0map_true)
@test isequal(parammap, parammap_true)
@test isequal(initial_assignment_map, initial_assignment_map_true)

p = SBML.Parameter(name = "k2", value = 1.0, constant = false)
m = SBML.Model(parameters = Dict("k2" => p),
    rules = SBML.Rule[SBML.RateRule("k2", KINETICMATH2)])
u0map, parammap, initial_assignment_map = SBMLToolkit.get_mappings(m)
u0map_true = [SBMLToolkit.create_var("k2", IV; isbcspecies = true) => 1.0]
@test isequal(u0map, u0map_true)
@test Catalyst.isbc(first(u0map[1]))

p = SBML.Parameter(name = "k2", value = nothing, constant = true)
ia = Dict("k2" => KINETICMATH1)
m = SBML.Model(parameters = Dict("k1" => PARAM1, "k2" => p),
    initial_assignments = ia)
u0map, parammap, initial_assignment_map = SBMLToolkit.get_mappings(m)
parammap_true = [k1 => 1.0, SBMLToolkit.create_var("k2") => k1]
initial_assignment_map_true = [SBMLToolkit.create_var("k2") => k1]
@test isequal(parammap, parammap_true)
@test isequal(initial_assignment_map, initial_assignment_map_true)

m = SBML.Model(species = Dict("s2" => SPECIES2),
    rules = SBML.Rule[SBML.AlgebraicRule(KINETICMATH2)])
u0map, parammap, initial_assignment_map = SBMLToolkit.get_mappings(m)
Catalyst.isbc(first(u0map[1]))

# Test get_netstoich
r = SBML.Reaction(reactants = [SBML.SpeciesReference(species = "s1", stoichiometry = 1.0)],
    kinetic_math = KINETICMATH1,
    reversible = false)
ns = SBMLToolkit.netstoich("s1", r)
@test isequal(ns, -1)

r = SBML.Reaction(reactants = [SBML.SpeciesReference(species = "s1", stoichiometry = 1.0)],
    products = [SBML.SpeciesReference(species = "s1", stoichiometry = 2.0)],
    kinetic_math = KINETICMATH1,
    reversible = false)
ns = SBMLToolkit.netstoich("s1", r)
@test isequal(ns, 1)

# Test checksupport_file
@test_nowarn SBMLToolkit.checksupport_file(sbmlfile)
@test_throws ErrorException SBMLToolkit.checksupport_file(joinpath("data",
    "unsupported.sbml"))

# Test checksupport_string
@test_nowarn SBMLToolkit.checksupport_string("all good <reaction")
@test_throws ErrorException SBMLToolkit.checksupport_string("contains </delay>")

# Test convenience functions
@test readSBML(sbmlfile, DefaultImporter()) isa SBML.Model
@test readSBML(sbmlfile, ReactionSystemImporter()) isa ReactionSystem
@test readSBML(sbmlfile, ODESystemImporter()) isa ODESystem
