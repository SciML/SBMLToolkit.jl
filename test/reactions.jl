using SBMLToolkit
using Catalyst, SBML, ModelingToolkit
using Test

cd(@__DIR__)
sbmlfile = joinpath("data", "reactionsystem_01.xml")
const IV = Catalyst.DEFAULT_IV
@parameters k1, c1
@species s1(IV), s2(IV), s1s2(IV)

COMP1 = SBML.Compartment("c1", true, 3, 2.0, "nl", nothing, nothing, nothing, nothing,
                         SBML.CVTerm[])
SPECIES1 = SBML.Species(name = "s1", compartment = "c1", initial_amount = 1.0,
                        substance_units = "substance", only_substance_units = true,
                        boundary_condition = false, constant = false)  # Todo: Maybe not support units in initial_concentration?
KINETICMATH1 = SBML.MathIdent("k1")
KINETICMATH2 = SBML.MathApply("*", SBML.Math[SBML.MathIdent("k1"), SBML.MathIdent("s2")])
REACTION1 = SBML.Reaction(products = [
                              SBML.SpeciesReference(species = "s1", stoichiometry = 1),
                          ],
                          kinetic_math = KINETICMATH1,
                          reversible = false)
PARAM1 = SBML.Parameter(name = "k1", value = 1.0, constant = true)
MODEL1 = SBML.Model(parameters = Dict("k1" => PARAM1),
                    compartments = Dict("c1" => COMP1),
                    species = Dict("s1" => SPECIES1),
                    reactions = Dict("r1" => REACTION1))  # PL: For instance in the compartments dict, we may want to enforce that key and compartment.name are identical

# Test get_reactions
reaction = SBMLToolkit.get_reactions(MODEL1)[1]
truereaction = Catalyst.Reaction(k1, nothing, [s1], nothing, [1])  # Todo: implement Sam's suggestion on mass action kinetics
@test isequal(reaction, truereaction)

km = SBML.MathTime("x")
reac = SBML.Reaction(reactants = [
                         SBML.SpeciesReference(species = "s1", stoichiometry = 1.0),
                     ],
                     kinetic_math = km,
                     reversible = false)
model = SBML.Model(parameters = Dict("k1" => PARAM1),
                   compartments = Dict("c1" => COMP1),
                   species = Dict("s1" => SPECIES1),
                   reactions = Dict("r1" => reac))
@test isequal(IV, SBMLToolkit.get_reactions(model)[1].rate)

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
@test isequal(SBMLToolkit.get_unidirectional_components(k1), (k1, Num(0)))

# Test add_reaction!
rxs = Catalyst.Reaction[]
SBMLToolkit.add_reaction!(rxs, k1, SBML.SpeciesReference[], SBML.SpeciesReference[], MODEL1)
@test isequal(rxs, Catalyst.Reaction[])

rxs = Catalyst.Reaction[]
SBMLToolkit.add_reaction!(rxs, k1 * s1,
                          [SBML.SpeciesReference(species = "s1", stoichiometry = 1.0)],
                          SBML.SpeciesReference[], MODEL1)
reaction_true = Catalyst.Reaction(k1, [s1], nothing, [1], nothing, only_use_rate = false)
@test isequal(rxs[1], reaction_true)

rxs = Catalyst.Reaction[]
SBMLToolkit.add_reaction!(rxs, k1 * s1,
                          [SBML.SpeciesReference(species = "s1", stoichiometry = 1.0)],
                          SBML.SpeciesReference[], MODEL1,
                          enforce_rate = true)
reaction_true = Catalyst.Reaction(k1, [s1], nothing, [1], nothing, only_use_rate = true)
@test isequal(rxs, [reaction_true])

# Test stoich_convert_to_ints
@test isnothing(SBMLToolkit.stoich_convert_to_ints(nothing))
t = typeof(SBMLToolkit.stoich_convert_to_ints([1.0])[1])
@test isequal(t, Int64)

t = typeof(SBMLToolkit.stoich_convert_to_ints([1.1, 1])[1])
@test isequal(t, Float64)

# Test get_reagents
@test isequal((nothing, [s1], nothing, [1.0]),
              SBMLToolkit.get_reagents(REACTION1.reactants, REACTION1.products, MODEL1))

s = SBML.Species(name = "s", compartment = "c1", boundary_condition = true,
                 initial_amount = 1.0, substance_units = "substance",
                 only_substance_units = true)
var = SBMLToolkit.create_var("s", IV)
r = SBML.Reaction(reactants = [SBML.SpeciesReference(species = "s", stoichiometry = 1.0)],
                  reversible = false)

m = SBML.Model(species = Dict("s" => s), reactions = Dict("r1" => r))
@test isequal(([var], [var], [1.0], [1.0]),
              SBMLToolkit.get_reagents(r.reactants, r.products, m))
@test isequal((nothing, nothing, nothing, nothing),
              SBMLToolkit.get_reagents(r.products, r.reactants, m))

# Test use_rate
@test isequal(SBMLToolkit.use_rate(k1 * s1, [s1], [1]), (k1, false))  # Case hOSU=true
@test isequal(SBMLToolkit.use_rate(k1 * s1 * s2 / (c1 + s2), [s1], [1]),
              (k1 * s1 * s2 / (c1 + s2), true))  # Case Michaelis-Menten kinetics

# Test get_massaction
@test isequal(SBMLToolkit.get_massaction(k1 * s1, [s1], [1]), k1)  # Case hOSU=true
@test isequal(SBMLToolkit.get_massaction(k1 * s1 / c1, [s1], [1]), k1 / c1)  # Case hOSU=false
@test isequal(SBMLToolkit.get_massaction(k1 + c1, nothing, nothing), k1 + c1)  # Case zero order kineticLaw
@test isnan(SBMLToolkit.get_massaction(k1 * s1 * s2 / (c1 + s2), [s1], [1]))  # Case Michaelis-Menten kinetics
@test isnan(SBMLToolkit.get_massaction(k1 * s1 * IV, [s1], [1]))  # Case kineticLaw with time

@test isnan(SBMLToolkit.get_massaction(k1 * s1 * s2, [s1], [1]))
@test isnan(SBMLToolkit.get_massaction(k1 + c1, [s1], [1]))
@test_throws ErrorException SBMLToolkit.get_massaction(k1, nothing, [1])
