using SBMLToolkit
using Catalyst, SBML, SBMLToolkitTestSuite
using Test

function readmodel(sbml)
    SBMLToolkit.readSBMLFromString(sbml, doc -> begin
                                       set_level_and_version(3, 2)(doc)
                                       convert_promotelocals_expandfuns(doc)
                                   end)
end

COMP1 = SBML.Compartment("C", true, 3, 2.0, "nl", nothing, nothing, nothing, nothing,
                         SBML.CVTerm[])
SPECIES1 = SBML.Species(name = "B", compartment = "C", initial_amount = 1.0,
                        substance_units = "substance", only_substance_units = true,
                        boundary_condition = false, constant = false)  # Todo: Maybe not support units in initial_concentration?
PARAM1 = SBML.Parameter(name = "D", value = 1.0, constant = true)
SPECIES2 = SBML.Species(name = "Bc", compartment = "C", initial_amount = 1.0,
                        substance_units = "substance", only_substance_units = true,
                        boundary_condition = false, constant = true)  # Todo: Maybe not support units in initial_concentration?
PARAM2 = SBML.Parameter(name = "Dv", value = 1.0, constant = false)
MODEL1 = SBML.Model(parameters = Dict("D" => PARAM1, "Dv" => PARAM2),
                    compartments = Dict("C" => COMP1),
                    species = Dict("B" => SPECIES1, "Bc" => SPECIES2)) 

const IV = Catalyst.DEFAULT_IV
@species s1(IV)
# Test symbolicsRateOf
rate = SBMLToolkit.symbolicsRateOf(s1)
rate_true = SBMLToolkit.D(s1)
@test isequal(rate, rate_true)

# Test map_symbolics_ident
sym = SBMLToolkit.map_symbolics_ident(SBML.MathIdent("B"), MODEL1)
@species B(IV) [irreducible = false; isbcspecies = false]
@test isequal(sym, B)

# Test interpret_as_num
@species B(IV)
@variables Dv(IV)
@parameters C D Bc

test = SBML.MathApply("*",
                      SBML.Math[SBML.MathApply("+",
                                               SBML.Math[SBML.MathApply("*",
                                                                        SBML.Math[SBML.MathAvogadro("A"),
                                                                                  SBML.MathIdent("B")]),
                                                         SBML.MathApply("-",
                                                                        SBML.Math[SBML.MathApply("*",
                                                                                                 SBML.Math[SBML.MathIdent("C"),
                                                                                                           SBML.MathIdent("D"),
                                                                                                           SBML.MathIdent("Dv"),
                                                                                                           SBML.MathIdent("Bc")])])]),
                                SBML.MathTime("Time")])

@test isequal(SBMLToolkit.interpret_as_num(test, MODEL1), IV * (6.02214076e23 * B - C * D * Dv * Bc))

# Test get_substitutions
sbml, _, _ = SBMLToolkitTestSuite.read_case("00001")
m = readmodel(sbml)
subsdict = SBMLToolkit.get_substitutions(m)
@parameters k1, compartment
@species S1(IV), S2(IV)
subsdict_true = Dict(Num(Symbolics.variable(:S1; T = Real)) => S1,
                     Num(Symbolics.variable(:S2; T = Real)) => S2,
                     Num(Symbolics.variable(:k1; T = Real)) => k1,
                     Num(Symbolics.variable(:compartment; T = Real)) => compartment)
@test isequal(subsdict, subsdict_true)

# Test create_var
var = SBMLToolkit.create_var("s2")
@species s2
@test isequal(var, s2)

var = SBMLToolkit.create_var("s2", IV, isbcspecies = true)
Catalyst.isbc(var)

# Test create_param
par = SBMLToolkit.create_param("k")
@parameters k
@test isequal(par, k)

par = SBMLToolkit.create_param("k", isconstantspecies = true)
@test Catalyst.isconstant(par)

# Test has_rule_type
sbml, _, _ = SBMLToolkitTestSuite.read_case("00031")  # rateRule
m = readmodel(sbml)
res = SBMLToolkit.has_rule_type("S1", m, SBML.RateRule)
@test res

res = SBMLToolkit.has_rule_type("nospecies", m, SBML.RateRule)
@test !res

sbml, _, _ = SBMLToolkitTestSuite.read_case("00039")  # algebraicRule
m = readmodel(sbml)
res = SBMLToolkit.has_rule_type("S1", m, SBML.AlgebraicRule)
@test res

res = SBMLToolkit.has_rule_type("nospecies", m, SBML.AlgebraicRule)
@test !res
