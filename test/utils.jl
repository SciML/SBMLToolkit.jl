using SBMLToolkit
using SBML, SBMLToolkitTestSuite
using Test

function readmodel(sbml)
    SBMLToolkit.readSBMLFromString(sbml, doc -> begin
                                       set_level_and_version(3, 2)(doc)
                                       convert_simplify_math(doc)
                                   end)
end

const IV = Catalyst.DEFAULT_IV
@variables s1(IV)
# Test symbolicsRateOf
rate = SBMLToolkit.symbolicsRateOf(s1)
rate_true = SBMLToolkit.D(s1)
@test isequal(rate, rate_true)

# Test map_sumbolics_ident
sym = SBMLToolkit.map_symbolics_ident(SBML.MathIdent("s1"))
@variables s1
@test isequal(sym, s1)

# Test interpret_as_num
@variables A B C D Time

test = SBML.MathApply("*",
                      SBML.Math[SBML.MathApply("+",
                                               SBML.Math[SBML.MathApply("*",
                                                                        SBML.Math[SBML.MathIdent("A"),
                                                                                  SBML.MathIdent("B")]),
                                                         SBML.MathApply("-",
                                                                        SBML.Math[SBML.MathApply("*",
                                                                                                 SBML.Math[SBML.MathIdent("C"),
                                                                                                           SBML.MathIdent("D")])])]),
                                SBML.MathTime("Time")])

@test isequal(SBMLToolkit.interpret_as_num(test), IV * (A * B - C * D))

# Test get_substitutions
sbml, _, _ = SBMLToolkitTestSuite.read_case("00001")
m = readmodel(sbml)
subsdict = SBMLToolkit.get_substitutions(m)
@parameters k1, compartment
@variables S1(IV), S2(IV)
subsdict_true = Dict(Num(Symbolics.variable(:S1; T = Real)) => S1,
                     Num(Symbolics.variable(:S2; T = Real)) => S2,
                     Num(Symbolics.variable(:k1; T = Real)) => k1,
                     Num(Symbolics.variable(:compartment; T = Real)) => compartment)
@test isequal(subsdict, subsdict_true)

# Test create_var
var = SBMLToolkit.create_var("s2")
@variables s2
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
