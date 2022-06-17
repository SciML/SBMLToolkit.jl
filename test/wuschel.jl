# this file tests that a big model simulates
# SBML.Model with 4139 reactions, 1265 species, and 522 parameters. (1012 equations)

sbml_url = "https://www.ebi.ac.uk/biomodels/model/download/MODEL1112100000.2?filename=MODEL1112100000_url.xml"
sbml = String(take!(Downloads.download(sbml_url, IOBuffer())))
m = readSBMLFromString(sbml, doc -> begin
        # set_level_and_version(3, 2)(doc) # fails on wuschel
        convert_simplify_math(doc)
    end)
sys = ODESystem(m)
@test length(equations(sys)) == 1012
@test length(states(sys)) == 1012
#ssys = structural_simplify(sys)  # Todo: Figure out why this complains about ExtraVariablesSystemException
prob = ODEProblem(sys, [], (0, 10.0))
solve(prob, Tsit5(), save_everystep = false)
