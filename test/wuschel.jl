# this file tests that a big model simulates
# SBML.Model with 4139 reactions, 1265 species, and 522 parameters. (1012 equations)

myread(s) = readSBMLFromString(s, doc -> begin
    # set_level_and_version(3, 2)(doc) # fails on wuschel
    convert_simplify_math(doc)
end)

url = "https://www.ebi.ac.uk/biomodels/model/download/MODEL1112100000.2?filename=MODEL1112100000_url.xml"
io = IOBuffer()
Downloads.download(url, io)
s = String(take!(io))
m = myread(s)

rs = ReactionSystem(m)
@test isequal(length(rs.eqs), 4139)
@test isequal(length(rs.species), 1265)

ssys = structural_simplify(ODESystem(m))
prob = ODEProblem(ssys, [], (0, 10.0))
solve(prob, Tsit5(), save_everystep = false)
