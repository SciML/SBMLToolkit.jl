using Test
using SBMLToolkit

using DifferentialEquations
using CSV, DataFrames
using Plots

function setup_settings_txt(fn)
    ls = readlines(fn)
    spls = split.(ls, ": ")
    filter!(x->length(x) == 2, spls)
    Dict(map(x -> x[1] => Meta.parse(x[2]), spls))
end

casedir = joinpath("..","SBMLBioModelsRepository.jl","data","sbml-test-suite", "semantic")
case = "00191"  # "00056" for reversible
url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$case/$case-sbml-l3v2.xml"
sbml_fn = joinpath(casedir,case,case*"-sbml-l3v2.xml")
settings_fn = joinpath(casedir,case,case*"-settings.txt")
results_fn = joinpath(casedir,case,case*"-results.csv")

if !isfile(sbml_fn)
    download(url, sbml_fn)
end
ml = readSBML(sbml_fn, doc -> begin
            set_level_and_version(3, 2)(doc)
            convert_simplify_math(doc)
        end)
# println("size")
# println(ml.compartments["compartment"].size)
rs = ReactionSystem(ml)
println("rate")
println(rs.eqs)
sys = ODESystem(ml)
n_dvs = length(states(sys))
n_ps = length(parameters(sys))
# structural_simplify(sys)
settings = setup_settings_txt(settings_fn)
println(settings)
results = CSV.read(results_fn, DataFrame)
ts = LinRange(settings["start"], settings["duration"], settings["steps"])
prob = ODEProblem(sys, Pair[], (settings["start"], Float64(settings["duration"])); saveat=ts)
println(prob)
sol = solve(prob, Tsit5())
solm = Array(sol)'
m = Matrix(results[1:end-1, 2:end])

# pyplot()
plt = plot(solm, linestyle=:dot)
plt = plot!(m)
gui(plt)
rm(joinpath("test","logs",),recursive=true)
mkdir(joinpath("test","logs"))
savefig(joinpath("test","logs",case*".png"))

# @test isapprox(solm, m; atol=1) ? true : false
# @test isapprox(solm, m; rtol=5e-2) ? true : false
