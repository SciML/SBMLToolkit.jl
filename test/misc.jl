const cases = ["00007", "00022", "00140"]

const algo = Dict("00862" => Rodas4,
                  "00863" => Rodas4,
                  "00864" => Rodas4,
                  "00882" => Rodas4)
                  
const logdir = joinpath(@__DIR__, "logs")
ispath(logdir) && rm(logdir,recursive=true)
mkdir(logdir)

function setup_settings_txt(text)
    ls = split(text, "\n")
    spls = split.(ls, ": ")
    filter!(x->length(x) == 2, spls)
    Dict(map(x -> x[1] => Meta.parse(x[2]), spls))
end

function to_concentrations(sol, ml, results)
    volumes = [1.]
    sol_df = DataFrame(sol)
    for sn in names(sol_df)[2:end]
        if haskey(ml.species, sn[1:3-end])
            spec = ml.species[sn[1:end-3]]
            comp = ml.compartments[spec.compartment]
            ic = spec.initial_concentration
            isnothing(ic) ? push!(volumes, 1.) : push!(volumes, comp.size)
        end
    end
    sol_df = sol_df./Array(volumes)'
    
    idx = [sol.t[i] in results[:, 1] ? true : false for i in 1:length(sol.t)]
    sol_df = sol_df[idx, :]
    sol_df
end

"plots the difference between the suites' reported solution and DiffEq's sol"
function verify_plot(case, sys, solm, resm, ts)
    open(joinpath(logdir, case*".txt"), "w") do file
        write(file, "Reactions:\n")
        write(file, repr(equations(sys))*"\n")
        write(file, "ODEs:\n")
        write(file, repr(equations(sys))*"\n")
    end
    plt = plot(ts, solm)
    plt = plot!(ts, resm, linestyle=:dot)
    savefig(joinpath(logdir, case*".png"))
end

function read_case(case)
    base_url = "https://raw.githubusercontent.com/sbmlteam/sbml-test-suite/master/cases/semantic/$case/$case"
    sbml_url = base_url*"-sbml-l3v2.xml"
    settings_url = base_url*"-settings.txt"
    results_url = base_url*"-results.csv"
    
    sbml = String(take!(Downloads.download(sbml_url, IOBuffer())))
    settings = String(take!(Downloads.download(settings_url, IOBuffer())))
    results = String(take!(Downloads.download(results_url, IOBuffer())))
    
    # Read results
    settings = setup_settings_txt(settings)
    results = CSV.read(IOBuffer(results), DataFrame)
    (sbml, settings, results)
end

function verify(case)
    # Read case
    sbml, settings, results = read_case(case)
   
    # Read SBML
    ml = readSBMLFromString(sbml, doc -> begin
                set_level_and_version(3, 2)(doc)
                convert_simplify_math(doc)
            end)

    rs = ReactionSystem(ml)

    sys = convert(ODESystem, rs; include_zero_odes = false)
    if length(ml.events) > 0
        sys = ODESystem(ml)
    end
    
    ssys = structural_simplify(sys)
    
    ts = results[:, 1]  # LinRange(settings["start"], settings["duration"], settings["steps"]+1)
    prob = ODEProblem(ssys, Pair[], (settings["start"], Float64(settings["duration"])); saveat=ts, check_length=false)

    algorithm = case in keys(algo) ? algo[case] : Tsit5()
    sol = solve(prob, algorithm; abstol=settings["absolute"], reltol=settings["relative"], saveat=ts)
    sol_df = to_concentrations(sol, ml, results)
    CSV.write(joinpath(logdir, "SBMLTk_"*case*".csv"), sol_df)

    cols = names(sol_df)[2:end]
    res_df = results[:, [c[1:end-3] for c in cols]]
    solm = Matrix(sol_df[:, cols])
    resm = Matrix(res_df)
    res = isapprox(solm, resm; atol=1e-9, rtol=3e-2)
    @test res
    atol = maximum(solm .- resm)
    verify_plot(case, sys, solm, resm, ts)
    nothing
end

for case in cases
    verify(case)
end
