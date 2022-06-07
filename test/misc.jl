# const case_ids = [7, 22, 140, 170, 679]
const case_ids = [200:300...]
const cases = map(x -> x[end-4:end], .*("0000", string.(case_ids)))

const algomap = Dict("00177" => Rodas4,
                  "00862" => Rodas4,
                  "00863" => Rodas4,
                  "00864" => Rodas4,
                  "00882" => Rodas4)

const ss_fail = ["00023", "00024"]
const logdir = joinpath(@__DIR__, "logs")
ispath(logdir) && rm(logdir,recursive=true)
mkdir(logdir)

const expected_errs = 
    ["Model contains no reactions.",
    "are not yet implemented.",
    "Please make reaction irreversible or rearrange kineticLaw to the form `term1 - term2`.",
    "BoundsError(String[], (1,))",  # Occurs wher no V3L2 file is available
    "COBREXA.jl",
    "no method matching length(::Nothing)", "MethodError(iterate, (nothing,),", # Occurs for insance in case 00029, where S1(t) = 7 is the only eqn.
    "Stoichiometry must be a non-negative integer.",
    "NaN result for non-NaN input.",
    "RequestError("]

function setup_settings_txt(text)
    ls = split(text, "\n")
    spls = split.(ls, ": ")
    filter!(x->length(x) == 2, spls)
    Dict(map(x -> x[1] => Meta.parse(x[2]), spls))
end

function to_concentrations(sol, ml, res_df, ia)
    volumes = [1.]
    sol_df = DataFrame(sol)
    for sn in names(sol_df)[2:end]
        if haskey(ml.species, sn[1:3-end])
            spec = ml.species[sn[1:end-3]]
            comp = ml.compartments[spec.compartment]
            ic = spec.initial_concentration
            isnothing(ic) || haskey(ia, sn[1:end-3]) ? push!(volumes, 1.) : push!(volumes, comp.size)
        else
            push!(volumes, 1.)
        end
    end
    sol_df = sol_df./Array(volumes)'
    
    idx = [sol.t[i] in res_df[:, 1] ? true : false for i in 1:length(sol.t)]
    sol_df = sol_df[idx, :]
    rename!(sol_df, "timestamp" => "time")
    rename!(sol_df, [rstrip(n, ['(', 't', ')']) for n in names(sol_df)])
    sol_df[:, [c for c in names(res_df) if c in names(sol_df)]]
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
    res_df = CSV.read(IOBuffer(results), DataFrame)
    (sbml, settings, res_df)
end

function verify_case(case; verbose=true)
    k = 0
    diffeq_retcode = :nothing
    expected_err = false
    res = false
    err = ""
    try
    # Read case
        sbml, settings, res_df = read_case(case)
    
        # Read SBML
        SBMLToolkit.checksupport(sbml)
        ml = readSBMLFromString(sbml, doc -> begin
            set_level_and_version(3, 2)(doc)
            convert_simplify_math(doc)
        end)
        ia = readSBMLFromString(sbml, doc -> begin
            set_level_and_version(3, 2)(doc)
        end)
        ia = ia.initial_assignments
        k = 1

        rs = ReactionSystem(ml)
        k = 2
        
        sys = convert(ODESystem, rs; include_zero_odes = true, combinatoric_ratelaws=false)
        if length(ml.events) > 0
            sys = ODESystem(ml)
        end
        k = 3
        
        # ssys = structural_simplify(sys)
        k = 4
        
        ts = res_df[:, 1]  # LinRange(settings["start"], settings["duration"], settings["steps"]+1)
        prob = ODEProblem(sys, Pair[], (settings["start"], Float64(settings["duration"])); saveat=ts)
        k = 5
        
        algorithm = get(algomap, case, Tsit5())
        sol = solve(prob, algorithm; abstol=settings["absolute"], reltol=settings["relative"])
        diffeq_retcode = sol.retcode
        k = diffeq_retcode == :Success ? 6 : k

        sol_df = to_concentrations(sol, ml, res_df, ia)
        CSV.write(joinpath(logdir, "SBMLTk_"*case*".csv"), sol_df)

        solm = Matrix(sol_df)
        resm = Matrix(res_df[:, [c for c in names(sol_df) if c in names(res_df)]])
        res = isapprox(solm, resm; atol=1e-9, rtol=3e-2)
        res || verify_plot(case, sys, solm, resm, ts)
    catch e
        err = string(e)
        expected_err = any(occursin.(expected_errs, err))
        if length(err) > 1000 # cutoff since I got ArgumentError: row size (9088174) too large 
            err = err[1:1000]
        end
    finally
        verbose && @info("Case $(case) done with a code $k and error msg: $err")
        return (case, expected_err, res, err, k, diffeq_retcode)
    end
end

function verify_all(cases; verbose=true)
    df = DataFrame(case=String[], expected_err=Bool[], res=Bool[],
                   error=String[], k=Int64[], diffeq_retcode=Symbol[])
    for case in cases
        ret = verify_case(case; verbose=verbose)
        verbose && @info ret 
        push!(df, ret)
    end
    verbose && print(df)
    fn = joinpath(logdir, "test_suite_$(cases[1])-$(cases[end]).csv")
    CSV.write(fn, df)
    df
end

df = verify_all(cases)
# sys = xx[]
# for i in 1:length(cases)
#     @test sum(Vector(df[i, ["expected_err", "res"]])) == 1
# end
