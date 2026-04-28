# this file tests that a big model simulates
# SBML.Model with 4139 reactions, 1265 species, and 522 parameters. (1012 equations)

using SBMLToolkit
using Downloads, ModelingToolkit, OrdinaryDiffEq
using Test

sbml_url = "https://www.biomodels.org/biomodels/services/download/get-files/MODEL1112100000/2/MODEL1112100000_url.xml"
sbml = try
    String(take!(Downloads.download(sbml_url, IOBuffer())))
catch e
    # BioModels' download endpoints intermittently 403 self-hosted CI runner
    # IPs (Cloudflare/WAF). The test exercises SBMLToolkit on a large model,
    # not the network, so skip rather than fail when the fetch is blocked.
    @warn "Wuschel test: SBML download failed, skipping" exception=(e, catch_backtrace())
    nothing
end

if sbml !== nothing
    m = readSBMLFromString(
        sbml, doc -> begin
            # set_level_and_version(3, 2)(doc) # fails on wuschel
            convert_promotelocals_expandfuns(doc)
        end
    )
    sys = ODESystem(m)
    @test length(equations(sys)) == 1012
    @test length(unknowns(sys)) == 1012
    #ssys = structural_simplify(sys)  # Todo: Figure out why this complains about ExtraVariablesSystemException
    prob = ODEProblem(sys, [], (0, 10.0))
    solve(prob, Tsit5(), save_everystep = false)
else
    @test_skip false
end
