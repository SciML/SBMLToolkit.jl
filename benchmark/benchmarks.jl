using SBMLToolkit, BenchmarkTools

const SUITE = BenchmarkGroup()

data_dir = joinpath(@__DIR__, "..", "test", "data")

# =============================================================================
# SBML reading → ReactionSystem/ODESystem
# =============================================================================

SUITE["read"] = BenchmarkGroup()

SUITE["read"]["model_26"] = @benchmarkable readSBML(
    joinpath($data_dir, "00026-sbml-l3v2.xml")
)
SUITE["read"]["model_31"] = @benchmarkable readSBML(
    joinpath($data_dir, "00031-sbml-l3v2.xml")
)
SUITE["read"]["model_31_odesys"] = @benchmarkable readSBML(
    joinpath($data_dir, "00031-sbml-l3v2.xml"), ODESystemImporter()
)
SUITE["read"]["model_31_rssys"] = @benchmarkable readSBML(
    joinpath($data_dir, "00031-sbml-l3v2.xml"), ReactionSystemImporter()
)
