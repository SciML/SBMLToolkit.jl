using SafeTestsets, Test

@safetestset "Quality Assurance" include("qa.jl")
@safetestset "Systems" include("systems.jl")
@safetestset "Reactions" include("reactions.jl")
@safetestset "Rules" include("rules.jl")
@safetestset "Events" include("events.jl")
@safetestset "Utils" include("utils.jl")
@safetestset "Simulation results" include("simresults.jl")
@safetestset "Wuschel" include("wuschel.jl")
