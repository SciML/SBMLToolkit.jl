using Pkg

using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "QA"
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include("qa/qa.jl")
else
    @safetestset "SBMLToolkit.jl" begin
        @safetestset "Systems" begin
            include("systems.jl")
        end
        @safetestset "Reactions" begin
            include("reactions.jl")
        end
        @safetestset "Rules" begin
            include("rules.jl")
        end
        @safetestset "Events" begin
            include("events.jl")
        end
        @safetestset "Utils" begin
            include("utils.jl")
        end
        @safetestset "Simulation results" begin
            include("simresults.jl")
        end
        @safetestset "Wuschel" begin
            include("wuschel.jl")
        end
    end
end
