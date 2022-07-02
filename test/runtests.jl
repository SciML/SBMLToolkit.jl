using SafeTestsets, Test

@safetestset "SBMLToolkit.jl" begin
    @safetestset "Systems" begin include("systems.jl") end
    @safetestset "Reactions" begin include("reactions.jl") end
    @safetestset "Rules" begin include("rules.jl") end
    @safetestset "Events" begin include("events.jl") end
    @safetestset "Utils" begin include("utils.jl") end
    # @safetestset "Misc" begin include("misc.jl") end
    @safetestset "Wuschel" begin include("wuschel.jl") end
end
