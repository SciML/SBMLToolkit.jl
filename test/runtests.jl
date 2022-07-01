using SafeTestsets, Test

@safetestset "SBMLToolkit.jl" begin
    @safetestset "Model to MTK conversions" begin include("reactionsystem.jl") end
    @safetestset "Rules" begin include("rules.jl") end
    @safetestset "Events" begin include("events.jl") end
    @safetestset "Misc" begin include("misc.jl") end
    @safetestset "Wuschel" begin include("wuschel.jl") end
end
