using SBMLToolkit
using Test

using SBML
using Catalyst
using Downloads
using ModelingToolkit
using OrdinaryDiffEq

@testset "SBMLToolkit.jl" begin
    @testset "Model to MTK conversions" begin include("reactionsystem.jl") end
    @testset "Rules" begin include("rules.jl") end
    @testset "Events" begin include("events.jl") end
    @testset "Wuschel" begin include("wuschel.jl") end
end
