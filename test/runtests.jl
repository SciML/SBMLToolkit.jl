using SBMLToolkit
using Test

using Catalyst
using OrdinaryDiffEq
using ModelingToolkit
using SBML

@testset "SBMLToolkit.jl" begin
    @testset "Model to MTK conversions" begin include("reactionsystem.jl") end
end
