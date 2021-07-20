using SBMLToolkit
using Test

using SBML
using Catalyst
using OrdinaryDiffEq
using SBML

@testset "SBMLToolkit.jl" begin
    @testset "Model to MTK conversions" begin include("reactionsystem.jl") end
end
