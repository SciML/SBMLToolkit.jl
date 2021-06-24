using SBMLToolkit
using Test

using SBML
using Catalyst
using DifferentialEquations

@testset "SBMLToolkit.jl" begin
    include("reactionsystem.jl")
end
