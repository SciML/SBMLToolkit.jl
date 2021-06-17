using SBMLToolkit
using Test

using Catalyst
using OrdinaryDiffEq
using SBML

@testset "SBMLToolkit.jl" begin
    include("reactionsystem.jl")
end
