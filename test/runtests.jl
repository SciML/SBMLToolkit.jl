using SBMLToolkit
using Test

using SBML
using Catalyst
using OrdinaryDiffEq
using SBML

@testset "SBMLToolkit.jl" begin
    include("reactionsystem.jl")
end
