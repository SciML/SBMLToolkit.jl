using SBMLToolkit
using Test

@testset "Precompile workload" begin
    @test isnothing(SBMLToolkit._precompile_workload())
end
