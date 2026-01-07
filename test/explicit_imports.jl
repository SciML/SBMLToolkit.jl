using ExplicitImports
using SBMLToolkit
using Test

@testset "ExplicitImports" begin
    @testset "No implicit imports" begin
        @test check_no_implicit_imports(SBMLToolkit) === nothing
    end

    @testset "No stale explicit imports" begin
        @test check_no_stale_explicit_imports(SBMLToolkit) === nothing
    end
end
