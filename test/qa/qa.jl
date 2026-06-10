using SBMLToolkit, Aqua, ExplicitImports, Test

@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(SBMLToolkit)
    Aqua.test_ambiguities(SBMLToolkit, recursive = false)
    Aqua.test_deps_compat(SBMLToolkit)
    Aqua.test_piracies(
        SBMLToolkit,
        treat_as_own = [SBMLToolkit.SBML.Model]
    )
    Aqua.test_project_extras(SBMLToolkit)
    Aqua.test_stale_deps(SBMLToolkit)
    Aqua.test_unbound_args(SBMLToolkit)
    Aqua.test_undefined_exports(SBMLToolkit)
end

@testset "ExplicitImports" begin
    @testset "No implicit imports" begin
        @test check_no_implicit_imports(SBMLToolkit) === nothing
    end

    @testset "No stale explicit imports" begin
        # setmetadata is needed at runtime by @species macro expansion but not directly referenced in source
        @test check_no_stale_explicit_imports(SBMLToolkit; ignore = (:setmetadata,)) === nothing
    end
end
