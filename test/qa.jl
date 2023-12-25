using SBMLToolkit, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(SBMLToolkit)
    Aqua.test_ambiguities(SBMLToolkit, recursive = false)
    Aqua.test_deps_compat(SBMLToolkit)
    Aqua.test_piracies(SBMLToolkit,
        treat_as_own = [SBMLToolkit.SBML.Model])
    Aqua.test_project_extras(SBMLToolkit)
    Aqua.test_stale_deps(SBMLToolkit)
    Aqua.test_unbound_args(SBMLToolkit)
    Aqua.test_undefined_exports(SBMLToolkit)
end
