using SBMLToolkit, Test
using SciMLTesting

run_qa(
    SBMLToolkit;
    explicit_imports = true,
    api_docs_kwargs = (; rendered = true),
    aqua_kwargs = (;
        ambiguities = (; recursive = false),
        piracies = (; treat_as_own = [SBMLToolkit.SBML.Model]),
    ),
    ei_kwargs = (;
        # `setmetadata` is needed at runtime by the `@species` macro expansion but is
        # not directly referenced in source.
        no_stale_explicit_imports = (; ignore = (:setmetadata,)),
        # `Symbolics.unwrap` re-exports `SymbolicUtils.unwrap`; on Julia 1.10 `Base.which`
        # attributes the name to its owner `SymbolicUtils`, so `via_owners` flags the
        # `Symbolics.unwrap` accesses (it does not flag on >=1.11). The source must keep
        # the `Symbolics.unwrap` spelling: at the compat floor the wrapper-aware form is
        # required for `@species`/`@parameters` results, whereas `SymbolicUtils.unwrap`
        # is a no-op there.
        all_qualified_accesses_via_owners = (; ignore = (:unwrap,)),
        # Names accessed qualified from their owner package but not (yet) declared
        # `public`/exported there. Ignored per source package; drop each as the
        # upstream library marks the name public.
        all_qualified_accesses_are_public = (;
            ignore = (
                # SBML (non-SciML dep; names not declared public upstream)
                :AlgebraicRule, :AssignmentRule, :Math, :MathApply, :MathAvogadro,
                :MathConst, :MathIdent, :MathTime, :MathVal, :Model, :RateRule,
                :Rule, :SpeciesReference, :default_constants, :default_function_mapping,
                :extensive_kinetic_math, :initial_amounts, :interpret_math, :isfreein,
                :seemsdefined,
                # Symbolics: not yet declared public upstream.
                :fixpoint_sub,
            ),
        ),
    ),
)
