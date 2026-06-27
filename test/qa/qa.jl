using SBMLToolkit, Test
using SciMLTesting

run_qa(
    SBMLToolkit;
    explicit_imports = true,
    aqua_kwargs = (;
        ambiguities = (; recursive = false),
        piracies = (; treat_as_own = [SBMLToolkit.SBML.Model]),
    ),
    ei_kwargs = (;
        # `setmetadata` is needed at runtime by the `@species` macro expansion but is
        # not directly referenced in source.
        no_stale_explicit_imports = (; ignore = (:setmetadata,)),
        # `unwrap`/`value`/`fixpoint_sub` are written as `Symbolics.X` in source, but
        # ExplicitImports' qualifier resolution differs across Julia versions: on 1.10
        # it attributes `Symbolics.value`/`Symbolics.fixpoint_sub` to ModelingToolkit
        # (which re-exports the `Symbolics` module name) and reports `Symbolics.unwrap`
        # as owned by SymbolicUtils, so `via_owners` flags all three on lts but none on
        # >=1.11. `unwrap` in particular must stay `Symbolics.unwrap` (not
        # `SymbolicUtils.unwrap`): at the compat floor the wrapper-aware `Symbolics`
        # form is required for `@species` results, whereas `SymbolicUtils.unwrap(::Num)`
        # is a no-op there. No single spelling satisfies the owner check on both 1.10
        # and >=1.11 without breaking floor correctness.
        all_qualified_accesses_via_owners = (;
            ignore = (:fixpoint_sub, :unwrap, :value),
        ),
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
                # SymbolicUtils / Symbolics: not yet declared public upstream.
                :fixpoint_sub, :symtype,
            ),
        ),
    ),
)
