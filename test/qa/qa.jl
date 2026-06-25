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
        # Names accessed qualified from their owner package but not (yet) declared
        # `public`/exported there. Ignored per source package; drop each as the
        # upstream library marks the name public.
        all_qualified_accesses_are_public = (;
            ignore = (
                # SBML
                :AlgebraicRule, :AssignmentRule, :Math, :MathApply, :MathAvogadro,
                :MathConst, :MathIdent, :MathTime, :MathVal, :Model, :RateRule,
                :Rule, :SpeciesReference, :default_constants, :default_function_mapping,
                :extensive_kinetic_math, :initial_amounts, :interpret_math, :isfreein,
                :seemsdefined,
                # SymbolicUtils / Symbolics. `unwrap` is accessed via both: the
                # wrapper-aware `Symbolics.unwrap` (utils/events) and
                # `SymbolicUtils.unwrap` (systems); they differ at the compat floor.
                :BasicSymbolic, :isadd, :ismul, :symtype,
                :fixpoint_sub, :get_variables, :unwrap, :value,
            ),
        ),
    ),
)
