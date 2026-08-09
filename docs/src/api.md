```@meta
CurrentModule = SBMLToolkit
```

# API

SBMLToolkit provides importer selectors and extends the system constructors owned by SBML,
Catalyst, and ModelingToolkit.

```@autodocs
Modules = [SBMLToolkit]
Pages = ["SBMLToolkit.jl", "systems.jl"]
Private = false
```

## Compatibility reexports

SBMLToolkit keeps these dependency-owned bindings available for compatibility. New code
should import each binding from its owning package.

- `ReactionSystem`: the `Catalyst.ReactionSystem` type and constructor.
- `ODESystem`: the `ModelingToolkit.ODESystem` type and constructor.
- `readSBML`: the `SBML.readSBML` loader extended by the importer selectors above.
- `readSBMLFromString`: the `SBML.readSBMLFromString` loader for in-memory XML.
- `set_level_and_version`: the corresponding SBML document transformation from SBML.jl.
- `convert_promotelocals_expandfuns`: the SBML.jl normalization transformation.
- `convert_simplify_math`: SBMLToolkit's compatibility alias for
  `convert_promotelocals_expandfuns`.
