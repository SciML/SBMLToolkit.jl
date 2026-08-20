module SBMLToolkit

using Catalyst: Catalyst, @species, Reaction, ReactionSystem, default_t, default_time_deriv,
    isspecies
import ModelingToolkit: ModelingToolkit, @parameters, Equation, ODESystem, complete
using Symbolics: Symbolics, Num
using SBML: SBML, convert_promotelocals_expandfuns, readSBML, readSBMLFromString,
    set_level_and_version
using SymbolicUtils: SymbolicUtils, expand, setmetadata, simplify, substitute

include("drafts.jl")
include("systems.jl")
include("reactions.jl")
include("rules.jl")
include("events.jl")
include("utils.jl")

"""
    convert_simplify_math(args...)

Compatibility alias for `SBML.convert_promotelocals_expandfuns`.

All arguments are forwarded unchanged. New code should import and call
`SBML.convert_promotelocals_expandfuns` directly.

# Examples

```julia
normalized_document = convert_simplify_math(document)
```
"""
const convert_simplify_math = convert_promotelocals_expandfuns

export ReactionSystem, ODESystem
export readSBML, readSBMLFromString, set_level_and_version, convert_simplify_math,
    convert_promotelocals_expandfuns, checksupport_file
export DefaultImporter, ReactionSystemImporter, ODESystemImporter

include("precompile.jl")

end
