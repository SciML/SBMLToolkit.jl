module SBMLToolkit

using Catalyst
using SBML
using SymbolicUtils
using Setfield

include("drafts.jl")
include("systems.jl")
include("reactions.jl")
include("rules.jl")
include("events.jl")
include("utils.jl")
include("nameswap.jl")

@deprecate convert_simplify_math convert_promotelocals_expandfuns

export ReactionSystem, ODESystem
export readSBML, readSBMLFromString, set_level_and_version, convert_simplify_math,
       convert_promotelocals_expandfuns
export DefaultImporter, ReactionSystemImporter, ODESystemImporter

end
