# SBMLToolkit

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/SBMLToolkit/stable/)

[![codecov](https://codecov.io/gh/SciML/SBMLToolkit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/SciML/SBMLToolkit.jl)
[![Build Status](https://github.com/SciML/SBMLToolkit.jl/workflows/CI/badge.svg)](https://github.com/SciML/SBMLToolkit.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

SBMLToolkit.jl is a lightweight tool to import models specified in the Systems Biology Markup Language (SBML) into the Julia SciML ecosystem. There are multiple ways to specify the same model in SBML. Please help us improving SBMLToolkit.jl by creating a GitHub issue if you experience errors when converting your SBML model.

SBMLToolkit uses the [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wrapper of the [libSBML](https://model.caltech.edu/software/libsbml/) library to lower dynamical SBML models into dynamical systems. If you would like to import BioNetGen models use the `writeSBML()` export function or import the `.net` file with [ReactionNetworkImporters.jl](https://github.com/SciML/ReactionNetworkImporters.jl). For constrained-based modeling, please have a look at [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl). We also recommend trying [SBMLImporter.jl](https://github.com/sebapersson/SBMLImporter.jl) (Pros: respects directionality of events, outputs concentrations instead of amounts. Cons: imports ODESystems, i.e. does not interface with Catalyst).


## Installation

To install SBMLToolkit.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("SBMLToolkit")
```

## Tutorial

SBML models can be simulated with the following steps (note that `sol` is in absolute quantities rather than concentrations):

```julia
using SBMLToolkit, ModelingToolkit, OrdinaryDiffEq

SBMLToolkit.checksupport_file("my_model.xml")
mdl = readSBML("my_model.xml", doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_promotelocals_expandfuns(doc)
end)

rs = ReactionSystem(mdl)  # If you want to create a reaction system
odesys = convert(ODESystem, rs)  # Alternatively: ODESystem(mdl)

tspan = (0.0, 1.0)
prob = ODEProblem(odesys, [], tspan, [])
sol = solve(prob, Tsit5())
```

Alternatively, SBMLToolkit also provides more concise methods to import `SBML.Models`, `Catalyst.ReactionSystems`, and `ModelingToolkit.ODESystems`.

```julia
mdl = readSBML(sbmlfile, DefaultImporter())
rs = readSBML(sbmlfile, ReactionSystemImporter())
odesys = readSBML(sbmlfile, ODESystemImporter())
```

## License

The package is released under the [MIT license](https://github.com/SciML/SBMLToolkit.jl/blob/main/LICENSE).

## Development team

This package was developed by [Paul F. Lang](https://www.linkedin.com/in/paul-lang-7b54a81a3/) at the University of Oxford, UK and [Anand Jain](https://github.com/anandijain) at the University of Chicago, USA.

## Questions and comments

Please use GitHub issues, the #sciml-sysbio channel in the [Julia Slack workspace](https://julialang.org/slack/) or email [Paul F. Lang](mailto:plang@biosim.ai) with any questions or comments.
