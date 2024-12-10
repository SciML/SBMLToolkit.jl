# SBMLToolkit

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/SBMLToolkit/stable/)

[![codecov](https://codecov.io/gh/SciML/SBMLToolkit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/SciML/SBMLToolkit.jl)
[![Build Status](https://github.com/SciML/SBMLToolkit.jl/workflows/CI/badge.svg)](https://github.com/SciML/SBMLToolkit.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

SBMLToolkit.jl is a lightweight tool to import models specified in the Systems Biology Markup Language (SBML) into the Julia SciML ecosystem. There are multiple ways to specify the same model in SBML. Please help us improving SBMLToolkit.jl by creating a GitHub issue if you experience errors when converting your SBML model.

SBMLToolkit uses the [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wrapper of the [libSBML](https://model.caltech.edu/software/libsbml/) library to lower dynamical SBML models into completed dynamical systems. If you would like to import BioNetGen models use the `writeSBML()` export function or import the `.net` file with [ReactionNetworkImporters.jl](https://github.com/SciML/ReactionNetworkImporters.jl). For constrained-based modeling, please have a look at [COBREXA.jl](https://github.com/COBREXA/COBREXA.jl). We also recommend trying [SBMLImporter.jl](https://github.com/sebapersson/SBMLImporter.jl). While SBMLToolkit.jl has a slightly cleaner interface, SBMLImporter.jl respects directionality of events, can output concentrations in addition to amounts, and provides better simulation performance for models including time-triggered events and SBML `piecewise` expressions. If you are an experienced SBML user and interested if SBMLToolkit supports certain (SBML test suite)[https://github.com/sbmlteam/sbml-test-suite] cases, you can download the logs from the latest CI run of the [SBMLToolkitTestSuite](https://github.com/SciML/SBMLToolkitTestSuite.jl/actions/workflows/CI.yml).

## Installation

To install SBMLToolkit.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("SBMLToolkit")
```

## Tutorial

SBML models can be simulated with the following steps (note that `sol` is in absolute quantities rather than concentrations):

```julia
using SBMLToolkit, OrdinaryDiffEq

odesys = readSBML("my_model.xml", ODESystemImporter())

tspan = (0.0, 1.0)
prob = ODEProblem(odesys, [], tspan, [])
sol = solve(prob, Tsit5())
```

While this imports an `ODESystem` directly, you can also import a Catalyst.jl `ReactionSystem`:

```julia
using SBMLToolkit

rs = readSBML("my_model.xml", ReactionSystemImporter())
```

One common case where this is useful is if you want to run stochastic instead of ODE simulations.

In the very unlikely case that you need fine-grained control over the SBML import, you can create an SBML.jl `Model` (we strongly recommend manually running `checksupport_file("my_model.xml")` before)

```julia
using SBML

mdl = readSBML("my_model.xml", doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_promotelocals_expandfuns(doc)
end)
```

The conversion to SBML level 3 version 2 is necessary, because older versions are not well supported in SBMLToolkit. `convert_promotelocals_expandfuns` basically flattens the SBML before the import. Once you have obtained the `Model`, you can convert it to a `ReactionSystem` and `ODESystem`.

```julia
using SBMLToolkit

rs = ReactionSystem(mdl)
odesys = convert(ODESystem, rs)
```

## License

The package is released under the [MIT license](https://github.com/SciML/SBMLToolkit.jl/blob/main/LICENSE).

## Questions and comments

Please use GitHub issues and the #sciml-sysbio channel in the [Julia Slack workspace](https://julialang.org/slack/) with any questions or comments.

# Citation

If you use SBMLToolkit.jl in your research, please cite [this paper](https://www.degruyter.com/document/doi/10.1515/jib-2024-0003/html):

```
@article{lang_sbmltoolkitjl_2024,
	title = {{SBMLToolkit}.jl: a {Julia} package for importing {SBML} into the {SciML} ecosystem},
	doi = {10.1515/jib-2024-0003},
	journal = {Journal of Integrative Bioinformatics},
	author = {Lang, Paul F. and Jain, Anand and Rackauckas, Christopher},
	year = {2024},
}
```
