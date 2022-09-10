# SBMLToolkit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://paulflang.github.io/SBMLToolkit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://paulflang.github.io/SBMLToolkit.jl/dev)
[![Build Status](https://github.com/paulflang/SBMLToolkit.jl/workflows/CI/badge.svg)](https://github.com/paulflang/SBMLToolkit.jl/actions)
[![Coverage](https://codecov.io/gh/paulflang/SBMLToolkit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/paulflang/SBMLToolkit.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

SBMLToolkit.jl is a lightweight tool to import models specified in the Systems Biology Markup Language (SBML) into the Julia SciML ecosystem. More specifically, SBMLToolkit.jl extracts reactions, initial conditions and parameter values from SBML files. Events, Rules and several other SBML components are not yet supported. For details on support of SBML test suite cases, please refer to the [SBML Test-Suite Results Tracker](https://github.com/SciML/SBMLToolkit.jl/issues/45). There are multiple ways to specify the same model in SBML. Please help us improving SBMLToolkit.jl by creating a GitHub issue if you experience errors when converting your SBML model.

SBMLToolkit uses the [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wrapper of the [libSBML](https://model.caltech.edu/software/libsbml/) library to lower dynamical SBML models into dynamical systems. If you are interested in constrained-based modelling please have a look at [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl).

## Installation
SBMLToolkit.jl is available on the Julia package managing system. To install SBMLToolkit run the following in the REPL:
  ```
  ]add SBMLToolkit
  ```

## Tutorial
SBML models can be simulated with the following steps (note that `sol` is in absolute quantities rather than concentrations):
  ```julia
  using SBMLToolkit, ModelingToolkit, OrdinaryDiffEq

  SBMLToolkit.checksupport_file("my_model.xml")
  mdl = readSBML("my_model.xml", doc -> begin
      set_level_and_version(3, 2)(doc)
      convert_simplify_math(doc)
  end)
  rs = ReactionSystem(mdl)
  odesys = convert(ODESystem, rs)

  tspan = [0., 1.]
  prob = ODEProblem(odesys, [], tspan, [])
  sol = solve(prob, Tsit5())
  ```


## License
The package is released under the [MIT license](https://github.com/paulflang/SBMLToolkit.jl/blob/main/LICENSE).


## Development team
This package was developed by [Paul F. Lang](https://www.linkedin.com/in/paul-lang-7b54a81a3/) at the University of Oxford, UK and [Anand Jain](https://github.com/anandijain) at the University of Chicago, USA.


## Questions and comments
Please contact [Paul F. Lang](mailto:paul.lang@wolfson.ox.ac.uk) or [Anand Jain](mailto:anandj@uchicago.edu) with any questions or comments.
