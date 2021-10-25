# SBMLToolkit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://paulflang.github.io/SBMLToolkit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://paulflang.github.io/SBMLToolkit.jl/dev)
[![Build Status](https://github.com/paulflang/SBMLToolkit.jl/workflows/CI/badge.svg)](https://github.com/paulflang/SBMLToolkit.jl/actions)
[![Coverage](https://codecov.io/gh/paulflang/SBMLToolkit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/paulflang/SBMLToolkit.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

SBMLToolkit.jl is a lightweight tool to import models specified in the Systems Biology Markup Language (SBML) into the Julia SciML ecosystem. More specifically, SBMLToolkit.jl extracts reactions, initial conditions and parameter values from SBML files. Events, Rules and several other SBML components are not yet supported. There are multiple ways to specify the same model in SBML. Please help us improving SBMLToolkit.jl by creating a GitHub issue if you experience errors when converting your SBML model.

SBMLToolkit uses the [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wrapper of the [libSBML](https://model.caltech.edu/software/libsbml/) library to lower dynamical SBML models into dynamical systems. If you are interested in constrained-based modelling please have a look at [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl).

## Installation
SBMLToolkit.jl is available on the Julia package managing system. To install SBMLToolkit run the following in the REPL:
  ```
  ]add SBMLToolkit
  ```

## Tutorial
SBML models can be simulated with the following steps (note that `sol` is in absolute quantities rather than concentrations):
  ```julia
  using SBMLToolkit, ModelingToolkit

  SBMLToolkit.checksupport("my_model.xml")
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

## SBML test suite report
SBMLtoolkit.jl is regularly tested against the [SBML test suite](https://github.com/sbmlteam/sbml-test-suite/tree/master/cases/semantic). The last test was run on October 24, 2021 on cases 00001 to 01664. Detailed results can be found in [test_suite_report.csv](./test_suite_report.csv). Out of the 1664 cases 1219 expectedly failed, as
* they contain SBML constraints, rules or events, or
* do not contain any reactions, or
* cannot be converted to SBML Version 3 Level 2, or
* contain kineticLaws in reversible reactions with ambiguous forward and reverse part, or
* are flux balance models.

Of the remaining 445 cases
* 422 can be converted to a `ModelingToolkit.ReactionSystem`, of which
* 416 simulate without returning an error, of which
* 416 yield the correct simulation results.

## License
The package is released under the [MIT license](https://github.com/paulflang/SBMLToolkit.jl/blob/main/LICENSE).


## Development team
This package was developed by [Paul F. Lang](https://www.linkedin.com/in/paul-lang-7b54a81a3/) at the University of Oxford, UK and [Anand Jain](https://github.com/anandijain) at the University of Chicago, USA.


## Questions and comments
Please contact [Paul F. Lang](mailto:paul.lang@wolfson.ox.ac.uk) or [Anand Jain](mailto:anandj@uchicago.edu) with any questions or comments.
