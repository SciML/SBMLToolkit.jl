```@meta
CurrentModule = SBMLToolkit
```

# SBMLToolkit

SBMLToolkit.jl is a lightweight tool to import models specified in the Systems Biology Markup Language (SBML) into the Julia SciML ecosystem. There are multiple ways to specify the same model in SBML. Please help us improve SBMLToolkit.jl by creating a GitHub issue if you experience errors when converting your SBML model.

SBMLToolkit uses the [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl) wrapper of the [libSBML](https://sbml.org/software/libsbml/) library to lower dynamical SBML models into dynamical systems. If you would like to import BioNetGen models, use the `writeSBML()` export function or import the `.net` file with [ReactionNetworkImporters.jl](https://github.com/SciML/ReactionNetworkImporters.jl). For constrained-based modeling, please have a look at [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl).

## Installation

To install SBMLToolkit.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("SBMLToolkit")
```

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Tutorial

SBML models can be simulated with the following steps (note that `sol` is in absolute quantities rather than concentrations):

```julia
using SBMLToolkit, ModelingToolkit, OrdinaryDiffEq

SBMLToolkit.checksupport_file("my_model.xml")
mdl = readSBML("my_model.xml", doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_promotelocals_expandfuns(doc)
end)
rs = complete(ReactionSystem(mdl))  # If you want to create a reaction system
odesys = convert(ODESystem, rs)  # Alternatively: ODESystem(mdl)

tspan = (0.0, 1.0)
prob = ODEProblem(odesys, [], tspan, [])
sol = solve(prob, Tsit5())
```

SBMLToolkit also provides the following convenience functions to import `SBML.Models`, `Catalyst.ReactionSystems` and `ModelingToolkit.ODESystems`:

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

Please use GitHub issues, the #sciml-sysbio channel in the [Julia Slack workspace](https://julialang.org/slack/) or email [Paul F. Lang](mailto:paul.lang@juliacomputing.com) or [Anand Jain](mailto:anandj@uchicago.edu) with any questions or comments.

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
