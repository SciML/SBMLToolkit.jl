var documenterSearchIndex = {"docs":
[{"location":"api/#API-documentation","page":"API documentation","title":"API documentation","text":"","category":"section"},{"location":"api/","page":"API documentation","title":"API documentation","text":"Modules = [SBMLToolkit]\nPages = [\"systems.jl\"]","category":"page"},{"location":"api/#Catalyst.ReactionSystem-Tuple{SBML.Model}","page":"API documentation","title":"Catalyst.ReactionSystem","text":"ReactionSystem(model::SBML.Model; kwargs...)\n\nCreate a ReactionSystem from an SBML.Model.\n\nSee also ODESystem.\n\n\n\n\n\n","category":"method"},{"location":"api/#ModelingToolkit.ODESystem-Tuple{SBML.Model}","page":"API documentation","title":"ModelingToolkit.ODESystem","text":"ODESystem(model::SBML.Model; include_zero_odes = true, kwargs...)\n\nCreate an ODESystem from an SBML.Model.\n\nSee also ReactionSystem.\n\n\n\n\n\n","category":"method"},{"location":"api/#SBMLToolkit.checksupport_file-Tuple{String}","page":"API documentation","title":"SBMLToolkit.checksupport_file","text":"checksupport_file(filename::String)\n\nCheck if SBML file is supported by SBMLToolkit.jl.\n\n\n\n\n\n","category":"method"},{"location":"api/#SBMLToolkit.checksupport_string-Tuple{String}","page":"API documentation","title":"SBMLToolkit.checksupport_string","text":"checksupport_string(filename::String)\n\nCheck if SBML passed as string is supported by SBMLToolkit.jl.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SBMLToolkit","category":"page"},{"location":"#SBMLToolkit","page":"Home","title":"SBMLToolkit","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SBMLToolkit.jl is a lightweight tool to import models specified in the Systems Biology Markup Language (SBML) into the Julia SciML ecosystem. There are multiple ways to specify the same model in SBML. Please help us improve SBMLToolkit.jl by creating a GitHub issue if you experience errors when converting your SBML model.","category":"page"},{"location":"","page":"Home","title":"Home","text":"SBMLToolkit uses the SBML.jl wrapper of the libSBML library to lower dynamical SBML models into dynamical systems. If you would like to import BioNetGen models, use the writeSBML() export function or import the .net file with ReactionNetworkImporters.jl. For constrained-based modelling, please have a look at COBREXA.jl.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install SBMLToolkit.jl, use the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"SBMLToolkit\")","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please refer to the SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages for guidance on PRs, issues, and other matters relating to contributing to SciML.\nSee the SciML Style Guide for common coding practices and other style decisions.\nThere are a few community forums:\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Slack\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Zulip\nOn the Julia Discourse forums\nSee also SciML Community page","category":"page"},{"location":"#Tutorial","page":"Home","title":"Tutorial","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SBML models can be simulated with the following steps (note that sol is in absolute quantities rather than concentrations):","category":"page"},{"location":"","page":"Home","title":"Home","text":"using SBMLToolkit, ModelingToolkit, OrdinaryDiffEq\n\nSBMLToolkit.checksupport_file(\"my_model.xml\")\nmdl = readSBML(\"my_model.xml\", doc -> begin\n                   set_level_and_version(3, 2)(doc)\n                   convert_promotelocals_expandfuns(doc)\n               end)\nrs = ReactionSystem(mdl)  # If you want to create a reaction system\nodesys = convert(ODESystem, rs)  # Alternatively: ODESystem(mdl)\n\ntspan = (0.0, 1.0)\nprob = ODEProblem(odesys, [], tspan, [])\nsol = solve(prob, Tsit5())","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is released under the MIT license.","category":"page"},{"location":"#Development-team","page":"Home","title":"Development team","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package was developed by Paul F. Lang at the University of Oxford, UK and Anand Jain at the University of Chicago, USA.","category":"page"},{"location":"#Questions-and-comments","page":"Home","title":"Questions and comments","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please use GitHub issues, the #sciml-sysbio channel in the Julia Slack workspace or email Paul F. Lang or Anand Jain with any questions or comments.","category":"page"},{"location":"#Reproducibility","page":"Home","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status(; mode = PKGMODE_MANIFEST) # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can also download the \n<a href=\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nversion = TOML.parse(read(\"../../Project.toml\", String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\", String))[\"name\"]\nlink = \"https://github.com/SciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n       \"/assets/Manifest.toml\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"\">manifest</a> file and the\n<a href=\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nversion = TOML.parse(read(\"../../Project.toml\", String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\", String))[\"name\"]\nlink = \"https://github.com/SciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n       \"/assets/Project.toml\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"\">project</a> file.","category":"page"}]
}
