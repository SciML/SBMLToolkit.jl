using SBMLToolkit
using Documenter

DocMeta.setdocmeta!(SBMLToolkit, :DocTestSetup, :(using SBMLToolkit); recursive = true)

makedocs(;
         modules = [SBMLToolkit],
         authors = "paulflang, anandijain",
         repo = "https://github.com/paulflang/SBMLToolkit.jl/blob/{commit}{path}#{line}",
         sitename = "SBMLToolkit.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://paulflang.github.io/SBMLToolkit.jl",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
         ])

deploydocs(;
           repo = "github.com/paulflang/SBMLToolkit.jl")
