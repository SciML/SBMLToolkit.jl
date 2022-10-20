using SBMLToolkit
using Documenter

DocMeta.setdocmeta!(SBMLToolkit, :DocTestSetup, :(using SBMLToolkit); recursive = true)

makedocs(;
         modules = [SBMLToolkit],
         authors = "paulflang, anandijain",
         sitename = "SBMLToolkit.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://docs.sciml.ai/SBMLToolkit/stable/",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
             "API documentation" => "api.md",
         ])

deploydocs(;
           repo = "github.com/SciML/SBMLToolkit.jl")
