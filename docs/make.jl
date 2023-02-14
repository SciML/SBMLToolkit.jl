using SBMLToolkit
using Documenter

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

DocMeta.setdocmeta!(SBMLToolkit, :DocTestSetup, :(using SBMLToolkit); recursive = true)

makedocs(;
         modules = [SBMLToolkit],
         authors = "paulflang, anandijain",
         sitename = "SBMLToolkit.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://docs.sciml.ai/SBMLToolkit/stable/",
                                  assets = ["assets/favicon.ico"]),
         clean = true, doctest = false, linkcheck = true,
         linkcheck_ignore = ["https://www.linkedin.com/in/paul-lang-7b54a81a3/"],
         strict = [
             :doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # Other available options are
             # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
         ],
         pages = [
             "Home" => "index.md",
             "API documentation" => "api.md",
         ])

deploydocs(;
           repo = "github.com/SciML/SBMLToolkit.jl")
