using Documenter, NBodySimulator

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "NBodySimulator.jl",
    authors = "Chris Rackauckas",
    modules = [NBodySimulator],
    clean = true, doctest = false, linkcheck = true,
    warnonly = [:docs_block, :missing_docs],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/NBodySimulator/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/NBodySimulator.jl.git";
    push_preview = true)
