using Documenter, NBodySimulator

include("pages.jl")

makedocs(; sitename = "NbodySimulator.jl",
         authors = "Chris Rackauckas",
         modules = [NBodySimulator],
         clean = true, doctest = false,
         format = Documenter.HTML(; analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://nbodysimulator.sciml.ai/stable/"),
         pages = pages)

deploydocs(; repo = "github.com/SciML/NBodySimulator.jl.git",
           push_preview = true)
