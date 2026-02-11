using MMFNEXUS
using Documenter

DocMeta.setdocmeta!(MMFNEXUS, :DocTestSetup, :(using MMFNEXUS); recursive=true)

makedocs(;
    modules=[MMFNEXUS],
    authors="Bram D. Alferink <bramalferink@gmail.com> and contributors",
    sitename="MMFNEXUS.jl",
    format=Documenter.HTML(;
        canonical="https://BDAlferink.github.io/MMFNEXUS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Tutorial" => "tutorial.md",
        "citations" => "citations.md",
    ],
)

deploydocs(;
    repo="github.com/BDAlferink/MMFNEXUS.jl",
    push_preview = true, # not sure why, if problems, check out
    devbranch="main",
)
