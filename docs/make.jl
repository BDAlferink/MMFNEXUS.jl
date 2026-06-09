using MMFNexus
using Documenter

DocMeta.setdocmeta!(MMFNexus, :DocTestSetup, :(using MMFNexus); recursive=true)

makedocs(;
    modules=[MMFNexus],
    authors="Bram D. Alferink <bramalferink@gmail.com> and contributors",
    sitename="MMFNexus.jl",
    format=Documenter.HTML(;
        canonical="https://BDAlferink.github.io/MMFNexus.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Tutorial" => "tutorial.md",
        "citations" => "citations.md",
    ],
    checkdocs = :none # disable check for now

)

deploydocs(;
    repo="github.com/BDAlferink/MMFNexus.jl",
    push_preview = true, # not sure why, if problems, check out
    devbranch="main",
)
