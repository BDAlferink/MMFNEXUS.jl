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
    ],
)

deploydocs(;
    repo="github.com/BDAlferink/MMFNEXUS.jl",
    devbranch="main",
)
