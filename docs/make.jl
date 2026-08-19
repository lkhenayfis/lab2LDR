using lab2ldr
using Documenter

DocMeta.setdocmeta!(lab2ldr, :DocTestSetup, :(using lab2ldr); recursive=true)

makedocs(;
    modules=[lab2ldr],
    authors="lucas-khenayfis",
    sitename="lab2ldr.jl",
    format=Documenter.HTML(;
        canonical="https://lkhenayfis.github.io/lab2ldr.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lkhenayfis/lab2ldr.jl",
    devbranch="main",
)
