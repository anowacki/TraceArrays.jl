using TraceArrays
using Documenter

DocMeta.setdocmeta!(TraceArrays, :DocTestSetup, :(using TraceArrays); recursive=true)

makedocs(;
    modules=[TraceArrays],
    authors="Andy Nowacki <a.nowacki@leeds.ac.uk> and contributors",
    repo="https://github.com/anowacki/TraceArrays.jl/blob/{commit}{path}#{line}",
    sitename="TraceArrays.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://anowacki.github.io/TraceArrays.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/anowacki/TraceArrays.jl",
    devbranch="main",
)
