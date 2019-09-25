include("../src/PiecewisePolynomials.jl")
using .PiecewisePolynomials
using Documenter

makedocs(
    sitename = "PiecewisePolynomials.jl",
    format = Documenter.HTML(prettyurls = true),
    pages = Any[
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/loooj58/PiecewisePolynomials.jl.git",
    forcepush = true
)
