using Pkg
Pkg.activate(".")
using Documenter
using CBFToolbox

push!(LOAD_PATH, "../src/")
makedocs(;
    sitename="CBFToolbox.jl",
    pages=[
        "Home" => "index.md",
        "Systems" => "systems.md",
        "Barrier Functions" => "barriers.md",
        "Controllers" => "controllers.md",
        "Index" => "api.md",
    ],
)

# deploydocs(; repo="github.com/maxhcohen/CBFToolbox.jl.git", devbranch="burn-it-down")
