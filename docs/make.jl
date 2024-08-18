push!(LOAD_PATH, "../src/")
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Documenter
using ControlBarrierFunctions

makedocs(;
    sitename="ControlBarrierFunctions.jl",
    pages=[
        "Home" => "index.md",
        "Systems" => "systems.md",
        "Barrier Functions" => "barriers.md",
        "Controllers" => "controllers.md",
        "Index" => "api.md",
    ],
)

deploydocs(; repo="github.com/maxhcohen/ControlBarrierFunctions.jl", devbranch="main")
