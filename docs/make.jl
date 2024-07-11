using Documenter, StructuredLight

DocMeta.setdocmeta!(StructuredLight, :DocTestSetup, :(using StructuredLight); recursive=true)

makedocs(
    modules=[StructuredLight],
    sitename="StructuredLight.jl",
    pages=[
        "index.md",
        "initial_profiles.md",
        "visualization.md",
        "propagation.md",
        "miscellany.md",
        "examples.md"]
)

deploydocs(
    repo="https://github.com/marcsgil/StructuredLight.jl.git",
)