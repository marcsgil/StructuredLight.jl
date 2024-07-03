using Documenter, StructuredLight

makedocs(
    sitename="StructuredLight.jl",
    pages=[
        "index.md",
        "initial_profiles.md",
        "propagation.md",
        "miscellany.md",
        "examples.md"]
)

deploydocs(
    repo="https://github.com/marcsgil/StructuredLight.jl.git",
)