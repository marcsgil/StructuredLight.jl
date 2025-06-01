using Documenter, StructuredLight

DocMeta.setdocmeta!(StructuredLight, :DocTestSetup, :(using StructuredLight, CairoMakie); recursive=true)

makedocs(
    sitename="StructuredLight.jl",
    pages=[
        "index.md",
        "quick_start.md",
        "beam_profiles.md",
        "visualization.md",
        "propagation.md",
        "holograms.md",
        "aberration_correction.md",
        "gpu_support.md",
        "miscellany.md",
        "examples.md",
    ],
    draft=false
)

deploydocs(
    repo="https://github.com/marcsgil/StructuredLight.jl.git",
)