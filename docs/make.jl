using Documenter, StructuredLight

makedocs(
    sitename="StructuredLight.jl",
    pages = [
        "index.md",
        "initial_profiles.md",
        "visualization.md",
        "propagation.md",
        "miscellany.md",
        "Examples" => [ "astig_conversion.md", "kerr_diffraction_rings.md" ]
    ]
    )