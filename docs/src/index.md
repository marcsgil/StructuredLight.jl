# Introduction

This package provides tools to simulate the propagation of paraxial light beams. This includes the calculation of [Laguerre-Gauss](@ref) and [Hermite-Gauss](@ref) beam profiles, the propagation in free space ([`free_propagation`](@ref)) as well as in Kerr media ([`kerr_propagation`](@ref)). Both the propagation methods can be run on Nvidia GPUs (check [CUDA support](@ref) for details). We also provide methods that help the [Visualization](@ref) of such beams.

## Minimal Example

The following code is a minimal working example for this package:

```@example
using StructuredLight

rs = LinRange(-5,5,256) #Define a linear grid of points

ψ₀ = lg(rs,rs) #Calculates the fundamental Laguerre-Gaussian mode

ψ = free_propagation(ψ₀,rs,rs,1) #Propagates the mode through a distance of z=1

using CairoMakie # for the visualization

fig = Figure(size = (1600,800), figure_padding=0)
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
hidedecorations!(ax1)
hidedecorations!(ax2)

heatmap!(ax1, abs2.(ψ₀), colormap = :jet) #visualizes the initial mode
heatmap!(ax2, abs2.(ψ), colormap = :jet) #visualizes the propagated mode

fig
```

This illustrates the basic idea of the package: first, you construct a matrix representing the mode you want to propagate, and then one calls a propagation function, in this case `free_propagation`. Then, we visualize the propagated beam.

## Contents

```@contents
Pages = [
        "index.md",
        "initial_profiles.md",
        "visualization.md",
        "propagation.md",
        "miscellany.md",
        "examples.md"
    ]
```
