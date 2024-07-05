# Introduction

This package provides tools to simulate the propagation of paraxial light beams. This includes the calculation of [Laguerre-Gauss](@ref) and [Hermite-Gauss](@ref) beam profiles, the action of [Lenses](@ref), the propagation in free space ([`free_propagation`](@ref)) as well as in Kerr media ([`kerr_propagation`](@ref)). Both the propagation methods can be run on Nvidia GPUs (check [CUDA support](@ref) for details). We also provide methods that help the [Visualization](@ref) of such beams.

## Minimal Example

The following code is a minimal working example for this package:

```julia
using StructuredLight #Loads the package

rs = LinRange(-5,5,256) #Define a linear grid of points

ψ₀ = lg(rs,rs) #Calculates the fundamental Laguerre-Gaussian mode

visualize(abs2.(ψ₀)) #visualizes the mode

ψ = free_propagation(ψ₀,rs,rs,1) #Propagates the mode through a distance of z=1

visualize(abs2.(ψ₀)) #visualizes the evolved mode
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