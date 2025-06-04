# StructuredLight.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://marcsgil.github.io/StructuredLight.jl)
[![CI](https://github.com/marcsgil/StructuredLight.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/marcsgil/StructuredLight.jl/actions/workflows/CI.yml)

A comprehensive Julia package for simulating, manipulating, and visualizing spatially structured light beams. StructuredLight.jl provides efficient tools for paraxial beam propagation, including specialized functions for generating structured modes, simulating optical elements, and creating computer-generated holograms for spatial light modulators.

## Features

- **Beam Profiles**: Generate Laguerre-Gauss and Hermite-Gauss modes, diagonal Hermite-Gauss modes, and various aperture shapes
- **Free Space Propagation**: Simulate beam propagation using fast Fourier transforms
- **Nonlinear Propagation**: Kerr effect simulation for nonlinear media
- **Computer-Generated Holograms**: Create phase patterns for spatial light modulators with multiple algorithms
- **Phase Modulation**: Apply lenses, tilted lenses, and aberration correction using Zernike polynomials
- **Linear Combinations**: Efficiently compute superpositions of structured modes
- **GPU Acceleration**: Full CUDA support for high-performance computations
- **Visualization**: Comprehensive plotting tools via Makie.jl extension

## Installation

```julia
using Pkg
Pkg.add("StructuredLight")
```

## Quick Start

```julia
using StructuredLight

# Create a spatial grid
rs = LinRange(-5, 5, 256)

# Generate fundamental Gaussian mode
E0 = lg(rs, rs)

# Visualize the beam intensity
visualize(abs2.(E0))

# Propagate through free space
E_propagated = free_propagation(E0, rs, rs, 1.0)
```

## Advanced Examples

### GPU Acceleration
```julia
using CUDA

# Transfer to GPU for faster computation
E0_gpu = lg(rs, rs) |> cu
E_propagated = free_propagation(E0_gpu, rs, rs, 1.0)
```

### Computer-Generated Holograms
```julia
# Create a hologram for an SLM
target_field = lg(rs, rs, p=1, l=2)  # LG mode
hologram = generate_hologram(target_field, 255, 20, 20)
```

## Compatibility

- **Julia**: 1.10 - 1.12
- **GPU Support**: CUDA.jl for NVIDIA GPUs
- **Visualization**: Makie.jl ecosystem
- **Documentation**: [Full documentation available](https://marcsgil.github.io/StructuredLight.jl)

## Related Packages

- [SpatialLightModulator.jl](https://github.com/marcsgil/SpatialLightModulator.jl): Interface with physical SLM devices
