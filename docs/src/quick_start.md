# Quick Start

## What is StructuredLight.jl?

StructuredLight.jl is a Julia package for simulating structured optical beams - light fields with complex spatial patterns like orbital angular momentum, customizable intensity profiles, and controlled phase distributions. The package enables:

- Generation of standard beam modes (Laguerre-Gaussian, Hermite-Gaussian)
- Beam propagation through free space and nonlinear media
- Application of optical elements (lenses, apertures, phase masks)
- Creation of computer-generated holograms for producing these beams in the lab
- Visualization and analysis of beam evolution

## Installation

To install the `StructuredLight` package, you can use the Julia package manager. Open the Julia REPL type `]` to enter the package manager mode, and then run:

```julia-repl
add StructuredLight
```

## Basic Examples

### Example 1: Creating and Visualizing a Gaussian Beam

Let's start with the simplest case - a fundamental Gaussian beam:

```@example
using StructuredLight

# Define spatial coordinates (in arbitrary units)
# This creates a 2D grid from -5 to +5 with 256 points in each direction
spatial_range = LinRange(-3, 3, 256)

# Create a fundamental Gaussian beam (Laguerre-Gaussian with l=0, p=0)
# The beam is evaluated at the focal plane (z=0)
gaussian_beam = lg(spatial_range, spatial_range)

# Visualize the intensity profile
using CairoMakie
visualize(abs2.(gaussian_beam))
```

### Example 2: Beam Propagation

Now let's see how this beam evolves as it propagates through space:

```@example
using StructuredLight, CairoMakie

# Define spatial grid
spatial_range = LinRange(-5, 5, 256)

# Create initial Gaussian beam at focus (z=0)
# We use the same spatial range for both x and y coordinates
initial_beam = lg(spatial_range, spatial_range)

# Propagate the beam to distance z=1
# Here we assume that the wavenumber k=1 for simplicity. This can be adjusted with a keyword argument `k`.
propagated_beam = free_propagation(initial_beam, spatial_range, spatial_range, 1)

# Compare initial and propagated beams side by side
fig = Figure(size = (800, 400))
ax1 = Axis(fig[1,1], title="At focus (z=0)")
ax2 = Axis(fig[1,2], title="Propagated (z=1)")

heatmap!(ax1, abs2.(initial_beam), colormap = :jet)
heatmap!(ax2, abs2.(propagated_beam), colormap = :jet)

fig
```

### Example 3: Vortex Beam with Orbital Angular Momentum

One of the most studied structured light beams are the ones carrying orbital angular momentum:

```@example
using StructuredLight, CairoMakie

spatial_range = LinRange(-5, 5, 256)

# Create a vortex beam with topological charge l=3
# This beam carries 3ℏ orbital angular momentum per photon
vortex_beam = lg(spatial_range, spatial_range, l=3)

# Show both intensity and phase
fig = Figure(size = (800, 400))
ax1 = Axis(fig[1,1], title="Intensity |ψ|²")
ax2 = Axis(fig[1,2], title="Phase arg(ψ)")

heatmap!(ax1, abs2.(vortex_beam), colormap = :jet)
heatmap!(ax2, angle.(vortex_beam), colormap = :hsv)

fig
```

### Example 4: Beam Shaping with Apertures

You can shape beams using various apertures and masks:

```@example
using StructuredLight, CairoMakie

spatial_range = LinRange(-4, 4, 256)

# Create a vortex beam
beam = lg(spatial_range, spatial_range, l=2)

# Apply different apertures
triangular_aperture = triangle(spatial_range, spatial_range, 3.0)
square_aperture = square(spatial_range, spatial_range, 2.0)

# Shape the beam with apertures
triangular_beam = beam .* triangular_aperture
square_beam = beam .* square_aperture

# Visualize all variants
fig = Figure(size = (1200, 400))
ax1 = Axis(fig[1,1], title="Original vortex")
ax2 = Axis(fig[1,2], title="Triangular aperture")
ax3 = Axis(fig[1,3], title="Square aperture")

heatmap!(ax1, abs2.(beam), colormap = :jet)
heatmap!(ax2, abs2.(triangular_beam), colormap = :jet)
heatmap!(ax3, abs2.(square_beam), colormap = :jet)

fig
```

## Coordinate System and Grid Setup

### Understanding Spatial Grids

StructuredLight.jl uses **Cartesian coordinate systems** where:

- **Coordinates**: (x, y) represent transverse spatial coordinates  
- **Units**: All spatial quantities use consistent units (usually arbitrary units, but can represent physical distances)

### Physical Interpretation

```julia
# Arbitrary units - could represent any length scale
spatial_range = LinRange(-2.5, 2.5, 256)

# Physical units example (e.g., millimeters):
x_mm = LinRange(-5.0, 5.0, 256)  # 10mm total width
```

### Common Conventions

- **Beam waist**: Usually placed at z=0 (focal plane)
- **Propagation direction**: Positive z-axis by convention
- **Angular coordinates**: For `lg()` modes, azimuthal angle measured from positive x-axis

## Key Concepts

- **Spatial grid**: Define your coordinate system using some sort of range to specify the spatial extent and resolution. Here we use `LinRange`
- **Beam modes**: Use `lg()` for Laguerre-Gaussian modes (circular symmetry, orbital angular momentum), `hg()` for Hermite-Gaussian modes (rectangular symmetry), or any other custom modes
- **Propagation**: Use `free_propagation()` to evolve beams through free space
- **Visualization**: The `visualize()` function provides quick plotting, and requires loading `CairoMakie`

## Next Steps

- Explore [Beam Profiles](@ref) for detailed mode generation
- Dive into [Visualization](@ref) for learning some predefined plotting techniques
- Learn about [Propagation](@ref) for advanced beam evolution
- Check [Holograms](@ref) for creating computer-generated holograms in order to produce structured light beams in the lab
- Visit [Phase Modulation](@ref) for optical element simulation
- See the [examples page](examples.md) for complete applications and research reproductions