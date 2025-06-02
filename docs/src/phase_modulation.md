# Phase Modulation

The phase modulation tools in StructuredLight.jl provide a flexible framework for simulating optical elements, correcting aberrations, and implementing custom phase patterns commonly used in adaptive optics, spatial light modulator applications, and structured light generation. This includes simulating optical elements like lenses, correcting aberrations, and applying custom phase masks. The phase modulation framework is general-purpose and can be used for any application requiring spatially-varying phase modifications.

## Optical Elements

### Lenses

The package provides functions to simulate the phase modulation introduced by lenses:

```@docs
lens
lens!
tilted_lens
tilted_lens!
```

#### Basic Lens Example

```@example
using StructuredLight, CairoMakie

# Create coordinate grids
x = LinRange(-4, 4, 512)
y = LinRange(-4, 4, 512)
z = LinRange(0, 0.5, 128)

# Generate a Hermite-Gaussian beam
ψ₀ = hg(x, y, m=1, n=1)

# Apply a lens with focal length f=10
f = 0.5
lens_phase = lens(x, y, f, f)

# The beam is focused after propagation through the lens
ψ_focused = free_propagation(ψ₀ .* lens_phase, x, y, z)

save_animation(abs2.(ψ_focused), "focused_beam.mp4")
return nothing # hide
```
![](focused_beam.mp4)

#### Astigmatic Lens Example

For creating astigmatic focusing with different focal lengths in x and y:

```@example
using StructuredLight # hide
x = LinRange(-4, 4, 256) # hide
y = LinRange(-4, 4, 256) # hide
ψ₀ = hg(x, y, m=1, n=1) # hide
k = 1.0 # hide

# Different focal lengths in x and y directions
fx, fy = 8.0, 12.0
astigmatic_lens = lens(x, y, fx, fy; k=k)
ψ_astigmatic = ψ₀ .* astigmatic_lens

nothing # hide
```

#### Tilted Lens Example

For simulating a tilted spherical lens:

```@example
using StructuredLight # hide
x = LinRange(-4, 4, 256) # hide
y = LinRange(-4, 4, 256) # hide
ψ₀ = hg(x, y, m=1, n=1) # hide
k = 1.0 # hide

# Lens tilted by 15 degrees
f = 10.0
ϕ = π/12  # 15 degrees in radians
tilted_lens_phase = tilted_lens(x, y, f, ϕ; k=k)
ψ_tilted = ψ₀ .* tilted_lens_phase

nothing # hide
```

## Zernike Polynomials

Zernike polynomials are the standard basis for describing optical aberrations. They form a complete orthogonal set over the unit circle.

```@docs
zernike_polynomial
zernike_polynomial!
```

### Understanding Zernike Indices

Zernike polynomials are characterized by two indices:
- **n**: radial order (non-negative integer)
- **m**: azimuthal frequency (integer with |m| ≤ n and n-|m| even)

Common aberrations correspond to specific Zernike modes:
- **(n=1, m=±1)**: Tilt
- **(n=2, m=0)**: Defocus  
- **(n=2, m=±2)**: Astigmatism
- **(n=3, m=±1)**: Coma
- **(n=4, m=0)**: Spherical aberration

### Zernike Examples

```@example
using StructuredLight, CairoMakie

# Create coordinates
ρ = LinRange(-1, 1, 256)

function reject_outside_disk!(dest, ρ)
    for n ∈ axes(dest, 2), m ∈ axes(dest, 1)
        if ρ[n]^2 + ρ[m]^2 > 1
            dest[n, m] = NaN  # Set outside disk to NaN for visualization
        end
    end
end

# Defocus aberration (n=2, m=0)
Z₂₀ = zernike_polynomial(ρ, ρ, 0, 2)

# Astigmatism (n=2, m=2)
Z₂₂ = zernike_polynomial(ρ, ρ, 2, 2)

# Coma (n=3, m=1)
Z₃₁ = zernike_polynomial(ρ, ρ, 1, 3)

# Spherical aberration (n=4, m=0)
Z₄₀ = zernike_polynomial(ρ, ρ, 0, 4)

for Z in [Z₂₀, Z₂₂, Z₃₁, Z₄₀]
    reject_outside_disk!(Z, ρ)
end

# Visualize the Zernike polynomials
fig = Figure(size=(800, 750))
ax1 = Axis(fig[1, 1], title="Defocus (Z₂₀)", aspect=1)
ax2 = Axis(fig[1, 2], title="Astigmatism (Z₂₂)", aspect=1)
ax3 = Axis(fig[2, 1], title="Coma (Z₃₁)", aspect=1)
ax4 = Axis(fig[2, 2], title="Spherical (Z₄₀)", aspect=1)

heatmap!(ax1, ρ, ρ, Z₂₀, colormap=:bluesreds)
heatmap!(ax2, ρ, ρ, Z₂₂, colormap=:bluesreds)
heatmap!(ax3, ρ, ρ, Z₃₁, colormap=:bluesreds)
heatmap!(ax4, ρ, ρ, Z₄₀, colormap=:bluesreds)

fig
```

## General Phase Application

The main phase modulation function allows you to apply arbitrary phase corrections defined as linear combinations of basis functions:

```@docs
apply_phase!
```


### Practical Examples

This is the same example as above, but using the `apply_phase!` to update a lens phase in place:

```@example
using StructuredLight, CairoMakie

# Create coordinate grids
x = LinRange(-4, 4, 512)
y = LinRange(-4, 4, 512)
z = LinRange(0, 0.5, 128)

# Generate a Hermite-Gaussian beam
ψ₀ = hg(x, y, m=1, n=1)

# Apply a lens with focal length f=10
funcs = (args -> - (args[1]^2 / 0.5 + args[2]^2 / 0.5) / 2,)
coeffs = (1,)

apply_phase!(ψ₀, funcs, coeffs, (x, y))

# The beam is focused after propagation through the lens
ψ_focused = free_propagation(ψ₀, x, y, z)

save_animation(abs2.(ψ_focused), "focused_beam2.mp4")

nothing # hide
```

![](focused_beam2.mp4)

Here's a complete example showing how to apply multiple phase corrections:

```@example
using StructuredLight

# Create coordinate system
rs = LinRange(-3, 3, 256)
grid = (rs, rs)

# Generate initial beam
ψ₀ = lg(rs, rs, p=1, l=2)  # Laguerre-Gaussian beam

# Define aberration correction functions using Zernike polynomials
# Each function takes a tuple of coordinates (x, y)
defocus(args) = zernike_polynomial(args..., 0, 2)     # Defocus
astigmatism(args) = zernike_polynomial(args..., 2, 2) # Astigmatism  
coma(args) = zernike_polynomial(args..., 1, 3)        # Coma

# Correction coefficients (in radians of phase)
correction_funcs = (defocus, astigmatism, coma)
correction_coeffs = (0.5, -0.3, 0.2)  # Phase corrections in radians

# Apply corrections to the beam
apply_phase!(ψ₀, correction_funcs, correction_coeffs, grid)

nothing # hide
```

### Custom Phase Functions

You can also define custom phase modulation functions:

```@example
using StructuredLight # hide
rs = LinRange(-3, 3, 256) # hide
grid = (rs, rs) # hide
ψ_corrected = lg(rs, rs, p=1, l=2) # hide

# Custom quadratic phase (similar to a lens)
function custom_lens(args)
    x, y = args
    f = 10.0  # focal length
    return -(x^2 + y^2) / (2f)  # Quadratic phase
end

# Apply custom phase modulation
custom_correction = (custom_lens,)
custom_coeffs = (1.0,)
apply_phase!(ψ_corrected, custom_correction, custom_coeffs, grid)

nothing # hide
```

## Advanced Example: Adaptive Optics Simulation

Here's a more complex example simulating an adaptive optics correction:

```@example
using StructuredLight

# System parameters
D = 0.01        # aperture diameter (m)
rs = LinRange(-D/2, D/2, 512)

# Generate initial beam with aberrations
ψ₀ = lg(rs, rs, p=0, l=1)  # Vortex beam

# Simulate atmospheric turbulence using multiple Zernike modes
turbulence_modes = [
    (args) -> zernike_polynomial((args ./ D)..., 0, 2),   # Defocus
    (args) -> zernike_polynomial((args ./ D)..., 2, 2),   # Astigmatism
    (args) -> zernike_polynomial((args ./ D)..., 1, 3),   # Coma X
    (args) -> zernike_polynomial((args ./ D)..., -1, 3),  # Coma Y  
    (args) -> zernike_polynomial((args ./ D)..., 0, 4),   # Spherical
]

# Random turbulence strengths
turbulence_coeffs = 0.1 .* ntuple(n->rand(), 5) 

# Apply turbulence
ψ_turbulent = copy(ψ₀)
apply_phase!(ψ_turbulent, turbulence_modes, turbulence_coeffs, (rs, rs))

# Measure aberrations (in practice, this would be done with a wavefront sensor)
measured_coeffs = .-turbulence_coeffs  # Perfect measurement for demo

# Apply correction
ψ_corrected = copy(ψ_turbulent)
apply_phase!(ψ_corrected, turbulence_modes, measured_coeffs, (rs, rs))

nothing # hide
```

