# Aberration Correction

Optical systems often introduce aberrations that distort beam profiles. StructuredLight.jl provides comprehensive tools for simulating optical elements and correcting aberrations using phase modulation techniques.

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

```julia
using StructuredLight

# Create coordinate grids
x = LinRange(-5, 5, 256)
y = LinRange(-5, 5, 256)

# Generate a Hermite-Gaussian beam
ψ₀ = hg(x, y, m=1, n=1)

# Apply a lens with focal length f=10
f = 10.0
k = 1.0  # wavenumber
lens_phase = lens(x, y, f, f; k=k)

# The focused beam after the lens
ψ_focused = ψ₀ .* lens_phase
```

#### Astigmatic Lens Example

For creating astigmatic focusing with different focal lengths in x and y:

```julia
# Different focal lengths in x and y directions
fx, fy = 8.0, 12.0
astigmatic_lens = lens(x, y, fx, fy; k=k)
ψ_astigmatic = ψ₀ .* astigmatic_lens
```

#### Tilted Lens Example

For simulating a tilted spherical lens:

```julia
# Lens tilted by 15 degrees
f = 10.0
ϕ = π/12  # 15 degrees in radians
tilted_lens_phase = tilted_lens(x, y, f, ϕ; k=k)
ψ_tilted = ψ₀ .* tilted_lens_phase
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
- **(n=1, m=±1)**: Tip/tilt
- **(n=2, m=0)**: Defocus  
- **(n=2, m=±2)**: Astigmatism
- **(n=3, m=±1)**: Coma
- **(n=4, m=0)**: Spherical aberration

### Zernike Examples

```julia
# Create normalized coordinates (Zernike polynomials are defined on unit circle)
ρ = LinRange(-1, 1, 256)

# Defocus aberration (n=2, m=0)
Z₂₀ = zernike_polynomial(ρ, ρ, 0, 2)

# Astigmatism (n=2, m=2)
Z₂₂ = zernike_polynomial(ρ, ρ, 2, 2)

# Coma (n=3, m=1)  
Z₃₁ = zernike_polynomial(ρ, ρ, 1, 3)

# Spherical aberration (n=4, m=0)
Z₄₀ = zernike_polynomial(ρ, ρ, 0, 4)
```

## Aberration Correction

The main aberration correction function allows you to apply arbitrary phase corrections defined as linear combinations of basis functions:

```@docs
aberration_correction!
```

### Practical Aberration Correction

Here's a complete example showing how to correct for multiple aberrations:

```julia
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
ψ_corrected = copy(ψ₀)
aberration_correction!(ψ_corrected, correction_funcs, correction_coeffs, grid)
```

### Custom Aberration Functions

You can also define custom phase correction functions:

```julia
# Custom quadratic phase (similar to a lens)
function custom_lens(args)
    x, y = args
    f = 10.0  # focal length
    return -(x^2 + y^2) / (2f)  # Quadratic phase
end

# Apply custom correction
custom_correction = (custom_lens,)
custom_coeffs = (1.0,)
aberration_correction!(ψ_corrected, custom_correction, custom_coeffs, grid)
```

## Advanced Example: Adaptive Optics Simulation

Here's a more complex example simulating an adaptive optics correction:

```julia
using StructuredLight

# System parameters
λ = 633e-9      # wavelength (m)
k = 2π/λ        # wavenumber
D = 0.01        # aperture diameter (m)
rs = LinRange(-D/2, D/2, 512)

# Generate initial beam with aberrations
ψ₀ = lg(rs, rs, p=0, l=1)  # Vortex beam

# Simulate atmospheric turbulence using multiple Zernike modes
turbulence_modes = [
    (args) -> zernike_polynomial(args..., 0, 2),   # Defocus
    (args) -> zernike_polynomial(args..., 2, 2),   # Astigmatism
    (args) -> zernike_polynomial(args..., 1, 3),   # Coma X
    (args) -> zernike_polynomial(args..., -1, 3),  # Coma Y  
    (args) -> zernike_polynomial(args..., 0, 4),   # Spherical
]

# Random turbulence strengths (in waves RMS)
turbulence_coeffs = 0.1 * λ * randn(5)  # Convert to radians

# Apply turbulence
ψ_turbulent = copy(ψ₀)
aberration_correction!(ψ_turbulent, turbulence_modes, turbulence_coeffs, (rs, rs))

# Measure aberrations (in practice, this would be done with a wavefront sensor)
measured_coeffs = -turbulence_coeffs  # Perfect measurement for demo

# Apply correction
ψ_corrected = copy(ψ_turbulent)
aberration_correction!(ψ_corrected, turbulence_modes, measured_coeffs, (rs, rs))
```

## Tips for Effective Aberration Correction

1. **Coordinate Normalization**: Zernike polynomials are defined on the unit circle. Make sure your coordinate system is properly normalized.

2. **Phase Units**: Correction coefficients are in radians. To convert from waves RMS to radians: `phase_rad = waves_rms * 2π`.

3. **GPU Acceleration**: All aberration correction functions support GPU acceleration via KernelAbstractions.jl backends.

4. **Memory Efficiency**: Use the in-place versions (`lens!`, `aberration_correction!`) when working with large arrays to minimize memory allocation.

5. **Orthogonality**: Zernike polynomials are orthogonal over the unit circle, making them ideal for decomposing arbitrary phase aberrations.

The aberration correction tools in StructuredLight.jl provide a flexible framework for simulating realistic optical systems and implementing phase correction algorithms commonly used in adaptive optics and spatial light modulator applications.