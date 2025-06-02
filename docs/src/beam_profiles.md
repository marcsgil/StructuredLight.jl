# Beam Profiles

This package implements the Laguerre-Gauss modes ([`lg`](@ref)), the Hermite-Gauss modes ([`hg`](@ref)) and its diagonal version ([`diagonal_hg`](@ref)).

In all cases, one must specify the points or the grids over which the mode is calculated. The other beam parameters are specified through keyword arguments.

## Laguerre-Gauss

```@docs
lg
normalization_lg
```

### Examples

Fundamental mode a distance z=.5 away from focus:
```@example
using StructuredLight, CairoMakie
rs = LinRange(-3,3,256)
ψ₀ = lg(rs,rs,.5)
visualize(abs2.(ψ₀))
```

Analytical propagation of the fundamental mode:
```@example
using StructuredLight, CairoMakie
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = lg(rs,rs,zs) 
save_animation(abs2.(ψs),"prop_fundamental.mp4")
return nothing # hide
```

![](prop_fundamental.mp4)

Analytical propagation of a more complicated Laguerre-Gauss:
```@example
using StructuredLight, CairoMakie
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = lg(rs,rs,zs,p=1,l=-2)
save_animation(abs2.(ψs),"prop_p=1,l=-2.mp4")
return nothing # hide
```

![](prop_p=1,l=-2.mp4)

## Hermite-Gauss


```@docs
hg
normalization_hg
```

### Examples

Analytical propagation of a Hermite-Gauss with adjusted waist:
```@example
using StructuredLight, CairoMakie
rs = LinRange(-6,6,256) 
zs = LinRange(0,.5,32)
ψs = hg(rs,rs,zs,m=1,n=2,w=.4)
save_animation(abs2.(ψs),"prop_adjusted_waist.mp4")
```

![](prop_adjusted_waist.mp4)

Superposition of Hermite-Gausss. Note that we can propagate backwards using negative values of `z`.
```@example
using StructuredLight, CairoMakie
rs = LinRange(-4,4,256) 
zs = LinRange(-.5,.7,128)
ψ₀ = hg(rs,rs) + hg(rs,rs,n=2)/√2
ψs = free_propagation(ψ₀,rs,rs,zs)
save_animation(abs2.(ψs),"prop_backwards.mp4")
```

![](prop_backwards.mp4)

## Diagonal Hermite-Gauss

```@docs
diagonal_hg
```

### Examples

Analytical propagation of a Hermite-Gauss with adjusted wavenumber:
```@example
using StructuredLight, CairoMakie
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = diagonal_hg(rs,rs,zs,m=4,n=3,k=2)
save_animation(abs2.(ψs),"prop_adjusted_k.mp4")
```

![](prop_adjusted_k.mp4)

## Aperture Functions & Beam Shaping

StructuredLight.jl provides several aperture functions for beam shaping and diffraction studies. These functions return boolean arrays indicating which points are within the aperture geometry.

```@docs
rectangular_apperture
square
single_slit
double_slit
pupil
triangle
```

### Examples

Basic rectangular aperture:
```@example
using StructuredLight, CairoMakie
xs = LinRange(-2, 2, 256)
ys = LinRange(-2, 2, 256)
aperture = rectangular_apperture(xs, ys, 1.5, 1.0)
visualize(aperture)
```

Single slit diffraction setup:
```@example
using StructuredLight, CairoMakie
xs = LinRange(-3, 3, 256)
ys = LinRange(-3, 3, 256)
slit = single_slit(xs, ys, 0.5)
visualize(slit)
```

Double slit interference:
```@example
using StructuredLight, CairoMakie
xs = LinRange(-4, 4, 256)
ys = LinRange(-4, 4, 256)
double_slits = double_slit(xs, ys, 0.3, 2.0)
visualize(double_slits)
```

Circular pupil (common in optical systems):
```@example
using StructuredLight, CairoMakie
xs = LinRange(-2, 2, 256)
ys = LinRange(-2, 2, 256)
circular_aperture = pupil(xs, ys, 1.0)
visualize(circular_aperture)
```

Beam shaping with apertures - applying aperture to Gaussian beam:
```@example
using StructuredLight, CairoMakie
xs = LinRange(-3, 3, 256)
ys = LinRange(-3, 3, 256)
gaussian_beam = hg(xs, ys)
square_aperture = square(xs, ys, 2.0)
shaped_beam = gaussian_beam .* square_aperture
visualize(abs2.(shaped_beam))
```

## Advanced Beam Manipulation

For creating complex beam patterns through superposition of multiple functions:

```@docs
linear_combination
grid_linear_combination
grid_linear_combination!
```

### Examples

Creating a superposition of Hermite-Gaussian modes:
```@example
using StructuredLight, CairoMakie
rs = LinRange(-3, 3, 100)
grid = (rs, rs)

# Define component functions
f1(args) = hg(args..., m=1)
f2(args) = hg(args..., n=1)

funcs = (f1, f2)
coeffs = (1/√2, im/√2)

# This creates a Laguerre-Gaussian mode through HG superposition
superposition = grid_linear_combination(funcs, coeffs, grid)
visualize(abs2.(superposition), scaling=4)
```

Creating custom beam patterns with multiple components:
```@example
using StructuredLight, CairoMakie
rs = LinRange(-4, 4, 128)
grid = (rs, rs)

# Define multiple beam components
f1(args) = lg(args..., p=0, l=1)  # Vortex beam
f2(args) = hg(args..., m=2, n=0)  # Hermite-Gaussian
f3(args) = hg(args..., m=0, n=2)  # Another HG mode

funcs = (f1, f2, f3)
coeffs = (0.6, 0.3, 0.3)

complex_beam = grid_linear_combination(funcs, coeffs, grid)
visualize(abs2.(complex_beam), scaling=4)
```