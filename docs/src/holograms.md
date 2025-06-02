# Holograms

This section covers the generation of computer-generated holograms (CGHs) for spatial light modulators (SLMs). The StructuredLight.jl package provides algorithms to convert complex beam profiles into phase patterns that can be displayed on phase-only SLMs. In order to display the holograms in the SLM, we suggest the package [SpatialLightModulator.jl](https://github.com/marcsgil/SpatialLightModulator.jl)

## Overview

Computer-generated holograms encode both amplitude and phase information of a desired optical field into a phase-only pattern. This is essential for spatial light modulators that can only modulate the phase of light.

The basic idea is to take a complex field (which can be a superposition of different modes) and convert it into a phase pattern that can be displayed on an SLM. The hologram generation process involves calculating the phase shift required for each pixel to reconstruct the desired field at a specific distance.

Here's a simple example of how to generate a hologram for a Laguerre-Gaussian beam:
```@example
using StructuredLight, CairoMakie

# Create coordinate arrays
# Should have the resolution of the SLM
x = LinRange(-5, 5, 512)
y = LinRange(-5, 5, 512)

# Generate a Laguerre-Gaussian beam with topological charge l=2
beam = lg(x, y, l=2)

# Create hologram with the default method (BesselJ1)
# We specify the two_pi_modulation, x_period, and y_period
hologram = generate_hologram(beam, 255, 20, 20)

# Visualize the hologram
# For displaying in the SLM, you would typically use a package like SpatialLightModulator.jl
visualize(hologram, colormap = :grays)
```

One can see the distinct fork structure of the vortex beam in the hologram.

## Description

The main functions for hologram generation are:

```@docs
generate_hologram
generate_hologram!
```

These functions take a relative field, as well as extra parameters and output phase pattern suitable for SLMs. The relative field is the desired field divided by the input field at the hologram plane. When the input field is a large Gaussian beam, we can approximate it as a plane wave, and the relative field becomes simply the desired field.

The other parameters are:

- **`two_pi_modulation`**: Integer value between 0 and 255 that produces 2π phase modulation (check your SLM specifications)
- **`x_period`, `y_period`**: Period, in pixels, of the diffraction grating for the x and y directions


The package implements two well-established hologram generation methods:

- **BesselJ1 Method**: Based on the inverse of the Bessel J₁ function, providing high-quality amplitude encoding
```@docs
BesselJ1
```
- **Simple Method**: A straightforward amplitude modulation approach with direct phase encoding
```@docs
Simple
```

### Example: Hermite-Gaussian Mode Hologram

```@example
using StructuredLight, CairoMakie

# Create coordinate arrays
# Should have the resolution of the SLM
x = LinRange(-5, 5, 512)
y = LinRange(-5, 5, 512)

# Generate a Hermite-Gaussian beam with m=2, n=2
beam = hg(x, y, m=2, n=2)

# Create hologram with the Simple() method
# We specify the two_pi_modulation, x_period, and y_period
hologram = generate_hologram(beam, 255, 20, 20, Simple())

# Visualize the hologram
visualize(hologram, colormap = :grays)
```

### Creating Superposition States

```@example
using StructuredLight, CairoMakie

# Create superposition of LG modes
x = LinRange(-8, 8, 512)
y = LinRange(-8, 8, 512)

# Individual modes
lg_0_1 = lg(x, y, p=0, l=1, w=3.0)
lg_0_minus1 = lg(x, y, p=0, l=-1, w=3.0)

# Superposition (equal weights)
superposition = lg_0_1 + 0.5 * lg_0_minus1

# Generate hologram
super_hologram = generate_hologram(superposition, 255, 30, 30)
visualize(super_hologram, colormap=:grays)
```

### Gaussian input
```@example
using StructuredLight, CairoMakie

# Create coordinate arrays
x = LinRange(-5, 5, 512)
y = LinRange(-5, 5, 512)

# We wish to produce a Laguerre-Gaussian beam with l=3 
desired = lg(x, y, l = 3)

# But the input is a gaussian beam
input = lg(x, y, w=2)

# Define the relative field
relative_field = desired ./ input

# Generate hologram for the relative field
hologram_relative = generate_hologram(relative_field, 255, 20, 20)

visualize(hologram_relative, colormap = :grays)
```

### GPU-Accelerated Hologram Generation

```julia
using StructuredLight, CUDA

# Create beam on GPU
# We use Float32 for better performance on GPU
x = LinRange(-5f0, 5f0, 1024)
y = LinRange(-5f0, 5f0, 1024)
beam_gpu = lg(x, y, p=0, l=3, w=2.5, backend=CUDABackend())

# Generate hologram on GPU
hologram_gpu = generate_hologram(beam_gpu, 255, 20, 20)

# Transfer back to CPU if needed
hologram_cpu = Array(hologram_gpu)
```