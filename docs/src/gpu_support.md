# GPU Support

StructuredLight.jl provides GPU acceleration through [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl), enabling high-performance computations for large-scale optical simulations.

StructuredLight.jl provides GPU support through two complementary approaches:

1. **Backend parameter**: Functions like `lg()`, `hg()`, and `zernike_polynomial()` accept a `backend` keyword argument (e.g., `CUDABackend()`) to generate arrays directly on the specified device.

2. **Automatic dispatch**: Functions like `free_propagation()` and `generate_hologram()`, or functions that operate in place automatically detect GPU arrays in their inputs and execute on the same device, providing seamless GPU execution without explicit backend specification.

## GPU-Accelerated Functions

All major computational functions support GPU execution:

- **Beam Profiles**: `hg()`, `lg()`, `diagonal_hg()` - Support `backend` parameter
- **Beam Manipulation**: `grid_linear_combination()` - Automatic GPU dispatch
- **Propagation**: `free_propagation()`, `kerr_propagation()` - Automatic GPU dispatch  
- **Holograms**: `generate_hologram()` - Both `BesselJ1()` and `Simple()` methods with automatic GPU dispatch for in-place version and `backend` parameter for out-of-place version
- **Phase Modulation**: `zernike_polynomial()` - Support `backend` parameter

## Supported Backends

**Fully Tested and Supported:**
- **CUDA**: NVIDIA GPUs - Works flawlessly with comprehensive testing

**Known Issues:**
- **Metal**: Apple Silicon GPUs - Tested but currently not working due to known compatibility issues

**Theoretical Support (Untested):**
- **ROCm**: AMD GPUs - Untested, may work through KernelAbstractions
- **oneAPI**: Intel GPUs - Untested, may work through KernelAbstractions



## Basic CUDA Usage

### Method 1: Using Backend Parameter

```julia
using StructuredLight, CUDA

# Generate beam profile directly on GPU
rs = LinRange(-5, 5, 1024)
ψ₀ = lg(rs, rs, backend=CUDABackend())

# The result is automatically a CuArray
typeof(ψ₀) # CuArray{ComplexF64, 2, ...}
```

### Method 2: Automatic Dispatch

```julia
using StructuredLight, CUDA

rs = LinRange(-5, 5, 1024)
ψ₀ = lg(rs, rs, backend=CUDABackend())

# Propagation automatically runs on GPU via multiple dispatch
zs = LinRange(0, 1, 32)
ψs = free_propagation(ψ₀, rs, rs, zs)
```

## Advanced Examples

### Large-Scale Beam Propagation

```julia
using StructuredLight, CUDA, CairoMakie

# High-resolution grid for detailed simulation
xs = LinRange(-12, 12, 1024)
ys = LinRange(-12, 12, 1024)
zs = LinRange(0, 0.2, 100)

# Generate complex beam pattern on GPU
initial_beam = lg(xs, ys, p=2, l=3, backend=CUDABackend())

# Apply beam shaping
aperture = pupil(xs, ys, 1.0)
shaped_beam = initial_beam .* CuArray(aperture)

# High-resolution propagation 
propagated = free_propagation(shaped_beam, xs, ys, zs)

# Transfer back to CPU for visualization
save_animation(abs2.(Array(propagated)), "gpu_propagation.mp4")
```

### GPU Beam Superposition

```julia
using StructuredLight, CUDA

rs = LinRange(-4, 4, 512)
grid = (rs, rs)

# Define component functions
f1(args) = hg(args..., m=1, n=0)
f2(args) = hg(args..., m=0, n=1) 
f3(args) = lg(args..., p=0, l=1)

funcs = (f1, f2, f3)
coeffs = (0.5, 0.5, 0.7im)

# Create superposition on GPU
superposition = grid_linear_combination(funcs, coeffs, grid, backend=CUDABackend())
```

### GPU Hologram Generation

```julia
using StructuredLight, CUDA

# Create target beam on GPU
x = LinRange(-5, 5, 256)
y = LinRange(-5, 5, 256)
target_beam = lg(x, y, l=2) |> cu

# Generate hologram on GPU - both methods supported
hologram_bessel = generate_hologram(target_beam, 255, 2, 2, BesselJ1())
hologram_simple = generate_hologram(target_beam, 255, 2, 2, Simple())

# Results are CuArrays
typeof(hologram_bessel) # CuArray{UInt8, 2, ...}
```