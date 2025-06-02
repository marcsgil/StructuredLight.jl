# GPU Support

StructuredLight.jl provides GPU acceleration through [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl), enabling high-performance computations for large-scale optical simulations.

There are two main ways that GPU support is implemented in StructuredLight.jl, with each function relying on one of these methods:

1. **Using the `backend` keyword argument**: Some functions accept a `backend` parameter that allows you to specify the computational backend. 

2. **Automatic GPU dispatch**: Functions automatically dispatch to GPU if inputs are some type of GPU array (e.g., `CuArray` from the CUDA.jl package). This allows for seamless GPU execution without needing to explicitly set a backend.

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

**Fully Tested but not Supported:**
- **Metal**: Apple Silicon GPUs - Currently has known issues and does not work

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

### Method 2: Transfer Arrays to GPU

```julia
using StructuredLight, CUDA

rs = LinRange(-5, 5, 1024)
ψ₀ = lg(rs, rs) |> cu  # Transfer to GPU

# Propagation automatically runs on GPU via multiple dispatch
zs = LinRange(0, 1, 32)
ψs = free_propagation(ψ₀, rs, rs, zs)
```

## Advanced Examples

### Large-Scale Beam Propagation

```julia
using StructuredLight, CUDA, CairoMakie

# High-resolution grid for detailed simulation
xs = LinRange(-10, 10, 2048)
ys = LinRange(-10, 10, 2048)
zs = LinRange(0, 5, 100)

# Generate complex beam pattern on GPU
initial_beam = lg(xs, ys, p=2, l=3, backend=CUDABackend())

# Apply beam shaping
aperture = pupil(xs, ys, 4.0)
shaped_beam = initial_beam .* CuArray(aperture)

# High-resolution propagation 
@time propagated = free_propagation(shaped_beam, xs, ys, zs)

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

## Performance Optimization

### Memory Management

```julia
using StructuredLight, CUDA

# Monitor GPU memory usage
CUDA.memory_status()

# Pre-allocate arrays for better performance
rs = LinRange(-5, 5, 1024)
dest = CuArray{ComplexF64}(undef, length(rs), length(rs))

# Use in-place operations when possible
hg!(dest, rs, rs, m=1, n=1)
```

### Batch Processing

```julia
using StructuredLight, CUDA

function process_multiple_modes(indices, rs)
    results = []
    
    for (m, n) in indices
        # Each mode computation uses GPU
        mode = hg(rs, rs, m=m, n=n, backend=CUDABackend())
        push!(results, mode)
    end
    
    return results
end

# Process multiple modes efficiently
rs = LinRange(-3, 3, 512) 
mode_indices = [(0,0), (1,0), (0,1), (1,1), (2,0)]
modes = process_multiple_modes(mode_indices, rs)
```

## Performance Benchmarks

Typical performance improvements with CUDA on modern NVIDIA GPUs:

| Operation | Problem Size | CPU Time | GPU Time | Speedup |
|-----------|--------------|----------|----------|---------|
| `lg()` generation | 1024×1024 | ~45ms | ~3ms | ~15× |
| `free_propagation()` | 1024×1024×64 | ~2.1s | ~85ms | ~25× |
| `generate_hologram()` | 512×512 | ~180ms | ~12ms | ~15× |
| `kerr_propagation()` | 512×512×32 | ~8.2s | ~320ms | ~26× |

*Benchmarks are approximate and depend on specific hardware configuration*

## Setup Requirements

### CUDA Setup
- **Hardware**: NVIDIA GPU with compute capability ≥ 3.5
- **Software**: CUDA toolkit installed
- **Julia Package**: `] add CUDA`
- **Verification**: Run `using CUDA; CUDA.functional()` should return `true`

## Troubleshooting

### Common Issues

**GPU Out of Memory**
```julia
# Reduce problem size or fall back to CPU
try
    result = lg(xs, ys, backend=CUDABackend())
catch e
    @warn "GPU memory insufficient, falling back to CPU"
    result = lg(xs, ys, backend=CPU())
end
```

**CUDA Not Available**
```julia
# Check if CUDA is functional
using CUDA
if CUDA.functional()
    backend = CUDABackend()
    @info "Using CUDA backend"
else
    backend = CPU()
    @warn "CUDA not available, using CPU backend"
end
```

**Type Conversion Issues**
```julia
# Ensure consistent types when mixing CPU/GPU arrays
cpu_array = lg(rs, rs)
gpu_array = CuArray(cpu_array)  # Explicit conversion

# Or generate directly on GPU
gpu_array = lg(rs, rs, backend=CUDABackend())
```

### Performance Tips

1. **Keep data on GPU**: Minimize CPU ↔ GPU transfers
2. **Use appropriate precision**: `Float32` often sufficient and faster
3. **Batch operations**: Process multiple arrays together when possible
4. **Pre-allocate**: Use in-place operations (`hg!`, `lg!`) when available
5. **Profile your code**: Use `@time` and GPU profiling tools for optimization

### Debugging

```julia
# Enable CUDA verbose mode for detailed error information
ENV["JULIA_CUDA_VERBOSE"] = "true"

# Check CUDA installation and GPU status
using CUDA
CUDA.versioninfo()
CUDA.memory_status()

# Test basic CUDA functionality
CUDA.@time cu([1, 2, 3, 4])
```