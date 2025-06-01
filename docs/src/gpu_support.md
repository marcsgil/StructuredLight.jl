# GPU Support

Both `free_propagation` and `kerr_propagation` can be run on Nvidia GPUs, which will greatly improve the performance of these functions. If you have one, you simply need to convert your initial profile to a `CuArray` and pass this converted array to the propagation methods (check the [CUDA.jl documentation](https://cuda.@examplegpu.org/stable/) for more details). Then, [multiple dispatch](https://docs.@examplelang.org/en/v1/manual/methods/#Methods) will do its magic!

Here is an example:
```julia
using StructuredLight, CairoMakie
using CUDA #It is necessary to load the CUDA package

ψ₀ = lg(rs,rs) |> cu #Transfers array to GPU

ψ = free_propagation(ψ₀,rs,rs,zs) #This is running on the GPU!
```