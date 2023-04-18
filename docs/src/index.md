# Introduction

This package provides tools to simulate the propagation of paraxial light beams. This includes the calculation of Laguerre-Gauss and Hermite-Gauss beam profiles, the action of lenses, the propagation in free space as well as in Kerr media. We also provide methods that help the visualization of such beams.

The package is not yet registered, but you can use it by hitting `]` to enter the Pkg REPL-mode and then typying

```
add https://github.com/marcsgil/StructuredLight.jl
```
## Example

The following code is a minimal working exemple for this package:

```julia
using StructuredLight #Loads the package

rs = LinRange(-5,5,256) #Define a linear grid of points

E0 = lg(rs,rs) #Calculates the fundamental Laguerre-Gaussian mode

visualize(E0) #visualizes the mode

E = free_propagation(E0,rs,rs,1) #Propagates the mode through a distance of z=1

visualize(E) #visualizes the evolved mode
```

This ilustrates the basic idea of the package: first, you construct a matrix representing the mode you want to propagate, and then one calls `free_propagation`.

## CUDA support

The propagation can be done on Nvidia GPUs. One only needs to pass the initial profile as a `CuArray`.

```julia
using CUDA

E0 = lg(rs,rs) |> cu #Transfers array to gpu

E = free_propagation(E0,rs,rs,zs)
```

