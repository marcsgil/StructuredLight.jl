# StructuredLight

This is a package that simulates the propagation of a paraxial light beam. Beyond that, it also defines a few important initial profiles, that frequently show up when working with structured light.

It is not yet registered, but you can use it by hitting `]` to enter the Pkg REPL-mode and then typying 

```
add https://github.com/marcsgil/StructuredLight.jl
```
# Example

The following code is a minimal working exemple for this package:

```julia
using FreeParaxialPropagation #Loads the package

rs = LinRange(-5,5,256) #Define a linear grid of points

E0 = lg(rs,rs) #Calculates the fundamental Laguerre-Gaussian mode

vizualize(E0) #Vizualizes the mode

E = free_propagation(E0,rs,rs,1) #Propagates the mode through a distance of z=1

vizualize(E) #Vizualizes the evolved mode
```

This ilustrates the basic idea of the package: first, you construct a matrix representing the mode you want to propagate, and then one calls `free_propagation`.

One can also call `free_propagation` with a vector of z values:

```julia
zs = LinRange(0,1,32)

E = free_propagation(E0,rs,rs,zs)

show_animation(E) #Produces an animation representing the propagation.
```

## CUDA support

The propagation can be done on Nvidia GPUs. One only needs to pass the initial profile as a `CuArray`.

```julia
using CUDA

E0 = lg(rs,rs) |> cu #Transfers array to gpu

E = free_propagation(E0,rs,rs,zs)
```

# Initial Profiles

This package implements the Laguerre-Gaussian modes (`lg`), the Hermite-Gaussian modes (`hg`) and its diagonal version (`diagonal_hg`). A detailed explanation of all the possible parameters can be accessed by entering in the help mode of the REPL (by hitting `?`) and then typing `lg`,`hg` or `diagonal_hg`.

Here are a few examples:

```julia
rs = LinRange(-5,5,256)
E0 = lg(rs,rs,.5) # Fundamental mode a distance z=.5 away from focus
vizualize(E0)
```

```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
Es = lg(rs,rs,zs) #Analytical propagation of the fundamental mode
show_animation(Es)
```

```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
Es = lg(rs,rs,zs,p=1,l=-2) #Analytical propagation of a more complicated laguerre-gaussian
show_animation(Es)
```

```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
Es = hg(rs,rs,zs,m=1,n=2,w0=.5) #Analytical propagation of a hermite-gaussian with adjusted waist
show_animation(Es)
```

```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
Es = diagonal_hg(rs,rs,zs,m=4,n=3,k=2) #Analytical propagation of a hermite-gaussian with adjusted wavenumber
show_animation(Es)
```

```julia
rs = LinRange(-4,4,256) 
zs = LinRange(-.5,.7,128) #We can propagate backwards using negative values of z
E0 = hg(rs,rs) + hg(rs,rs,n=2)/√2  #Superposition of hermite-gaussians
Es = free_propagation(E0,rs,rs,zs)
show_animation(Es)
```
