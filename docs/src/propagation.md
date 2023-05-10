# Propagation

The heavy lifting in this package is done by the propagation methods. The basic principle is that you give an initial profile and the grid over which it was calculated, and we give back to you the propagated beam. For now, we have included propagation in free space (`free_propagation`) and in Kerr Media (`kerr_propagation`). The implementation of propagation in media with quadratic nonlinearity is on the horizon.

```@docs
free_propagation
```

Examples:

```julia
#The simplest usage would be the following:
using StructuredLight

rs = LinRange(-6,6,256)

#Here we define an initial profile that isn't invariant upon propagation.
ψ₀ = lg(rs,rs) + lg(rs,rs,p=1) 

ψ = free_propagation(ψ₀,rs,rs,1) #Then, we propagate it by a distance z=1.

visualize(ψ₀) # Here is the initial beam.
visualize(ψ) # And here is the propagated beam.
```

```julia
#We can also provide a collection of z values to produce an animation:
using StructuredLight

rs = LinRange(-6,6,256)
zs = LinRange(0,1/2,32)

#This is a gaussian modulated by a sine function.
ψ₀ = lg(rs,rs) .* map(r->sin(6*r[2]),Iterators.product(rs,rs))

#Now the propagation is performed for each z ∈ zs. The output is a 3D array.
ψs = free_propagation(ψ₀,rs,rs,zs)

show_animation(ψs) # Here is the initial beam.
```

```@docs
kerr_propagation
```

One can also call `free_propagation` with a vector of z values:

