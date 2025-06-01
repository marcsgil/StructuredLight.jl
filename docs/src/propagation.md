# Propagation

The heavy lifting in this package is done by the propagation methods. The basic principle is that you give an initial profile and the grid over which it was calculated, and we give back to you the propagated beam. For now, we have included propagation in free space (`free_propagation`) and in Kerr Media (`kerr_propagation`). The implementation of propagation in media with quadratic nonlinearity is on the horizon.

## `free_propagation`

```@docs
free_propagation
```

### Examples:

```@example
#The simplest usage would be the following:
using StructuredLight, CairoMakie

rs = LinRange(-6,6,256)

#Here we define an initial profile that isn't invariant upon propagation.
ψ₀ = lg(rs,rs) / 2 + lg(rs,rs,p=1)

ψ = free_propagation(ψ₀,rs,rs,1) #Then, we propagate it by a distance z=1.

#Here are the initial profiles the propagated beam, side by side.
visualize([ψ₀,ψ] |> stack .|> abs2)
```

```@example
#We can also provide a collection of z values to produce an animation:
using StructuredLight, CairoMakie

rs = LinRange(-12,12,256)
zs = LinRange(0,1,64)

#This is a laguerre-gaussian modulated by a sine function.
ψ₀ = lg(rs,rs,l=2,p=1) .* map(r->sin(6*r[2]),Iterators.product(rs,rs))

#Now the propagation is performed for each z ∈ zs. The output is a 3D array.
ψs = free_propagation(ψ₀,rs,rs,zs)

save_animation(abs2.(ψs), "lg_times_sin.mp4")
```

![](lg_times_sin.mp4)

```@example
using StructuredLight, CairoMakie

rs = LinRange(-4, 4, 256)
zs = LinRange(0.01, 1, 32)

#This is a gaussian mode
ψ₀ = lg(rs, rs)
scalings = @. √(1 + 4 * zs^2) #Here, we introduce the scalings given by w(z)/w0

#Now we propagate, including the scalings
ψs = free_propagation(ψ₀, rs, rs, zs, scalings)

#Note that the scalings compensate the diffraction of the beam.
#Therefore, the animation seems still.
save_animation(abs2.(ψs), "standing_still.mp4")
```

![](standing_still.mp4)

### References

The book "Schmidt, J. D. (2010). Numerical Simulation of Optical Wave Propagation with Examples in MATLAB. United States: SPIE." is a great resource to learn about numerical propagation of paraxial light beams

## `kerr_propagation`

```@docs
kerr_propagation
```

### Example 
```@example
using StructuredLight, CairoMakie

rs = LinRange(-2.5,2.5,256) #The transverse grid
zs = LinRange(0,.12,32) #The z grid

ψ₀ = lg(rs,rs) #Calculates the fundamental Laguerre-Gaussian mode

#We perform the propagation with a strong nonlinearity
ψ = kerr_propagation(ψ₀,rs,rs,zs,512,g=200)

save_animation(abs2.(ψ), "kerr.mp4") #The beam colapses due to the self focusing effect
```

![](kerr.mp4)


### References

The package [NonlinearSchrodinger.jl](https://github.com/oashour/NonlinearSchrodinger.jl/tree/master) has more available solvers for this equation, but, as far as I can see, it only works with one spatial dimensional. Its author has also written [a paper](https://arxiv.org/abs/2103.14469) that explains the theory that goes behind the numerical solution. I have also written the package [GeneralizedGrossPitaevskii.jl](https://github.com/marcsgil/GeneralizedGrossPitaevskii.jl) that implements the propagation described by a generalized Gross-Pitaevskii equation over arbitrary dimensions.