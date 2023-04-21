# Propagation

One can also call `free_propagation` with a vector of z values:

```julia
zs = LinRange(0,1,32)

ψ = free_propagation(ψ₀,rs,rs,zs)

show_animation(ψ) #Produces an animation representing the propagation.
```