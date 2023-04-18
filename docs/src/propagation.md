# Propagation

One can also call `free_propagation` with a vector of z values:

```julia
zs = LinRange(0,1,32)

E = free_propagation(E0,rs,rs,zs)

show_animation(E) #Produces an animation representing the propagation.
```