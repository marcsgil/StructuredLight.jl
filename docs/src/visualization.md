# Visualization

This package offer tools that help the visualization of the generated beams. They are built on top of the [CairoMakie](https://docs.makie.org/stable/) package, which is a powerful plotting library in Julia. It needs to be installed separately and loaded with `using CairoMakie` before using the visualization functions.

## `visualize` 

The simplest of such tools is the `visualize` function, which displays the beam as a static image. It accepts two, three and four dimensional arrays. The first two dimensions are interpreted as the spatial coordinates, the third dimension is interpreted as different rows of images, and the fourth dimension is interpreted as different columns of images.

```@docs
visualize
```

### Example:

```@example
using StructuredLight, CairoMakie

rs = LinRange(-4,4,256) 

ψ₁ = stack( [diagonal_hg(rs,rs,m=3), diagonal_hg(rs,rs,n=3)] )

ψ₂ = stack( [hg(rs,rs,m=3), hg(rs,rs,n=3)] )

visualize(abs2.(ψ₁)) #Displays ψ₁ 

visualize(abs2.(ψ₂)) #Displays ψ₂

ψ₃ = stack([ψ₁,ψ₂])

visualize(abs2.(ψ₃)) #Displays ψ₁ and ψ₂ in a row.
``` 

## `save_animation`

With three-dimensional arrays, which are interpreted as a sequence of 2D images, we can also form animations. `save_animation` will save the animation on the given path. This is useful for visualizing the evolution of beams as they propagate through space or interact with different apertures.

```@docs
save_animation
```

### Example:

```@example
using StructuredLight, CairoMakie

rs = LinRange(-6,6,256) 
zs = LinRange(0,1,32)

ψ₀ = hg(rs,rs) + hg(rs,rs,n=2)/√2

ψs = free_propagation(ψ₀,rs,rs,zs)

save_animation(abs2.(ψs),"animation.mp4",framerate=12)

nothing # hide
```

![](animation.mp4)