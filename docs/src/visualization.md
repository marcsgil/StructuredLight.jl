# Visualization

## `visualize` 

This package offer tools that help the visualization of the generated beams. The simplest of such tools is the `visualize` function, which displays the beam as a static image. It accepts two, three and four dimensional arrays.

```@docs
visualize
```

Example:

```julia
using StructuredLight

rs = LinRange(-4,4,256) 

ψ₁ = stack( [diagonal_hg(rs,rs,m=3), diagonal_hg(rs,rs,n=3)] )

ψ₂ = stack( [hg(rs,rs,m=3), hg(rs,rs,n=3)] )

visualize(ψ₁)

visualize(ψ₂)

ψ₃ = stack([ψ₁,ψ₂])

visualize(ψ₃)
``` 

 With three-dimensional arrays, which are interpreted as a sequence of 2D images, we can also form different kinds of animations. 

## `show_animation`

`show_animation` displays a gif of the sequence of images. It will be automatically displayed on VS Code.

```@docs
show_animation
```

Example:

```julia
using StructuredLight

rs = LinRange(-6,6,256) 
zs = LinRange(0,1,32)

ψ₀ = hg(rs,rs) + hg(rs,rs,n=2)/√2

ψs = free_propagation(ψ₀,rs,rs,zs)

show_animation(ψs,fps=12)
```

## `save_animation`

`save_animation` will save the animation on the given path.

```@docs
save_animation
```

Example:

```julia
using StructuredLight

rs = LinRange(-6,6,256) 
zs = LinRange(0,1,32)

ψ₀ = hg(rs,rs) + hg(rs,rs,n=2)/√2

ψs = free_propagation(ψ₀,rs,rs,zs)

save_animation(ψs,"animation.gif",fps=12)
save_animation(ψs,"animation.mp4",fps=12)
```