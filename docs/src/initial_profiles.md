# Initial Profiles

This package implements the Laguerre-Gauss modes (`lg`), the Hermite-Gauss modes (`hg`) and its diagonal version (`diagonal_hg`).

In all cases, one must specify the points or the grids over which the mode is calculated. The other beam parameters are specified through keyword arguments.

## Rotated Hermite-Gauss

```@docs
rotated_hg
```

## Laguerre-Gauss

```@docs
lg
```

### Examples

Fundamental mode a distance z=.5 away from focus:
```@example
using StructuredLight
rs = LinRange(-3,3,256)
ψ₀ = lg(rs,rs,.5)
visualize(abs2.(ψ₀),ratio=2)
```

Analytical propagation of the fundamental mode:
```@example
using StructuredLight
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = lg(rs,rs,zs) 
show_animation(abs2.(ψs),ratio=2)
```

Analytical propagation of a more complicated Laguerre-Gauss:
```@example
using StructuredLight
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = lg(rs,rs,zs,p=1,l=-2)
show_animation(abs2.(ψs),ratio=2)
```

## Hermite-Gauss


```@docs
hg
```

### Examples

Analytical propagation of a Hermite-Gauss with adjusted waist:
```@example
using StructuredLight
rs = LinRange(-6,6,256) 
zs = LinRange(0,.5,32)
ψs = hg(rs,rs,zs,m=1,n=2,w=.5)
show_animation(abs2.(ψs),ratio=2)
```

Superposition of Hermite-Gausss. Note that we can propagate backwards using negative values of `z`.
```@example
using StructuredLight
rs = LinRange(-4,4,256) 
zs = LinRange(-.5,.7,128)
ψ₀ = hg(rs,rs) + hg(rs,rs,n=2)/√2
ψs = free_propagation(ψ₀,rs,rs,zs)
show_animation(abs2.(ψs),ratio=2)
```

## Diagonal Hermite-Gauss

```@docs
diagonal_hg
```

### Examples

Analytical propagation of a Hermite-Gauss with adjusted wavenumber:
```@example
using StructuredLight
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = diagonal_hg(rs,rs,zs,m=4,n=3,k=2)
show_animation(abs2.(ψs),ratio=2)
```

## Lenses

We also implement lenses:

```@docs
lens
tilted_lens
```