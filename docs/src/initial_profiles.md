# Initial Profiles

This package implements the Laguerre-Gaussian modes (`lg`), the Hermite-Gaussian modes (`hg`) and its diagonal version (`diagonal_hg`).

In all cases, one must specify the grids `xs` and `ys` over which the mode is calculated. One may also specify a distance `z` away from the focus, which defaults to `0`, or a collection of distances. The other beam parameters are specified through keyword arguments.

## Laguerre-Gaussian

```@docs
lg
```

### Examples

Fundamental mode a distance z=.5 away from focus:
```julia
rs = LinRange(-5,5,256)
ψ₀ = lg(rs,rs,.5)
visualize(ψ₀)
```

Analytical propagation of the fundamental mode:
```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = lg(rs,rs,zs) 
show_animation(ψs)
```

Analytical propagation of a more complicated laguerre-gaussian:
```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = lg(rs,rs,zs,p=1,l=-2)
show_animation(ψs)
```

## Hermite-Gaussian


```@docs
hg
```

### Examples

Analytical propagation of a hermite-gaussian with adjusted waist:
```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = hg(rs,rs,zs,m=1,n=2,w0=.5)
show_animation(ψs)
```

Superposition of hermite-gaussians. Note that we can propagate backwards using negative values of `z`.
```julia
rs = LinRange(-4,4,256) 
zs = LinRange(-.5,.7,128)
ψ₀ = hg(rs,rs) + hg(rs,rs,n=2)/√2
ψs = free_propagation(ψ₀,rs,rs,zs)
show_animation(ψs)
```

## Diagonal Hermite-Gaussian

```@docs
diagonal_hg
```

### Examples

Analytical propagation of a hermite-gaussian with adjusted wavenumber:
```julia
rs = LinRange(-5,5,256) 
zs = LinRange(0,1,32)
ψs = diagonal_hg(rs,rs,zs,m=4,n=3,k=2)
show_animation(ψs)
```

## Lenses

We also implement lenses:

```@docs
lens
tilted_lens
```