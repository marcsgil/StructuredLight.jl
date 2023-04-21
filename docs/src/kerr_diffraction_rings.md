# Diffraction Rings in Kerr Media

Here, we reproduce the results of [this work](https://opg.optica.org/oe/fulltext.cfm?uri=oe-18-21-22067&id=206115), where it is studied the appearance of rings in the far field profile of a gaussian beam after it crosses a nonlinear medium.

First we import the package, define a function that calculates the nonlinear phase term, and a function that calculates the three images.

```julia
using StructuredLight

function non_linear_phase(ψ,m,n)
    M = maximum(abs2,ψ)
    [ cis( π * n * (abs2(ψ[j,k])/M)^(m/2) ) for j ∈ axes(ψ,1), k ∈ axes(ψ,2) ]
end

function get_images(rs,ms,n,z₀,z,scalling)
    ψ = Array{ComplexF64}(undef,length(rs),length(rs),length(ms))

    ψ₀ = lg(rs,rs,z₀)

    for (i,m) in enumerate(ms)
        ψ₁ = ψ₀ .* non_linear_phase(ψ₀,m,n)
        ψ[:,:,i] = free_propagation(ψ₁,rs,rs,z,scalling)
    end

    ψ
end
```

We also define the values of `m` that are used throughout the article, and also the Rayleigh range `zᵣ`.

```julia
ms = (1,2,4)

zᵣ = 1/2;
```

This would be the figure 1:
```julia
rs = LinRange(-10,10,512)
visualize(get_images(rs,ms,2,-4zᵣ,15,10))
```

The other figures are just a variation of this one, by changing the distance from the waist and `n`.

Alternatively, we can solve the complete nonlinear Schrödinger equation to get the initial profile:

```julia
rs = LinRange(-10,10,512)
ψ₀ = lg(rs,rs,-4zᵣ)
M = maximum(abs2,ψ₀)

n = 2
ψ₁ = kerr_propagation(ψ₀,rs,rs,.01*zᵣ,512,g=400*n*π/M);
```

The output of this free propagation should then be equal to the previous case with `m=2`:
```julia
ψ = free_propagation(ψ₁,rs,rs,10,7)
visualize(ψ)
```