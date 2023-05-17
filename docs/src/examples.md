# Examples

## Astigmatic Conversion

In this example, we reproduce the results of [this work](https://www.sciencedirect.com/science/article/abs/pii/S0375960113001953?casa_token=4qY1zlrA1jAAAAAA:siRwxg9tPju8XHJkGtAjGVXJacg7pBbaZyFJUQscNaQplQ2ciYyoMQOlTexOlyaW9VSQBDViPph4), where it is shown that a tilted lens can "transform" a Laguerre-Gauss mode in a diagonal Hermite-Gauss mode.

```julia
#Here, we initialize the package and define the experimental parameters:
using StructuredLight
using CUDA

#All quantities have unit of (inverse) meter

w0 = 0.16e-3 #Waist
λ = 632.8e-9 #Wavelength
k = 2π/λ #Wavenumber
f = 50e-2 #Focal length of the lens

z₀ = 3.1 #Distance away from the focus where the beam encounters the lens

z_cr = z₀/(z₀/f-1) #Conversion distance

ξ = deg2rad(6) #Tilting angle

# Now, we set up our grid and the initial profile by including the action of a tilted lens:
rs = LinRange(-70w0,70w0,1024)
ψ₀ = lg(rs,rs,z₀,l=3,w0=w0,k=k) .* tilted_lens(rs,rs,f,ξ,k=k) |> cu

# Finally, we propagate. 
# Note that we introduce scalings, because, otherwise, the beam would be to small.
zs = z_cr .* LinRange(.97,1.03,64)
scalings = 0.015 .* vcat(LinRange(2.4,1,32),LinRange(1,2.4,32))
ψ = free_propagation(ψ₀,rs,rs,zs,k=k,scaling=scalings)
anim = show_animation(ψ,ratio=1/4)
```

By changing the initial angular momentum, one obtains different HG modes.

## Diffraction Rings in Kerr Media

Here, we reproduce the results of [this work](https://opg.optica.org/oe/fulltext.cfm?uri=oe-18-21-22067&id=206115), where it is studied the appearance of rings in the far field profile of a gaussian beam after it crosses a nonlinear medium.

```julia
#First we import the package.
using StructuredLight

#Define a function that calculates the nonlinear phase term.
function non_linear_phase(ψ,m,n)
    M = maximum(abs2,ψ)
    [ cis( π * n * (abs2(ψ[j,k])/M)^(m/2) ) for j ∈ axes(ψ,1), k ∈ axes(ψ,2) ]
end

#Define a function that calculates the three images.
function get_images(rs,ms,n,z₀,z,scalling)
    ψ = Array{ComplexF64}(undef,length(rs),length(rs),length(ms))

    ψ₀ = lg(rs,rs,z₀)

    for (i,m) in enumerate(ms)
        ψ₁ = ψ₀ .* non_linear_phase(ψ₀,m,n)
        ψ[:,:,i] = free_propagation(ψ₁,rs,rs,z,scalling)
    end

    ψ
end

#The values of `m` that are used throughout the article:
ms = (1,2,4)

#The Rayleigh range:
zᵣ = 1/2;

# This would be the figure 1:
rs = LinRange(-10,10,512)
visualize(get_images(rs,ms,2,-4zᵣ,15,10))
```

The other figures are just a variation of this one, by changing the distance from the waist and `n`.

Alternatively, we can solve the complete nonlinear Schrödinger equation to get the initial profile:

```julia
using StructuredLight
rs = LinRange(-10,10,512)
ψ₀ = lg(rs,rs,-4zᵣ)
M = maximum(abs2,ψ₀)

n = 2
ψ₁ = kerr_propagation(ψ₀,rs,rs,.01*zᵣ,512,g=400*n*π/M)

#The output of this free propagation should then be equal to the previous case with `m=2`:
ψ = free_propagation(ψ₁,rs,rs,10,7)
visualize(ψ)
```