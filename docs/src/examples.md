# Examples

## Mode Converter


In this example, we reproduce the results of [Beijersbergen, Marco W., et al. "Astigmatic laser mode converters and transfer of orbital angular momentum." Optics Communications 96.1-3 (1993): 123-132.](https://www.sciencedirect.com/science/article/pii/003040189390535D), where it is shown that a tilted lens can "transform" a Laguerre-Gauss mode in a diagonal Hermite-Gauss mode.

```@example
using StructuredLight
k = 1
f = √2
d = f / √2
w0 = √((2 + √2) * f / k)
rs = LinRange(-4w0, 4w0, 512)
zs = LinRange(0, 2d, 64)
##
ψ₀ = lg(rs, rs, -d; w0, k, l=1, p=0) .* lens(rs, rs, Inf, f; k)
ψ₁ = free_propagation(ψ₀, rs, rs, 2d; k)
ψs = free_propagation(ψ₀, rs, rs, zs; k)

visualize(ψ₀)
visualize(ψ₁)
```

## Astigmatic Conversion

In this example, we reproduce the results of [Pravin Vaity et al., "Measuring topological charge of optical vortices using a tilted convex lens," Phys. Lett. A, vol. 377, no. 15, pp. 1154-1156, 2013. DOI: 10.1016/j.physleta.2013.02.030](https://www.sciencedirect.com/science/article/abs/pii/S0375960113001953?casa_token=4qY1zlrA1jAAAAAA:siRwxg9tPju8XHJkGtAjGVXJacg7pBbaZyFJUQscNaQplQ2ciYyoMQOlTexOlyaW9VSQBDViPph4), where it is shown that a tilted lens can "transform" a Laguerre-Gauss mode in a diagonal Hermite-Gauss mode.

```@example
#Here, we initialize the package and define the experimental parameters:
using StructuredLight

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
ψ₀ = lg(rs,rs,z₀,l=3,w0=w0,k=k) .* tilted_lens(rs,rs,f,ξ,k=k)

# Finally, we propagate. 
# Note that we introduce scalings, because, otherwise, the beam would be to small.
zs = z_cr .* LinRange(.97,1.03,64)
scalings = 0.015 .* vcat(LinRange(2.4,1,32),LinRange(1,2.4,32))
ψ = free_propagation(ψ₀,rs,rs,zs,k=k,scalings)
anim = show_animation(ψ,ratio=1/2,fps=12)
```

By changing the initial angular momentum, one obtains different HG modes.

## Diffraction Rings in Kerr Media

Here, we reproduce the results of [E. V. Garcia Ramirez, M. L. Arroyo Carrasco, M. M. Mendez Otero, S. Chavez Cerda, and M. D. Iturbe Castillo, "Far field intensity distributions due to spatial self phase modulation of a Gaussian beam by a thin nonlocal nonlinear media," Opt. Express 18, 22067-22079 (2010)](https://opg.optica.org/oe/fulltext.cfm?uri=oe-18-21-22067&id=206115), where it is studied the appearance of rings in the far field profile of a Gaussian beam after it crosses a nonlinear medium.

```@example
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
visualize(get_images(rs,ms,2,-4zᵣ,15,10),ratio=2)
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