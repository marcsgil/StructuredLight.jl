# Examples

## Mode Converter

In this example, we reproduce the results of [Beijersbergen, Marco W., et al. "Astigmatic laser mode converters and transfer of orbital angular momentum." Optics Communications 96.1-3 (1993): 123-132.](https://www.sciencedirect.com/science/article/pii/003040189390535D), where it is shown that a tilted lens can "transform" a Laguerre-Gauss mode in a diagonal Hermite-Gauss mode.

```@example
using StructuredLight, CairoMakie
k = 1
f = √2
d = f / √2
w = √((2 + √2) * f / k)
rs = LinRange(-4w, 4w, 512)
zs = LinRange(0, 2d, 64)

ψ₀ = lg(rs, rs, -d; w, k, l=1, p=0) .* lens(rs, rs, Inf, f; k)
ψ₁ = free_propagation(ψ₀, rs, rs, 2d; k)
ψs = free_propagation(ψ₀, rs, rs, zs; k)

visualize(abs2.([ψ₀, ψ₁] |> stack))
```

## Astigmatic Conversion

In this example, we reproduce the results of [Pravin Vaity et al., "Measuring topological charge of optical vortices using a tilted convex lens," Phys. Lett. A, vol. 377, no. 15, pp. 1154-1156, 2013. DOI: 10.1016/j.physleta.2013.02.030](https://www.sciencedirect.com/science/article/abs/pii/S0375960113001953?casa_token=4qY1zlrA1jAAAAAA:siRwxg9tPju8XHJkGtAjGVXJacg7pBbaZyFJUQscNaQplQ2ciYyoMQOlTexOlyaW9VSQBDViPph4), where it is shown that a tilted lens can "transform" a Laguerre-Gauss mode in a diagonal Hermite-Gauss mode.

```@example
using StructuredLight, CairoMakie

#All quantities have unit of (inverse) meter

w = 0.16e-3 #Waist
λ = 632.8e-9 #Wavelength
k = 2π / λ #Wavenumber
f = 50e-2 #Focal length of the lens

z₀ = 3.1 #Distance away from the focus where the beam encounters the lens

z_cr = z₀ / (z₀ / f - 1) #Conversion distance

ξ = deg2rad(6) #Tilting angle

# Now, we set up our grid and the initial profile by including the action of a tilted lens:
rs = LinRange(-70w, 70w, 1024)
ψ₀ = lg(rs, rs, z₀, l=3; w, k) .* tilted_lens(rs, rs, f, ξ; k) #Applies the lens

# Finally, we propagate. 
# Note that we introduce scalings, because, otherwise, the beam would be to small.
zs = z_cr .* LinRange(0.97, 1.03, 64)
scalings = 0.015 .* vcat(LinRange(2.4, 1, 32), LinRange(1, 2.4, 32))
ψ = free_propagation(ψ₀, rs, rs, zs, k=k, scalings)
anim = save_animation(abs2.(ψ), "tilted_lens.mp4", framerate=12)
nothing # hide
```

![](tilted_lens.mp4)

By changing the initial angular momentum, one obtains different HG modes.

## Diffraction of vortices through a triangular aperture

See B. Pinheiro da Silva, G. H. dos Santos, A. G. de Oliveira, N. Rubiano da Silva, W. T. Buono, R. M. Gomes, W. C. Soares, A. J. Jesus-Silva, E. J. S. Fonseca, P. H. Souto Ribeiro, and A. Z. Khoury, "Observation of a triangular-lattice pattern in nonlinear wave mixing with optical vortices," Optica 9, 908-912 (2022)

```@example
using StructuredLight, CairoMakie
xs = LinRange(-16, 16, 512)
ys = LinRange(-16, 16, 512)
zs = LinRange(0, 0.6, 64)

# Create initial beam with triangular aperture
initial_beam = lg(xs, ys, l=3)
triangular_mask = triangle(xs, ys, 3.0)
shaped_initial = initial_beam .* triangular_mask

# Propagate the shaped beam
propagated = free_propagation(shaped_initial, xs, ys, zs)
save_animation(abs2.(propagated), "shaped_beam_prop.mp4")
nothing # hide
```

![](shaped_beam_prop.mp4)

## Diffraction Rings in Kerr Media

Here, we reproduce the results of [E. V. Garcia Ramirez, M. L. Arroyo Carrasco, M. M. Mendez Otero, S. Chavez Cerda, and M. D. Iturbe Castillo, "Far field intensity distributions due to spatial self phase modulation of a Gaussian beam by a thin nonlocal nonlinear media," Opt. Express 18, 22067-22079 (2010)](https://opg.optica.org/oe/fulltext.cfm?uri=oe-18-21-22067&id=206115), where it is studied the appearance of rings in the far field profile of a Gaussian beam after it crosses a nonlinear medium.

```@example
#First we import the package.
using StructuredLight, CairoMakie

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
visualize(abs2.(get_images(rs,ms,2,-4zᵣ,15,10)))
```

The other figures are just a variation of this one, by changing the distance from the waist and `n`.

Alternatively, we can solve the complete nonlinear Schrödinger equation to get the initial profile:

```@example
using StructuredLight, CairoMakie
rs = LinRange(-10,10,512)
zᵣ = 1/2
ψ₀ = lg(rs,rs,-4zᵣ)
M = maximum(abs2,ψ₀)

n = 2
ψ₁ = kerr_propagation(ψ₀,rs,rs,.01*zᵣ,512,g=400*n*π/M)

#The output of this free propagation should then be equal to the previous case with `m=2`:
ψ = free_propagation(ψ₁,rs,rs,10,7)
visualize(abs2.(ψ))
```