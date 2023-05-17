using StructuredLight
using CUDA

rs = LinRange(-4,4,256)
zs = LinRange(0,1,32)
scaling = fill(2,length(zs))

ψ₀ = lg(rs,rs);

ψ = free_propagation(ψ₀,rs,rs,zs,scaling=scaling);

visualize(ψ[:,:,1])
vi

@benchmark StructuredLight.free_propagation_old($ψ₀,$rs,$rs,$zs)
@benchmark free_propagation($ψ₀,$rs,$rs,$zs)
@benchmark CUDA.@sync free_propagation($ψ₀,$rs,$rs,zs)
##
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
ψ₀ = lg(rs,rs,z₀,l=3,w0=w0,k=k) .* tilted_lens(rs,rs,f,ξ,k=k) |> cu

# Finally, we propagate. 
# Note that we introduce scalings, because, otherwise, the beam would be to small.
zs = z_cr .* LinRange(.97,1.03,64)
scalings = 0.015 .* vcat(LinRange(2.4,1,32),LinRange(1,2.4,32))
ψ = free_propagation(ψ₀,rs,rs,zs,k=k,scaling=scalings)
show_animation(ψ,ratio=1/2)