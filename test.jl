using StructuredLight
using CUDA

rs = LinRange(-4,4,256)
zs = LinRange(0,1,32)
scaling = fill(2,length(zs))

ψ₀ = lg(rs,rs) |> cu;

ψ = free_propagation(ψ₀,rs,rs,zs,scaling=scaling)

visualize(ψ[:,:,2])
vi

x = rand(256,256)
y = similar(x)

@benchmark copy!($y,$x)

@benchmark StructuredLight.free_propagation_old($ψ₀,$rs,$rs,$zs)
@benchmark free_propagation($ψ₀,$rs,$rs,$zs)
@benchmark CUDA.@sync free_propagation($ψ₀,$rs,$rs,zs)