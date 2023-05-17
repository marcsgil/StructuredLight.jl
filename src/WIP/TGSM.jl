TGSM(x,y,x′,y′;σs=1,σg=1,u=1,k=1) = exp(-(x^2+y^2+x′^2+y′^2)/(4σs^2)-( (x-x′)^2+(y-y′)^2 )/(2σg^2) -k*u*im*(x*y′-x′*y))
TGSM(rs;σs=1,σg=1,u=1,k=1) = map(R->TGSM(R...,σs=σs,σg=σg,u=u,k=k),Iterators.product(rs,rs,rs,rs))


#Correlation function propagation
function get_intensity(W)
    @tullio I[x,y] := real(W[x,y,x,y])
end

function pure_state_density_operator(state)
    @tullio W[i,j,k,l] := state[i,j]*state[k,l]
end

function free_propagation(W::AbstractArray{T,4},grid,zs;k=1) where T<: Number
    cW = cu(W)

    plan1 = plan_ifft!(cW,(1,2))
    plan2 = plan_fft!(cW,(3,4))

    iplan1 = plan_fft!(cW,(1,2))
    iplan2 = plan_ifft!(cW,(3,4))

    ks = reciprocal_grid(grid)

    exponents = map(ks->ks[3]^2+ks[4]^2-ks[1]^2-ks[2]^2,Iterators.product(ks,ks,ks,ks)) |> ifftshift |> cu

    result = Array{Float64}(undef,size(W,1),size(W,2),length(zs))
    cache = similar(cW)
    
    for n in eachindex(zs)
        ifftshift!(cache,cW)
        plan2*(plan1*cache)
        factor = zs[n]/(2k)
        map!((x,k2)->(x*cis(factor*k2)),cache,cache,exponents)
        iplan2*(iplan1*cache)
        fftshift!(view(result,:,:,n),cache |> Array |> get_intensity)
    end

    result
end

using ParaxialBeamPropagation
##
k = 2π/(632e-6)
σs = .5
γ₀ = σs*√2
zr = γ₀^2
β = 1
η = 1
#σg=σs/β
#u=η*β^2/σs^2
σg=.1
u=.001
grid = DFTGrid(100,4γ₀)
rs = direct_grid(grid)
zs = LinRange(0,zr,10)

W = TGSM(rs,σs=σs,σg=σs,u=u,k=k)

ψs = Array{Float64}(undef,length(rs),length(rs),length(zs),2)
ψs[:,:,:,1] = abs2.(HG(rs,rs,zs,HGConfig(γ₀=γ₀)))
ψs[:,:,:,2] = free_propagation(W,grid,zs)

vizualize_beam(ψs,rs,rs,zs/zr,csmax=[1,1])
#get_animation(Array(ψs),rs,rs,zs,"test.mp4",20,csmax=[1,1])
##
ξ(a,b,u) = (2a+2b-√(4a^2+8a*b+(k*u)^2))/(2a+2b+√(4a^2+8a*b+(k*u)^2))
Γ₀(a,b,u) = (4a^2+8a*b+(k*u)^2)^(1/4)
t(b,u) = √((b+k*u/2)/(b-k*u/2))
λ(n,m,a,b,u) = (1-ξ(a,b,u))*ξ(a,b,u)^(abs(m)/2+n)*t(b,u)^m
Λ(n,m,a,b,u) = (t(b,u)^m)*(ξ(a,b,u)^(abs(m)/2+n))

a = 1/(4σs^2)
b = 1/(2σg^2)

ξ(a,b,u)

λ(0,-1,a,b,u)/λ(0,0,a,b,u)

Λ(1,,a,b,u)