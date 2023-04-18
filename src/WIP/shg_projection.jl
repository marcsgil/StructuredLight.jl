include("../fft_utils.jl")
include("../initial_profiles.jl")
include("../projection.jl")
include("shg_simulation.jl")
##
using Plots,LaTeXStrings,PlotThemes
Plots.theme(:dark)
default(fontfamily="Computer Modern",
	    linewidth=2, framestyle=:box, grid=true, minorticks=5,dpi=400)
scalefontsizes(1.3)
##
#Physical Parameters
w₀ = 0.16e-3
λ = 632.8e-9
γ₀ = w₀/√2
k = 2π/λ
zr = k*γ₀^2
grid = DFTGrid(512,5γ₀,2)
##
#SHG simulation parameters
k₁ = k
k₂ = k
k₃ = 2k
χ = 1e-10
χ₁ = χ
χ₂ = χ
χ₃ = χ
Δr = interval(grid)
α₁ = im/(2*k₁*Δr^2)
α₂ = im/(2*k₂*Δr^2)
α₃ = im/(2*k₃*Δr^2)
β₁ = im*χ₁*k₁/2
β₂ = im*χ₂*k₂/2
β₃ = im*χ₃*k₃/2
parameters = (α₁,α₂,α₃,β₁,β₂,β₃)
##
#Initial Conditions
l₁ = 2
l₂ = 2
l₃ = l₁ + l₂
ψ₁ = CuArray(evaluate_at_direct_grid(r->LG(r...,γ₀=γ₀,k=k₁,l=l₁),grid))
ψ₂ = CuArray(evaluate_at_direct_grid(r->LG(r...,γ₀=γ₀,k=k₂,l=l₂),grid))
ψ₃ = zero(ψ₁)
##
zs = LinRange(0,.5zr,100)
sol = simulate_shg(ψ₁,ψ₂,ψ₃,zs,parameters);
ψ3 = [Array(sol[j][:,:,3]) for j in eachindex(sol)]
sim = DiffEqArray(ψ3,zs)
##
heatmap(abs2.(ψ3[end]),color=:hot,aspect_ratio=:equal,
size=(450,400),xlims=(1,grid.N),ylims=(1,grid.N),xticks=false,yticks=false)
##
projections = project_result(0:5,l₃,sim,grid,γ₀=γ₀/√2,k=k₃)
cs = abs.(projections)
plot(projections.t/zr,cs,label=reshape(map(p->L"|c_{%$p}|", 0:5),1,6),xlabel=L"z/z_r",legend=:topleft,title=L"l_1=%$l₁;l_2=%$l₂")
#savefig("SHG/Projeções/RegimeForte/l1=$l₁;l2=$l₂.png")
##
lin_reg(x,y) = (x'x)\(x'y)
lin_reg(results.t,cs[1])/lin_reg(results.t,cs[2])