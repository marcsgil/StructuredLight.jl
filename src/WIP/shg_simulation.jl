using CUDA,DifferentialEquations

using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@init_parallel_stencil(CUDA, ComplexF64, 2)

@parallel function step!(dψ₁,dψ₂,dψ₃,ψ₁,ψ₂,ψ₃,α₁,α₂,α₃,β₁,β₂,β₃)
    @inn(dψ₁) = @. α₁*(@d2_xi(ψ₁) + @d2_yi(ψ₁)) + β₁*conj(@inn(ψ₂))*@inn(ψ₃)
    @inn(dψ₂) = @. α₂*(@d2_xi(ψ₂) + @d2_yi(ψ₂)) + β₂*conj(@inn(ψ₁))*@inn(ψ₃)
    @inn(dψ₃) = @. α₃*(@d2_xi(ψ₃) + @d2_yi(ψ₃)) + β₃*@inn(ψ₁)*@inn(ψ₂)
    return
end

function step!(du,u,(α₁,α₂,α₃,β₁,β₂,β₃),z)
    @views @parallel step!(du[:,:,1],du[:,:,2],du[:,:,3],u[:,:,1],u[:,:,2],u[:,:,3],α₁,α₂,α₃,β₁,β₂,β₃)
end

function simulate_second_order_nl(ψ₁,ψ₂,ψ₃,zs,parameters)
    u = CuArray{ComplexF64}(undef, (size(ψ₁)...,3))
    u[:,:,1] = CuArray(ψ₁)
    u[:,:,2] = CuArray(ψ₂)
    u[:,:,3] = CuArray(ψ₃)
    prob = ODEProblem(step!, u, (0, last(zs)), parameters)
    solve(prob,SSPRK43(),saveat=zs)
end

function format(sol::ODESolution)
    formated_sol = Array(sol)
    if ndims(formated_sol) == 4
        permutedims(formated_sol,(1,2,4,3))
    else
        formated_sol
    end
end