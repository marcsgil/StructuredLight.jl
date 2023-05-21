module StructuredLightCUDAExt

using StructuredLight, CUDA, CUDA.CUFFT, FourierTools
using Tullio, CUDAKernels, KernelAbstractions
using StructuredLight: reciprocal_grid, _free_propagation!, distribute, reciprocal_quadratic_phase, kerr_propagation_loop!

include("propagation_functions.jl")

function StructuredLight.free_propagation(ψ₀::CuArray,xs,ys,zs::AbstractArray;k=1)
    shifted_ψ₀ = ifftshift(ψ₀)
    qxs = reciprocal_grid(xs) |> ifftshift_view |> CuArray
    qys = reciprocal_grid(ys) |> ifftshift_view |> CuArray
    zs_cuda = CuArray(zs)

    ψ = _free_propagation!(shifted_ψ₀,qxs,qys,zs_cuda,k)

    fftshift(ψ,(1,2))
end

function StructuredLight.free_propagation(ψ₀::CuArray,xs,ys,zs::AbstractArray,scaling::AbstractArray;k=1)
    @assert length(zs) == length(scaling) "`zs` and `scaling` should have the same length"

    shifted_ψ₀ = ifftshift(ψ₀)

    shifted_xs = xs |> ifftshift_view |> CuArray
    shifted_ys = ys |> ifftshift_view |> CuArray
    zs_cuda = CuArray(zs)

    qxs = reciprocal_grid(xs) |> ifftshift_view |> CuArray
    qys = reciprocal_grid(ys) |> ifftshift_view |> CuArray

    scaling_cuda = CuArray(scaling)
    
    ψ₁ = _free_propagation!(shifted_ψ₀,shifted_xs,shifted_ys,zs_cuda,qxs,qys,scaling_cuda,k)

    fftshift(ψ₁,(1,2))
end

function StructuredLight.kerr_propagation(ψ₀::CuArray,xs,ys,zs::AbstractArray,total_steps;k=1,g=1)
    Zs = vcat(0,zs)

    steps = distribute(total_steps,Zs)

    result = similar(ψ₀,size(ψ₀)...,length(zs))

    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    ψ = ifftshift(ψ₀)

    qxs = reciprocal_grid(xs) |> ifftshift_view |> CuArray
    qys = reciprocal_grid(ys) |> ifftshift_view |> CuArray

    phases = reciprocal_quadratic_phase(qxs,qys,k)
    kernel = similar(ψ₀)

    for (i,divisions) in enumerate(steps)
        kerr_propagation_loop!(view(result,:,:,i),ψ,kernel,phases,Zs[i+1] - Zs[i],divisions,g,k,plan,iplan)
    end

    result
end
    
end