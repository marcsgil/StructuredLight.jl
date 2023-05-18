module StructuredLightCUDAExt

using StructuredLight, CUDA, CUDA.CUFFT

function StructuredLight.free_propagation(ψ₀::CUDA.CuArray,xs,ys,zs::AbstractArray;k=1,scaling::AbstractArray=ones(length(zs)))
    qxs = StructuredLight.reciprocal_grid(xs) |> StructuredLight.ifftshift_view
    qys = StructuredLight.reciprocal_grid(ys) |> StructuredLight.ifftshift_view
    
    kernel = StructuredLight.get_propagation_kernel(qxs,qys,zs,k) |> CuArray

    shifted_ψ₀ = ifftshift(ψ₀)
    fft!(shifted_ψ₀)

    kernel .*= shifted_ψ₀

    ifft!(kernel,(1,2))

    fftshift(kernel,(1,2))
end
    
end