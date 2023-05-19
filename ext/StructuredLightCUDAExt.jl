module StructuredLightCUDAExt

using StructuredLight, CUDA, CUDA.CUFFT, Tullio, FourierTools
using StructuredLight: reciprocal_grid

function StructuredLight.free_propagation(ψ₀::CuArray,xs,ys,zs::AbstractArray;k=1)
    shifted_ψ₀ = ifftshift(ψ₀)
    fft!(shifted_ψ₀)

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view

    @tullio cache[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    @tullio ψ₁[i,j,l] := cis( cache[i,j] * zs[l] )
    ψ₁ = CuArray(ψ₁)
    ψ₁ .*= shifted_ψ₀

    ifft!(ψ₁,(1,2))

    fftshift(ψ₁,(1,2))
end

function StructuredLight.free_propagation(ψ₀::CuArray,xs,ys,zs::AbstractArray,scaling::AbstractArray;k=1)
    shifted_ψ₀ = ifftshift(ψ₀)

    shifted_xs = xs |> ifftshift_view
    shifted_ys = ys |> ifftshift_view

    @tullio cache1[i,j] := k * ( shifted_xs[i]^2 + shifted_ys[j]^2 ) / 2
    @tullio cache2[i,j,l] := cis( cache1[i,j] * ( 1 - scaling[l] ) / zs[l] ) / scaling[l]
    ψ₁ = CuArray(cache2) .* shifted_ψ₀

    fft!(ψ₁,(1,2))

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view
    @tullio cache3[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    @tullio cache2[i,j,l] = cis( cache3[i,j] * zs[l] / scaling[l] )
    ψ₁ .*= CuArray(cache2)

    ifft!(ψ₁,(1,2))
    @tullio cache2[i,j,l] = cis( - cache1[i,j] * ( 1 - scaling[l] ) * scaling[l] / zs[l])
    ψ₁ .*= CuArray(cache2)

    fftshift(ψ₁,(1,2))
end
    
end