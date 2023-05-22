function StructuredLight._free_propagation!(ψ₀::CuArray,qxs::CuArray,qys::CuArray,zs::CuArray,k)
    fft!(ψ₀)

    @tullio phases[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    typeof(phases)
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( phases[i,j] * zs[l] )

    ifft!(ψ,(1,2))
end

function StructuredLight.free_propagation(ψ₀::CuArray,xs,ys,zs;k=1)
    shifted_ψ₀ = ifftshift(ψ₀)

    qxs = StructuredLight.reciprocal_grid(xs,shift=true) |> CuArray
    qys = StructuredLight.reciprocal_grid(ys,shift=true) |> CuArray
    
    ψ = StructuredLight._free_propagation!(shifted_ψ₀,qxs,qys,CuArray(zs),k)

    zs isa Number ? dropdims( fftshift(ψ,(1,2)),dims=3) : fftshift(ψ,(1,2))
end

function StructuredLight._free_propagation!(ψ₀::CuArray,xs::CuArray,ys::CuArray,zs::CuArray,qxs::CuArray,qys::CuArray,scaling::CuArray,k)
    @tullio direct_phases[i,j] := k * ( xs[i]^2 + ys[j]^2 ) / 2
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( direct_phases[i,j] * ( 1 - scaling[l] ) / zs[l] ) / scaling[l]

    fft!(ψ,(1,2))

    @tullio reciprocal_phases[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    @tullio ψ[i,j,l] *= cis( reciprocal_phases[i,j] * zs[l] / scaling[l] )

    ifft!(ψ,(1,2))
    @tullio ψ[i,j,l] *= cis( - direct_phases[i,j] * ( 1 - scaling[l] ) * scaling[l] / zs[l])
end

function StructuredLight.free_propagation(ψ₀::CuArray,xs,ys,zs,scaling;k=1)
    @assert length(zs) == length(scaling) "`zs` and `scaling` should have the same length"
    CUDA.@allowscalar @assert 0 ∉ zs "This method does not support `zs` containing `0`"

    shifted_ψ₀ = ifftshift(ψ₀)

    direct_xgrid = fftshift(xs) |> CuArray
    direct_ygrid = fftshift(ys) |> CuArray

    qxs = StructuredLight.reciprocal_grid(xs,shift=true) |> CuArray
    qys = StructuredLight.reciprocal_grid(ys,shift=true) |> CuArray
    
    ψ = StructuredLight._free_propagation!(shifted_ψ₀,direct_xgrid,direct_ygrid,CuArray(zs),qxs,qys,CuArray(scaling),k)

    zs isa Number ? dropdims( fftshift(ψ,(1,2)),dims=3) : fftshift(ψ,(1,2))
end