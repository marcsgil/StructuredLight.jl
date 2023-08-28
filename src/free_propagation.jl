"""
    free_propagation(ψ₀, xs, ys, z::Number [, scaling]; k=1)
    free_propagation(ψ₀,xs,ys,z::AbstractArray [, scaling]; k=1)

Propagate an inital profile `ψ₀`.

The propagation is the solution of `∇² ψ + 2ik ∂_z ψ = 0` at distance `z` under the initial condition `ψ₀`.

`xs` and `ys` are the grids over which `ψ₀` is calculated.

If `z` is an `AbstractArray`, the output is a 3D array representing the solution at every element of `z`.

The output at a distance `z[n]` is calculated on a scalled grid defined by `scaling[n] * xs` and `scaling[n] * ys`.

`k` is the wavenumber.
"""
function free_propagation(ψ₀,xs,ys,zs;k=1)
    FFTW.set_num_threads(8)

    shifted_ψ₀ = ifftshift(ψ₀)

    qxs = reciprocal_grid(xs,shift=true)
    qys = reciprocal_grid(ys,shift=true)
    
    ψ = _free_propagation!(shifted_ψ₀,qxs,qys,zs,k)

    zs isa Number ? dropdims( fftshift_view(ψ,(1,2)),dims=3) : fftshift_view(ψ,(1,2))
end

function _free_propagation!(ψ₀,qxs,qys,zs,k)
    fft!(ψ₀)

    @tullio phases[i,j] := - ( qxs[j]^2 + qys[i]^2 ) / 2k
    typeof(phases)
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( phases[i,j] * zs[l] )

    ifft!(ψ,(1,2))
end

function free_propagation(ψ₀,xs,ys,zs,scaling;k=1)
    @assert length(zs) == length(scaling) "`zs` and `scaling` should have the same length"
    @assert 0 ∉ zs "This method does not support `zs` containing `0`"

    FFTW.set_num_threads(8)

    shifted_ψ₀ = ifftshift_view(ψ₀)

    direct_xgrid = fftshift_view(xs)
    direct_ygrid = fftshift_view(ys)

    qxs = reciprocal_grid(xs,shift=true)
    qys = reciprocal_grid(ys,shift=true)
    
    ψ = _free_propagation!(shifted_ψ₀,direct_xgrid,direct_ygrid,zs,qxs,qys,scaling,k)

    zs isa Number ? dropdims( fftshift_view(ψ,(1,2)),dims=3) : fftshift_view(ψ,(1,2))
end

function _free_propagation!(ψ₀,xs,ys,zs,qxs,qys,scaling,k)
    @tullio direct_phases[i,j] := k * ( xs[j]^2 + ys[i]^2 ) / 2
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( direct_phases[i,j] * ( 1 - scaling[l] ) / zs[l] ) / scaling[l]

    fft!(ψ,(1,2))

    @tullio reciprocal_phases[i,j] := - ( qxs[j]^2 + qys[i]^2 ) / 2k
    @tullio ψ[i,j,l] *= cis( reciprocal_phases[i,j] * zs[l] / scaling[l] )

    ifft!(ψ,(1,2))
    @tullio ψ[i,j,l] *= cis( - direct_phases[i,j] * ( 1 - scaling[l] ) * scaling[l] / zs[l])
end

function non_spectral_propagation(ψ₀,xs,ys,z;k=1)
    Δx = xs[2] - xs[1]
    Δy = ys[2] - ys[1]
    α = √(k/(2im*z))
    kernel = similar(ψ₀)
    @tullio kernel[i,j] = π / 4 * (erf(α*(j+1)*Δx) - erf(α*j*Δx) ) * (erf(α*(i+1)*Δy) - erf(α*i*Δy) )
    (ifft(fft((ifftshift(ψ₀))) .* fft((ifftshift(kernel)))))
end