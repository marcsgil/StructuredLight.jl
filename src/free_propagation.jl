"""
    free_propagation(ψ₀,xs,ys,z;k=1,scaling=1)
    free_propagation(ψ₀,xs,ys,z::AbstractArray;k=1,scaling::AbstractArray=ones(length(z)))

Propagate an inital profile `ψ₀`.

The propagation is the solution of `∇² ψ + 2ik ∂_z ψ = 0` at distance `z` under the initial condition `ψ₀`.

`xs` and `ys` are the grids over which `ψ₀` is calculated.

If `z` is an `AbstractArray`, the output is a 3D array representing the solution at every element of `z`.

The output at a distance `z[n]` is calculated on a scalled grid defined by `scaling[n] * xs` and `scaling[n] * ys`.

`k` is the wavenumber.
"""
function free_propagation(ψ₀,xs,ys,zs::AbstractArray;k=1)
    FFTW.set_num_threads(8)

    shifted_ψ₀ = ifftshift(ψ₀)
    fft!(shifted_ψ₀)

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view

    @tullio cache[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    @tullio ψ₁[i,j,l] := shifted_ψ₀[i,j] * cis( cache[i,j] * zs[l] )

    ifft!(ψ₁,(1,2))

    fftshift_view(ψ₁,(1,2))
end

function free_propagation(ψ₀,xs,ys,z::Number,k=1)
    dropdims(free_propagation(ψ₀,xs,ys,[z];k=k),dims=3)
end

function free_propagation(ψ₀,xs,ys,zs::AbstractArray,scaling::AbstractArray;k=1)
    FFTW.set_num_threads(8)

    shifted_ψ₀ = ifftshift_view(ψ₀)

    shifted_xs = xs |> ifftshift_view
    shifted_ys = ys |> ifftshift_view

    @tullio cache1[i,j] := k * ( shifted_xs[i]^2 + shifted_ys[j]^2 ) / 2
    @tullio ψ₁[i,j,l] := shifted_ψ₀[i,j] * cis( cache1[i,j] * ( 1 - scaling[l] ) / zs[l] ) / scaling[l]

    fft!(ψ₁,(1,2))

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view
    @tullio cache2[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    @tullio ψ₁[i,j,l] *= cis( cache2[i,j] * zs[l] / scaling[l] )

    ifft!(ψ₁,(1,2))
    @tullio ψ₁[i,j,l] *= cis( - cache1[i,j] * ( 1 - scaling[l] ) * scaling[l] / zs[l])

    fftshift_view(ψ₁,(1,2))
end