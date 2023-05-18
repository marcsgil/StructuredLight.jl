function get_propagation_kernel(qxs,qys,zs,k)
    @tullio mini_kernel[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    @tullio kernel[i,j,l] := cis( mini_kernel[i,j] * zs[l] )
end

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
function free_propagation(ψ₀,xs,ys,zs::AbstractArray;k=1,scaling::AbstractArray=ones(length(zs)))
    FFTW.set_num_threads(8)
    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view

    kernel = get_propagation_kernel(qxs,qys,zs,k)

    shifted_ψ₀ = ifftshift(ψ₀)
    fft!(shifted_ψ₀)

    @tullio kernel[i,j,l] *= shifted_ψ₀[i,j]

    ifft!(kernel,(1,2))

    fftshift_view(kernel,(1,2))
end

function free_propagation(ψ₀,xs,ys,z::Number,k=1,scaling=1)
    dropdims(free_propagation(ψ₀,xs,ys,[z];k=k,scaling=[scaling]),dims=3)
end