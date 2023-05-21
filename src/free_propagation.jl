function _free_propagation!(ψ₀,qxs,qys,zs,k)
    fft!(ψ₀)

    reciprocal_phases = reciprocal_quadratic_phase(qxs,qys,k)
    ψ = apply_reciprocal_phases(ψ₀,reciprocal_phases,zs)

    ifft!(ψ,(1,2))
end

function _free_propagation!(ψ₀,xs,ys,zs,qxs,qys,scaling,k)
    direct_phases = direct_quadratic_phase(xs,ys,k)
    ψ = apply_direct_phases(ψ₀,direct_phases,zs,scaling)

    fft!(ψ,(1,2))

    reciprocal_phases = reciprocal_quadratic_phase(qxs,qys,k)
    apply_reciprocal_phases!(ψ,reciprocal_phases,zs,scaling)

    ifft!(ψ,(1,2))
    apply_direct_phases!(ψ,direct_phases,zs,scaling)
end

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

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view
    
    ψ = _free_propagation!(shifted_ψ₀,qxs,qys,zs,k)

    fftshift_view(ψ,(1,2))
end

function free_propagation(ψ₀,xs,ys,zs,scaling;k=1)
    @assert length(zs) == length(scaling) "`zs` and `scaling` should have the same length"

    FFTW.set_num_threads(8)

    shifted_ψ₀ = ifftshift_view(ψ₀)

    shifted_xs = xs |> ifftshift_view
    shifted_ys = ys |> ifftshift_view

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view
    
    ψ = _free_propagation!(shifted_ψ₀,shifted_xs,shifted_ys,zs,qxs,qys,scaling,k)

    fftshift_view(ψ,(1,2))
end

function free_propagation(ψ₀,xs,ys,z::Number;k=1)
    dropdims(free_propagation(ψ₀,xs,ys,[z];k=k),dims=3)
end

function free_propagation(ψ₀,xs,ys,z::Number,scaling;k=1)
    dropdims(free_propagation(ψ₀,xs,ys,[z],[scaling];k=k),dims=3)
end