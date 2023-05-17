function dispersion_step!(ψ,kernel,plan,iplan)
    plan*ψ
    map!(*,ψ,ψ,kernel)
    iplan*ψ
end

propagation_kernel(q,k,z) = @. cis( -z/2k * sum(abs2,q) )


"""
    free_propagation(ψ₀,xs,ys,z;k=1,scaling=1)
    free_propagation(ψ₀,xs,ys,z::AbstractArray;k=1,scaling::AbstractArray=ones(length(z)))

Propagate an inital profile `ψ₀`.

The propagation is the solution of `∇² ψ + 2ik ∂_z ψ = 0` at distance `z` under the initial condition `ψ₀`.

`xs` and `ys` are the grids over which `ψ₀` is calculated.

If `z` is an `AbstractArray`, the output is a 3D array representing the solution at every element of `z`.

If scaling isn't provided, the input and output grids are the same. Otherwise, the output at a distance `z[n]` is calculated on a scalled grid defined by `scaling[n] * xs` and `scaling[n] * ys`.

`k` is the wavenumber.
"""
function free_propagation(ψ₀,xs,ys,z::Number,k=1,scaling=1)
    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view

    cache = ifftshift(ψ₀)

    fft!(cache)
    @tullio cache[i,j] *= cis( - z / 2k * ( qxs[i]^2 + qys[j]^2 ) )
    ifft!(cache)

    fftshift_view(cache)
end


function free_propagation(ψ₀,xs,ys,zs::AbstractArray;k=1,scaling::AbstractArray=ones(length(zs)))
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    result = similar(ψ₀,size(ψ₀)...,length(zs))

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view

    Threads.@threads for n in eachindex(zs)
        cache = fftshift(ψ₀)
            
        plan*cache
        z = zs[n]
        @tullio cache[i,j] *= cis( - z / 2k * ( qxs[i]^2 + qys[j]^2 ) )
        iplan*cache

        fftshift!(view(result,:,:,n),cache)
    end

    result
end