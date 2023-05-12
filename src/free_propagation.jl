function fourier_propagation_kernel(κ,k,z)
    factor = -z/(2k)
    cis( factor * sum(abs2,κ) )
end

function inner_scaling_phase(r,k,z,scaling)
    if iszero(z) || isone(scaling)
        one(z)
    else
        cis( k * ( 1 - scaling ) * sum(abs2,r) / (2z))
    end
end

function outter_scaling_phase(r,k,z,scaling)
    if iszero(z) || isone(scaling)
        one(z)
    else
        cis( - k * ( 1 - scaling ) * scaling * sum(abs2,r) / (2z) )
    end
end

map_convert(ψ,f,args...) = convert(typeof(ψ),map(f,args...))

function free_propagation(ψ₀,xs,ys,z::Number,k,plan,iplan)
    kernel = convert(typeof(ψ₀), fourier_propagation_kernel.(reciprocal_grid(xs,ys),k,z)) |> ifftshift
    cache = ifftshift(ψ₀)

    fftshift(dispersion_step!(cache,kernel,plan,iplan))
end

function free_propagation(ψ₀,xs,ys,z::Number,k,plan,iplan,scaling)
    cache = ifftshift(ψ₀ .* map_convert(ψ₀, r -> inner_scaling_phase(r,k,z,scaling) / scaling, direct_grid(xs,ys)))

    kernel = map_convert(ψ₀, κ-> fourier_propagation_kernel(κ,k,z/scaling), reciprocal_grid(xs,ys)) |> ifftshift
    
    fftshift(dispersion_step!(cache,kernel,plan,iplan)) .* map_convert(ψ₀, r -> outter_scaling_phase(r,k,z,scaling), direct_grid(xs,ys))
end

"""
    free_propagation(ψ₀,xs,ys,z;k=1)
    free_propagation(ψ₀,xs,ys,z,scaling;k=1)
    free_propagation(ψ₀,xs,ys,z::AbstractArray;k=1)
    free_propagation(ψ₀,xs,ys,z::AbstractArray,scaling::AbstractArray;k=1)

Propagate an inital profile `ψ₀`.

The propagation is the solution of `∇² ψ + 2ik ∂_z ψ = 0` at distance `z` under the initial condition `ψ₀`.

`xs` and `ys` are the grids over which `ψ₀` is calculated.

If `z` is an `AbstractArray`, the output is a 3D array representing the solution at every element of `z`.

If scaling isn't provided, the input and output grids are the same. Otherwise, the output at a distance `z[n]` is calculated on a scalled grid defined by `scaling[n] * xs` and `scaling[n] * ys`.

`k` is the wavenumber.
"""
function free_propagation(ψ₀,xs,ys,z;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)
    free_propagation(ψ₀,xs,ys,z,k,plan,iplan)
end

function free_propagation(ψ₀,xs,ys,z,scaling;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)
    free_propagation(ψ₀,xs,ys,z,k,plan,iplan,scaling)
end

function free_propagation(ψ₀,xs,ys,z::AbstractArray;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    result = similar(ψ₀,size(ψ₀)...,length(z))

    ks = reciprocal_grid(xs,ys) |> collect |> ifftshift

    for (n,z) in enumerate(z)
        cache = ifftshift(ψ₀)
        phases = convert(typeof(ψ₀), fourier_propagation_kernel.(ks,k,z))
        dispersion_step!(cache,phases,plan,iplan)
        fftshift!(view(result,:,:,n),cache)
    end

    result
end

function free_propagation(ψ₀,xs,ys,z::AbstractArray,scaling::AbstractArray;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    result = similar(ψ₀,size(ψ₀)...,length(z))

    ks = reciprocal_grid(xs,ys) |> collect |> ifftshift

    for (n,z) in enumerate(z)
        cache = ifftshift(ψ₀ .* map_convert(ψ₀, r -> inner_scaling_phase(r,k,z,scaling[n]) / scaling[n], direct_grid(xs,ys)))
        kernel = map_convert(ψ₀, κ-> fourier_propagation_kernel(κ,k,z/scaling[n]), reciprocal_grid(xs,ys)) |> ifftshift
        dispersion_step!(cache,kernel,plan,iplan)
        fftshift!(view(result,:,:,n),cache)
        result[:,:,n] = view(result,:,:,n) .* map_convert(ψ₀, r -> outter_scaling_phase(r,k,z,scaling[n]), direct_grid(xs,ys))
    end

    result
end