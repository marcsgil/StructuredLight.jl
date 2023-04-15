function fourier_propagation_kernel(κ,k,z)
    factor = -z/(2k)
    cis( factor * sum(abs2,κ) )
end

function inner_scalling_phase(r,k,z,scalling)
    if iszero(z)
        one(z)
    else
        cis( k * ( 1 - scalling ) * sum(abs2,r) / (2z))
    end
end

function outter_scalling_phase(r,k,z,scalling)
    if iszero(z)
        one(z)
    else
        cis( - k * ( 1 - scalling ) * scalling * sum(abs2,r) / (2z) )
    end
end

map_convert(ψ,f,args...) = convert(typeof(ψ),map(f,args...))

function free_propagation(ψ₀,xs,ys,z::Number,k,plan,iplan)
    kernel = convert(typeof(ψ₀), fourier_propagation_kernel.(reciprocal_grid(xs,ys),k,z)) |> ifftshift
    cache = ifftshift(ψ₀)

    fftshift(dispersion_step!(cache,kernel,plan,iplan))
end

function free_propagation(ψ₀,xs,ys,z::Number,k,plan,iplan,scalling)
    cache = ifftshift(ψ₀ .* map_convert(ψ₀, r -> inner_scalling_phase(r,k,z,scalling) / scalling, direct_grid(xs,ys)))

    kernel = map_convert(ψ₀, κ-> fourier_propagation_kernel(κ,k,z/scalling), reciprocal_grid(xs,ys)) |> ifftshift
    

    fftshift(dispersion_step!(cache,kernel,plan,iplan)) .* map_convert(ψ₀, r -> outter_scalling_phase(r,k,z,scalling), direct_grid(xs,ys))
end

"""
    free_propagation(ψ₀,xs,ys,z;k=1)

Propagate an inital profile `ψ₀` over a distance `z`. 

`xs` and `ys` are the grids over which `ψ₀` is calculated.

`k` is the wavenumber.
"""
function free_propagation(ψ₀,xs,ys,z;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)
    free_propagation(ψ₀,xs,ys,z,k,plan,iplan)
end

function free_propagation(ψ₀,xs,ys,z,scalling;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)
    free_propagation(ψ₀,xs,ys,z,k,plan,iplan,scalling)
end

"""
    free_propagation(ψ₀,xs,ys,zs::AbstractArray;k=1)

Propagate an inital profile `ψ₀` over every distance in the array `zs`. 

`xs` and `ys` are the grids over which `ψ₀` is calculated.

`k` is the wavenumber.
"""
function free_propagation(ψ₀,xs,ys,zs::AbstractArray;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    result = similar(ψ₀,size(ψ₀)...,length(zs))

    ks = reciprocal_grid(xs,ys) |> collect |> ifftshift

    for (n,z) in enumerate(zs)
        cache = ifftshift(ψ₀)
        phases = convert(typeof(ψ₀), fourier_propagation_kernel.(ks,k,z))
        dispersion_step!(cache,phases,plan,iplan)
        fftshift!(view(result,:,:,n),cache)
    end

    result
end

function free_propagation(ψ₀,xs,ys,zs::AbstractArray,scallings;k=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    result = similar(ψ₀,size(ψ₀)...,length(zs))

    ks = reciprocal_grid(xs,ys) |> collect |> ifftshift

    for (n,z) in enumerate(zs)
        cache = ifftshift(ψ₀ .* map_convert(ψ₀, r -> inner_scalling_phase(r,k,z,scallings[n]) / scallings[n], direct_grid(xs,ys)))
        kernel = map_convert(ψ₀, κ-> fourier_propagation_kernel(κ,k,z/scallings[n]), reciprocal_grid(xs,ys)) |> ifftshift
        dispersion_step!(cache,kernel,plan,iplan)
        fftshift!(view(result,:,:,n),cache)
        result[:,:,n] = view(result,:,:,n) .* map_convert(ψ₀, r -> outter_scalling_phase(r,k,z,scallings[n]), direct_grid(xs,ys))
    end

    result
end