function dispersion_step!(ψ,kernel,plan,iplan)
    plan*ψ
    map!(*,ψ,ψ,kernel)
    iplan*ψ
end

evolve_phase(ψ,phase) = cis(phase) * ψ

propagation_kernel(κ,k,z) = @. cis( -z/2k * sum(abs2,κ) )

function scaling_phase(r,k,z,scaling)
    if iszero(z)
        1
    else
        @. cis( k * ( 1 - scaling ) * sum(abs2,r) / (2z))
    end   
end

function free_propagation(ψ₀,xs,ys,z::Number,k,plan,iplan,scaling)
    if isone(scaling)
        cache = ifftshift(ψ₀)
    else
        scaling_phases = oftype(ψ₀, scaling_phase(direct_grid(xs,ys),k,z,scaling))
        cache = ifftshift( ψ₀ .* scaling_phases ./ scaling )
    end  

    kernel = oftype(ψ₀, propagation_kernel(reciprocal_grid(xs,ys),k,z/scaling) ) |> ifftshift

    dispersion_step!(cache,kernel,plan,iplan)
    fftshift!(kernel,cache)

    if isone(scaling)
        kernel
    else
        kernel .= kernel.* conj(scaling_phases).^scaling
    end
end

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
function free_propagation(ψ₀,xs,ys,z;k=1,scaling=1)
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)
    free_propagation(ψ₀,xs,ys,z,k,plan,iplan,scaling)
end

function free_propagation(ψ₀,xs,ys,z::AbstractArray;k=1,scaling::AbstractArray=ones(length(z)))
    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    result = similar(ψ₀,size(ψ₀)...,length(z))

    κ = reciprocal_grid(xs,ys) |> collect |> ifftshift
    R = direct_grid(xs,ys) |> collect

    for (n,z) in enumerate(z)
        if isone(scaling[n])
            cache = ifftshift(ψ₀)
        else
            scaling_phases = oftype(ψ₀, scaling_phase(R,k,z,scaling[n]))
            cache = ifftshift( ψ₀ .* scaling_phases ./ scaling[n] )
        end  
        
        kernel = oftype(ψ₀, propagation_kernel(κ,k,z/scaling[n]) )
        
        dispersion_step!(cache,kernel,plan,iplan)
        fftshift!(view(result,:,:,n),cache)
        
        if !isone(scaling[n])
            result[:,:,n] = view(result,:,:,n) .* conj(scaling_phases).^scaling[n]
        end
    end

    result
end