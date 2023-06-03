function distribute(N,v)
    Δ = v[end] - v[1]
    @assert Δ > 0 "The total interval is null"

    δs = [ v[i] - v[i-1] for i in 2:length(v) ]   
    for δ in δs
        @assert δ ≥ 0 "The vector `v` must be in crescent order"
    end

    [ round(Int, N * δ / Δ, RoundUp) for δ in δs ]
end

function dispersion_step!(ψ,kernel,plan,iplan)
    plan*ψ
    @tullio ψ[i,j] *= kernel[i,j]
    iplan*ψ
end

function type_2A_step!(ψ,kernel,plan,iplan,factor,repetitions)
    if !iszero(repetitions)
        @tullio ψ[i,j] *= cis( factor * abs2(ψ[i,j]) / 2)
        for _ in 1:repetitions-1
            dispersion_step!(ψ,kernel,plan,iplan)
            @tullio ψ[i,j] *= cis( factor * abs2(ψ[i,j]) )
        end
        dispersion_step!(ψ,kernel,plan,iplan)
        @tullio ψ[i,j] *= cis( factor * abs2(ψ[i,j]) / 2)
    end
end

function kerr_propagation_loop!(dest,ψ₀,kernel,phases,z,divisions,g,k,plan,iplan)
    Δz = z / divisions
    phase_evolution_factor = g * Δz / 2k

    @tullio kernel[i,j] = cis( Δz * phases[i,j] )

    type_2A_step!(ψ₀,kernel,plan,iplan,phase_evolution_factor,divisions)

    fftshift!(dest, ψ₀)
end

"""
    kerr_propagation(ψ₀,xs,ys,zs,total_steps;k=1,g=1)

Solve `∇² ψ + 2ik ∂_z ψ = - g |ψ|² ψ` under the initial condition `ψ₀`.

`xs` and `ys` are the grids over which `ψ₀` is calculated.

The outputs are saved at every `zs`, which is a number or a collection of numbers.

`total_steps` is the number of steps over which we discretize the propagation. The larger the `total_steps`, the better the precision and the slower is the calculation.
"""
function kerr_propagation(ψ₀,xs,ys,zs,total_steps;k=1,g=1)
    Zs = vcat(0,Array(zs))

    steps = distribute(total_steps,Zs)

    result = similar(ψ₀,size(ψ₀)...,length(zs))

    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    ψ = ifftshift(ψ₀)

    qxs = reciprocal_grid(xs,shift=true)
    qys = reciprocal_grid(ys,shift=true)

    @tullio phases[i,j] := - ( qxs[j]^2 + qys[i]^2 ) / 2k
    kernel = similar(ψ₀)

    for (i,divisions) in enumerate(steps)
        kerr_propagation_loop!(view(result,:,:,i),ψ,kernel,phases,Zs[i+1] - Zs[i],divisions,g,k,plan,iplan)
    end

    result
end