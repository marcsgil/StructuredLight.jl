function type_2A_step!(ψ,phases,plan,iplan,factor,repetitions)
    for _ in 1:repetitions
        map!(ψ->evolve_phase(ψ,abs2(ψ)*factor/2),ψ,ψ)
        dispersion_step!(ψ,phases,plan,iplan)
        map!(ψ->evolve_phase(ψ,abs2(ψ)*factor/2),ψ,ψ)
    end
end

function partition(a::Integer,b::Integer)
    @assert a ≥ b

    result = fill(a ÷ b, b)

    for n in 1:mod(a,b)
        result[n] +=1
    end

    result
end

function kerr_propagation(ψ₀,xs,ys,zs,total_steps=0;k=1,g=1)
    #solves ∇² ψ + 2ik ∂_z ψ = - g |ψ|² ψ with initial condition ψ₀

    @assert iszero(first(zs))

    if iszero(total_steps) || total_steps < length(zs) - 1
        total_steps = length(zs)
    end

    steps = partition(total_steps,length(zs)-1)

    results = similar(ψ₀,size(ψ₀)...,length(zs))
    results[:,:,1] = ψ₀

    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    cache = ifftshift(ψ₀)

    ks = reciprocal_grid(xs,ys) |> collect |> ifftshift

    for (i,n) in enumerate(steps)
        Δz = (zs[i+1] - zs[i])/n
        phase_evolution_factor = g*Δz/(2k)

        dispersion_phases = get_dispersion_phases(Δz,k,ks,ψ₀)

        type_2A_step!(cache,dispersion_phases,plan,iplan,phase_evolution_factor,n)

        fftshift!(view(results,:,:,i+1), cache)
    end

    results
end