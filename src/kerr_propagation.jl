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

function distribute(N,v)
    δs = [ v[i] - v[i-1] for i in 2:length(v) ]

    for δ in δs
        @assert δ ≥ 0 "The points are not in crescent order"
    end

    Δ = v[end] - v[1]
    @assert Δ > 0 "The total interval is null"

    result = [ round(Int, N * δs[n] / Δ, RoundUp) for n in 1:length(δs) - 1 ]
    vcat(result, N - sum(result))
end

function kerr_propagation(ψ₀,xs,ys,zs,total_steps;k=1,g=1)
    #solves ∇² ψ + 2ik ∂_z ψ = - g |ψ|² ψ with initial condition ψ₀

    Zs = vcat(0,zs)

    steps = distribute(total_steps,Zs)

    results = similar(ψ₀,size(ψ₀)...,length(zs))

    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    cache = ifftshift(ψ₀)

    ks = reciprocal_grid(xs,ys) |> collect |> ifftshift

    for (i,n) in enumerate(steps)
        Δz = (Zs[i+1] - Zs[i])/n
        phase_evolution_factor = g*Δz/(2k)

        kernel = convert(typeof(ψ₀), fourier_propagation_kernel.(ks,k,Δz))

        type_2A_step!(cache,kernel,plan,iplan,phase_evolution_factor,n)

        fftshift!(view(results,:,:,i), cache)
    end

    if zs isa Number
        dropdims(results,dims=3)
    else
        results
    end
end
