function fourier_propagation_kernel(qs,k,z)
    cis( - z * sum(abs2,qs) / (2k) )
end

function dispersion_step!(ψ,kernel,plan,iplan)
    plan*ψ
    map!(*,ψ,ψ,kernel)
    iplan*ψ
end

evolve_phase(ψ,phase) = cis(phase) * ψ

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

    Δ = v[end] - v[1]
    @assert Δ > 0 "The total interval is null"

    result = [ round(Int, N * δs[n] / Δ, RoundUp) for n in 1:length(δs) - 1 ]
    vcat(result, N - sum(result))
end

"""
    kerr_propagation(ψ₀,xs,ys,zs,total_steps;k=1,g=1)

Solve `∇² ψ + 2ik ∂_z ψ = - g |ψ|² ψ` under the initial condition `ψ₀`.

`xs` and `ys` are the grids over which `ψ₀` is calculated.

The outputs are saved at every `zs`, which is a number or a collection of numbers.

`total_steps` is the number of steps over which we discretize the propagation. The larger the `total_steps`, the better the precision and the slower is the calculation.
"""
function kerr_propagation(ψ₀,xs,ys,zs,total_steps;k=1,g=1)
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