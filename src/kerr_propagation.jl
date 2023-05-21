function distribute(N,v)
    Δ = v[end] - v[1]
    @assert Δ > 0 "The total interval is null"

    δs = [ v[i] - v[i-1] for i in 2:length(v) ]   
    for δ in δs
        @assert δ ≥ 0 "The vector `v` must be in crescent order"
    end

    [ round(Int, N * δ / Δ, RoundUp) for δ in δs ]
end

function type_2A_step!(ψ,phases,plan,iplan,factor,repetitions)
    if !iszero(repetitions)
        apply_kerr_phase!(ψ,factor/2)
        for _ in 1:repetitions-1
            dispersion_step!(ψ,phases,plan,iplan)
            apply_kerr_phase!(ψ,factor)
        end
        dispersion_step!(ψ,phases,plan,iplan)
        apply_kerr_phase!(ψ,factor/2)
    end
end

function kerr_propagation_loop!(dest,ψ₀,kernel,phases,z,divisions,g,k,plan,iplan)
    Δz = z / divisions
    phase_evolution_factor = g * Δz / 2k

    fourier_propagation_kernel!(kernel,phases,Δz)

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
    Zs = vcat(0,zs)

    steps = distribute(total_steps,Zs)

    result = similar(ψ₀,size(ψ₀)...,length(zs))

    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    ψ = ifftshift(ψ₀)

    qxs = reciprocal_grid(xs) |> ifftshift_view
    qys = reciprocal_grid(ys) |> ifftshift_view

    phases = reciprocal_quadratic_phase(qxs,qys,k)
    kernel = similar(ψ₀)

    for (i,divisions) in enumerate(steps)
        kerr_propagation_loop!(view(result,:,:,i),ψ,kernel,phases,Zs[i+1] - Zs[i],divisions,g,k,plan,iplan)
    end

    result
end