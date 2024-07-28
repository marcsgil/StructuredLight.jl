function distribute(N, v)
    Δ = v[end] - v[1]
    @assert Δ > 0 "The total interval is null"

    δs = [v[i] - v[i-1] for i in 2:length(v)]
    for δ in δs
        @assert δ ≥ 0 "The vector `v` must be in crescent order"
    end

    [round(Int, N * δ / Δ, RoundUp) for δ in δs]
end

@kernel function mul_kernel!(ψ, ϕ)
    i, j = @index(Global, NTuple)
    ψ[i, j] *= ϕ[i, j]
end

@kernel function self_phase_kernel!(ψ, factor)
    i, j = @index(Global, NTuple)
    ψ[i, j] *= cis(factor * abs2(ψ[i, j]))
end

@kernel function phase_kernel!(phases, quadratic_qs, Δz)
    i, j = @index(Global, NTuple)
    phases[i, j] = cis(Δz * quadratic_qs[i, j])
end

function dispersion_step!(ψ, phases, plan, iplan, mul_func!)
    plan * ψ
    mul_func!(ψ, phases; ndrange=size(ψ))
    iplan * ψ
end

function type_2A_step!(ψ, phases, plan, iplan, factor, repetitions, mul_func!, self_phase_func!, ndrange)
    if !iszero(repetitions)
        self_phase_func!(ψ, factor / 2; ndrange)
        for _ in 1:repetitions-1
            dispersion_step!(ψ, phases, plan, iplan, mul_func!)
            self_phase_func!(ψ, factor; ndrange)
        end
        dispersion_step!(ψ, phases, plan, iplan, mul_func!)
        self_phase_func!(ψ, factor / 2; ndrange)
    end
end

function kerr_propagation_loop!(dest, ψ₀, phases, quadratic_qs, z, divisions, g, k, plan, iplan,
    mul_func!, self_phase_func!, phase_func!)
    Δz = z / divisions
    phase_evolution_factor = g * Δz / 2k

    ndrange = size(ψ₀)
    phase_func!(phases, quadratic_qs, Δz; ndrange)

    type_2A_step!(ψ₀, phases, plan, iplan, phase_evolution_factor, divisions,
        mul_func!, self_phase_func!, ndrange)

    fftshift!(dest, ψ₀)
end

"""
    kerr_propagation(ψ₀,xs,ys,zs,total_steps;k=1,g=1)

Solve `∇² ψ + 2ik ∂_z ψ = - g |ψ|² ψ` under the initial condition `ψ₀`.

`xs` and `ys` are the grids over which `ψ₀` is calculated.

The outputs are saved at every `zs`, which is a number or a collection of numbers.

`total_steps` is the number of steps over which we discretize the propagation. The larger the `total_steps`, the better the precision and the slower is the calculation.
"""
function kerr_propagation(ψ₀, xs, ys, zs, total_steps; k=1, g=1)
    Zs = vcat(0, zs)

    steps = distribute(total_steps, Zs)

    result = similar(ψ₀, size(ψ₀)..., length(zs))

    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    ψ = ifftshift(ψ₀)

    qxs = to_device(reciprocal_grid(xs, shift=true), ψ₀)
    qys = to_device(reciprocal_grid(ys, shift=true), ψ₀)

    quadratic_qs = @. -(qxs^2 + qys'^2) / 2k
    phases = similar(ψ₀)

    backend = get_backend(ψ)
    mul_func! = mul_kernel!(backend, 256)
    self_phase_func! = self_phase_kernel!(backend, 256)
    phase_func! = phase_kernel!(backend, 256)

    for (i, divisions) in enumerate(steps)
        kerr_propagation_loop!(view(result, :, :, i), ψ, phases, quadratic_qs, Zs[i+1] - Zs[i], divisions, g, k,
            plan, iplan, mul_func!, self_phase_func!, phase_func!)
    end

    result
end