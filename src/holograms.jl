function find_zero(f, x₁, x₂; tol=1e-6, maxiter=10^3)
    @assert f(x₁) * f(x₂) ≤ 0 "Error: f(x₁) and f(x₂) must have opposite signs"

    isapprox(f(x₁), 0; atol=tol) && return 2x₁ / 2
    isapprox(f(x₂), 0; atol=tol) && return 2x₂ / 2

    for _ ∈ 1:maxiter
        x = (x₁ + x₂) / 2
        y = f(x)

        if abs(y) < tol
            return x
        end

        if f(x) * f(x₁) < 0
            x₂ = x
        else
            x₁ = x
        end
    end
end

function inverse(f, y, x₁, x₂)
    find_zero(x -> f(x) - y, x₁, x₂)
end

function closest_index(number, grid)
    step = grid[2] - grid[1]
    idx = round(Int, (number - grid[1]) / step) + 1
    clamp(idx, 1, length(grid))
end

function zero_order_interpolation(x, xs, ys)
    i = closest_index(x, xs)
    ys[i]
end

abs_div(x, y) = abs(x / y)

function base_max_abs_div(x, y, M)
    for (x, y) ∈ zip(x, y)
        M = max(M, abs_div(x, y))
    end

    M
end

function max_abs_div(x::AbstractArray{T1}, y::AbstractArray{T2}, base_size,
    M=floatmin(float(promote_type(real(T1), real(T2))))) where {T1,T2}

    if length(x) < base_size || length(y) < base_size
        return base_max_abs_div(x, y, M)
    end

    x₁ = view(x, 1:length(x)÷2)
    x₂ = view(x, length(x)÷2+1:length(x))
    y₁ = view(y, 1:length(y)÷2)
    y₂ = view(y, length(y)÷2+1:length(y))

    M₁ = Threads.@spawn max_abs_div(x₁, y₁, base_size, M)
    M₂ = Threads.@spawn max_abs_div(x₂, y₂, base_size, M)

    M = max(fetch(M₁), fetch(M₂))
end

normalize_angle(angle) = mod2pi(angle + π) - π

abstract type HologramMethod end
struct BesselJ1 <: HologramMethod end
struct Simple <: HologramMethod end

const x_min_besselj1 = 0
const x_max_besselj1 = 0.5818
const y_min_besselj1 = 0
const y_max_besselj1 = 1.82337
const xs_besselj1 = LinRange(x_min_besselj1, x_max_besselj1, 1024)
const ys_besselj1 = inverse.(x -> besselj1(x), xs_besselj1, y_min_besselj1, y_max_besselj1)

inverse_besselj1(x) = zero_order_interpolation(x, xs_besselj1, ys_besselj1)

@kernel function hologram_kernel!(dest, desired, incoming, x, y, two_pi_modulation, M, kx, ky, ::Type{BesselJ1})
    i, j = @index(Global, NTuple)

    relative = desired[i, j] / incoming[i, j]
    ψ = (inverse_besselj1(x_max_besselj1 * abs(relative) / M)
         *
         sin((kx * x[i] + ky * y[j]) + angle(relative)))

    dest[i, j] = round(two_pi_modulation / 1.17 * (ψ / y_max_besselj1 + 1) / 2)
end

@kernel function hologram_kernel!(dest, desired, incoming, x, y, two_pi_modulation, M, kx, ky, ::Type{Simple})
    i, j = @index(Global, NTuple)

    relative = desired[i, j] / incoming[i, j]
    ψ = (inverse_besselj1(x_max_besselj1 * abs(relative) / M)
         *
         sin((kx * x[i] + ky * y[j]) + angle(relative)))

    ψ = normalize_angle()

    dest[i, j] = round(two_pi_modulation / 1.17 * (ψ / y_max_besselj1 + 1) / 2)
end

function generate_hologram!(dest, desired, incoming, x, y, max_modulation, x_period, y_period, method::Type{T}=BesselJ1) where {T<:HologramMethod}
    M = max_abs_div(desired, incoming, length(desired) ÷ (2 * Threads.nthreads()))

    kx = 2π / interval(x) / x_period
    ky = 2π / interval(y) / y_period

    backend = get_backend(dest)
    kernel! = hologram_kernel!(backend)
    kernel!(dest, desired, incoming, x, y, max_modulation, M, kx, ky, method; ndrange=size(dest))
end

function generate_hologram(desired, incoming, x, y, max_modulation, x_period, y_period, method::Type{T}=BesselJ1) where {T<:HologramMethod}
    dest = similar(desired)
    generate_hologram!(dest, desired, incoming, x, y, max_modulation, x_period, y_period, method)
    dest
end