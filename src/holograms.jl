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

const xs_besselj1 = LinRange(0, 0.5818, 512)
const ys_besselj1 = inverse.(x -> besselj(1, x), xs_besselj1, 0, 1.84)

inverse_besselj1(x) = zero_order_interpolation(x, xs_besselj1, ys_besselj1)

"""function standardize(holo, max_modulation)
    m, M = extrema(holo)
    @. round(UInt8, max_modulation * (holo - m) / (M - m))
end

function generate_hologram(desired, incoming, x, y, max_modulation, x_period, y_period)
    relative = desired ./ incoming
    M = maximum(abs, relative)

    kx = 2π / (x[begin+1] - x[begin]) / x_period
    ky = 2π / (y[begin+1] - y[begin]) / y_period

    ψ = @. inverse_besselj1(abs(relative) * 0.5818 / M)

    @tullio ψ[j, k] *= sin((kx * x[j] + ky * y[k]) + angle(relative[j, k]))

    standardize(ψ, round(max_modulation * 0.586))
end"""