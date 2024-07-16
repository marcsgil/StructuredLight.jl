function inverse(f, y, bracket)
    find_zero(x -> f(x) - y, bracket)
end

function interpolated_inverse(f, domain, N=512)
    xs = LinRange(f(domain[1]), f(domain[2]), N)
    ys = [inverse(f, x, domain) for x ∈ xs]

    linear_interpolation(xs, ys)
end

inverse_besselj1 = interpolated_inverse(besselj1, (0, 1.84))

function standardize(holo, max_modulation)
    m, M = extrema(holo)
    @. round(UInt8, max_modulation * (holo - m) / (M - m))
end

function ψ(a, ϕ)
    inverse_besselj1(a * 0.5818) * sin(ϕ)
end

function sin_ϕ(x, y, kx, ky, relative)
    s, c = sincos(kx * x + ky * y)
    s * real(relative) + c * imag(relative)
end


function generate_hologram(desired, incoming, x, y, max_modulation, x_period, y_period)
    relative = desired ./ incoming
    M = maximum(abs, relative)

    kx = 2π / (x[begin+1] - x[begin]) / x_period
    ky = 2π / (y[begin+1] - y[begin]) / y_period

    ψ = @. inverse_besselj1(abs(relative) * 0.5818 / M)

    @tullio ψ[j, k] *= sin((kx * x[k] + ky * y[j]) + angle(relative[j, k]))

    standardize(ψ, round(max_modulation * 0.586))
end