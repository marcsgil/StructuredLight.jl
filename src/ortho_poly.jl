"""
    laguerre(x, n::Integer, α::Real=0)

Evaluate the `n` th Laguerre polynomial at `x`.
"""
function laguerre(x, n::Integer, α::Real=0)
    T = typeof(x)
    α = convert(T, α)
    p0, p1 = one(T), -x + (α + 1)
    n == 0 && return p0
    for k = 1:n-1
        p1, p0 = ((2k + α + 1) / (k + 1) - x / (k + 1)) * p1 - (k + α) / (k + 1) * p0, p1
    end
    p1
end

"""
    hermite(x, n::Integer)

Evaluate the `n` th Hermite polynomial at `x`.
"""
function hermite(x, n::Integer)
    T = typeof(x)
    p0, p1 = one(T), 2x
    n == 0 && return p0
    for k = 1:n-1
        p1, p0 = 2x * p1 - 2k * p0, p1
    end
    p1
end

"""
    jacobi(x, n, a, b)

    Evaluate the `n` th Jacobi polynomial at `x` with parameters `a` and `b`.
"""
function jacobi(x, n, a, b)
    ox = one(x)
    zx = zero(x)
    if n == 0
        return ox
    elseif n == 1
        return ox / 2 * (a - b + (a + b + 2) * x)
    end

    p0 = ox
    p1 = ox / 2 * (a - b + (a + b + 2) * x)
    p2 = zx

    for i = 1:(n-1)
        a1 = 2 * (i + 1) * (i + a + b + 1) * (2 * i + a + b)
        a2 = (2 * i + a + b + 1) * (a * a - b * b)
        a3 = (2 * i + a + b) * (2 * i + a + b + 1) * (2 * i + a + b + 2)
        a4 = 2 * (i + a) * (i + b) * (2 * i + a + b + 2)
        p2 = ox / a1 * ((a2 + a3 * x) * p1 - a4 * p0)

        p0 = p1
        p1 = p2
    end

    return p2
end

"""
    radial_poly(ρ, m, n)

    Evaluate the radial polynomial of the `n` th Zernike polynomial at `ρ` with azimuthal order `m`.
"""
function radial_poly(ρ, m, n)
    δ = n - m
    !iseven(δ) && return zero(ρ)

    order = δ ÷ 2
    (-1)^order * ρ^m * jacobi(1 - 2ρ^2, order, m, 0)
end

"""
    zernike_polynomial(x, y, m, n)

    Evaluate the `n` th Zernike polynomial at `(x, y)` with azimuthal order `m`.
"""
function zernike_polynomial(x, y, m, n)
    ρ = sqrt(x^2 + y^2)
    θ = atan(y, x)

    trig_func = m ≥ 0 ? cos : sin

    radial_poly(ρ, abs(m), n) * trig_func(m * θ)
end