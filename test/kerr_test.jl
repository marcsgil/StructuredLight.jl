function Λ(η, p, q, l)
    L = abs(l)
    if q ≥ p
        (-1)^p * (2√η / (1 + η))^(L + 1) * √(prod(q+1:q+L) / prod(p+1:p+L)) *
        sum((-1)^n * binomial(q, n) * binomial(p + q + L - n, q + L)
            * ((1 - η) / (1 + η))^(p + q - 2n) for n in 0:p)
    else
        Λ(1 / η, q, p, l)
    end
end

function C(p, l)
    L = abs(l)
    2 * (-1)^p * binomial(2L, L - p) * √binomial(L + p, p) / (π * 3^((3L + 1) / 2))
end

function Φ(z, p)
    iszero(p) ? atan(z) : (1 - ((1 - im * z) / (1 + im * z))^p) / (2im * p)
end

function nonlinear_c(z, p, l)
    sum(C(q, l) * Φ(z, r) * Λ(3, r, p, l) * Λ(3, r, q, l) for q in 0:abs(l), r in 0:10^3) * im / 4
end

function calculate_projections(rs, zs, g, l)
    ψ₀ = lg(rs, rs; l)

    free_ψ = free_propagation(ψ₀, rs, rs, zs)

    δψ = (kerr_propagation(ψ₀, rs, rs, zs, 64, g=g) .- free_ψ) ./ g

    numerical_c = Matrix{real(eltype(δψ))}(undef, (length(zs), abs(l) + 3))

    analytic_c = [abs(nonlinear_c(z, p, l)) for z in zs, p in 0:abs(l)+2]

    for p in axes(numerical_c, 2)
        corrected_ψ = lg(rs, rs, zs, p=p - 1, l=l, γ=1 / √3)
        numerical_c[:, p] = overlap(corrected_ψ, δψ, rs, rs) .|> abs
    end

    analytic_c, numerical_c
end

rs = LinRange(-22, 22, 256)
zs = LinRange(0, 5, 16)
g = 0.01
l = 2
analytic_c, numerical_c = calculate_projections(rs, zs, g, l)

isapprox(analytic_c, numerical_c, rtol=1e-3)