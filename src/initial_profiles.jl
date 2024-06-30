function assert_hermite_indices(idxs...)
    for idx ∈ idxs
        @assert idx ≥ 0 "All indices must be non negative"
    end
end

float_type(args...) = float(typeof(sum(first, args)))

function normalization_hg(m, n, γ::T1=1) where {T1}
    #T = float(T1)
    #convert(T, prod(inv, 1:m) * prod(inv, 1:n) / γ / √(π * 2^(m + n)))

    oftype(float(γ), inv(γ * √(π * 2^(m + n) * factorial(n) * factorial(m))))
end

function _hg(x, y, s, c; m::Integer=0, n::Integer=0, γ=one(x))
    X = (x * c + y * s) / γ
    Y = (-x * s + y * c) / γ

    hermite(X, m) * hermite(Y, n) * exp((-X^2 - Y^2) / 2)
end

"""
    hg(x::Real, y::Real; m::Integer=0, n::Integer=0, w=one(x))

    hg(x, y; m::Integer=0, n::Integer=0, w=one(eltype(x)))

    hg(x::Real, y::Real, z::Real;
    m::Integer=0, n::Integer=0, w=one(x), k=one(x))

    hg(x, y, z::Real;
    m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

    hg(x, y, z;
    m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute the Hermite-Gaussian mode.

The optional keyword arguments are:

`m`: horizontal index

`n`: vertical index

`w0`: beam's waist

`k`: wavenumber

The different signatures evaluate points in a grid, e.g.

`hg(xs,ys,zs) = [hg(x,y,z) for x ∈ xs, y ∈ ys, z ∈ zs]`

but the left hand side is faster.
"""
function rotated_hg(x::Real, y::Real; θ, m::Integer=0, n::Integer=0, w=one(x))
    assert_hermite_indices(m, n)
    T = float_type(x, y, w)
    γ = convert(T, w / √2)
    s, c = sincos(θ)

    normalization_hg(m, n, γ) * _hg(x, y, s, c; m, n, γ)
end

function rotated_hg(x, y; θ, m::Integer=0, n::Integer=0, w=one(eltype(x)))
    assert_hermite_indices(m, n)
    T = float_type(x, y, w)
    γ = convert(T, w / √2)
    s, c = sincos(θ)

    N = normalization_hg(m, n, γ)

    @tullio _[j, k] := N * _hg(x[j], y[k], s, c; m, n, γ)
end

function _hg(x::Real, y::Real, α::Complex, s, c; m::Integer=0, n::Integer=0, γ=one(x))
    X = (x * c + y * s) / γ
    Y = (-x * s + y * c) / γ

    α * exp(α * (-X^2 - Y^2) / 2) * hermite(abs(α) * X, m) * hermite(abs(α) * Y, n)
end

function rotated_hg(x::Real, y::Real, z::Real;
    θ, m::Integer=0, n::Integer=0, w=one(x), k=one(x))
    assert_hermite_indices(m, n)

    T = float(typeof(sum((x, y, z, w, k))))
    γ = convert(T, w / √2)
    k = convert(T, k)
    s, c = sincos(θ)

    α = inv(1 + im * z / (k * γ^2))
    prefactor = normalization_hg(m, n, γ) * cis((m + n) * angle(α))
    prefactor * _hg(x, y, α, s, c; m, n, γ)
end

function rotated_hg(x, y, z::Real;
    θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))
    assert_hermite_indices(m, n)

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)
    s, c = sincos(θ)

    α = inv(1 + im * z / (k * γ^2))
    prefactor = normalization_hg(m, n, γ) * cis((m + n) * angle(α))

    @tullio _[j, k] := prefactor * _hg(x[j], y[k], α, s, c; m, n, γ)
end

function rotated_hg(x, y, z;
    θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))
    assert_hermite_indices(m, n)

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)
    s, c = sincos(θ)

    function prefactor(z)
        α = inv(1 + im * z / (k * γ^2))
        cis((m + n) * angle(α))
    end
    N = normalization_hg(m, n, γ)

    @tullio _[j, k, l] := N * prefactor(z[l]) * _hg(x[j], y[k], α, s, c; m, n, γ)
end

hg(x, y; kwargs...) = rotated_hg(x, y; θ=zero(first(x)), kwargs...)
hg(x, y, z; kwargs...) = rotated_hg(x, y, z; θ=zero(first(x)), kwargs...)

diagonal_hg(x, y; kwargs...) = rotated_hg(x, y; θ=oftype(float(first(x)), π / 4), kwargs...)
diagonal_hg(x, y, z; kwargs...) = rotated_hg(x, y, z; θ=oftype(float(first(x)), π / 4), kwargs...)

function core_lg(x, y, α, γ₀, p::Integer, l::Integer)
    r2 = (x^2 + y^2) / γ₀^2
    α * exp(-α * r2 / 2) * (abs(α) * (x + im * sign(l) * y) / γ₀)^abs(l) * laguerre(p, abs(l), abs2(α) * r2)
end

"""
    normalization_lg(;p,l,γ₀=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(; p, l, γ₀=1)
    convert(float(eltype(γ₀)), √inv(prod(p+1:p+abs(l)) * π) / γ₀)
end

