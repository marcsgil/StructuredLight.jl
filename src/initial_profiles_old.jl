"""
    laguerre_coefficients(n::Integer,α=0)

Compute the coefficients of the nth generalized Laguerre Polynomial.
"""
function laguerre_coefficients(n, α=0)
    ntuple(i -> -(-1)^i * binomial(n + α, n - i + 1) / factorial(i - 1), n + 1)
end

"""
    normalization_lg(;p,l,γ₀=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(; p, l, γ₀=1)
    convert(float(eltype(γ₀)), √inv(prod(p+1:p+abs(l)) * π) / γ₀)
end

function laguerre(n::Int, α::Number, x)
    T = typeof(x)
    α = convert(T, α)
    p0, p1 = one(T), -x + (α + 1)
    n == 0 && return p0
    for k = 1:n-1
        p1, p0 = ((2k + α + 1) / (k + 1) - x / (k + 1)) * p1 - (k + α) / (k + 1) * p0, p1
    end
    p1
end

function core_lg(x, y, α, γ₀, p::Integer, l::Integer)
    r2 = (x^2 + y^2) / γ₀^2
    α * exp(-α * r2 / 2) * (abs(α) * (x + im * sign(l) * y) / γ₀)^abs(l) * laguerre(p, abs(l), abs2(α) * r2)
end

"""
    lg(x::Real,y::Real,z::Real=0;
        p::Integer=0,l::Integer=0,w0::Real=1,k::Real=1) where T<: Real

    lg(x::AbstractVector{T},y::AbstractVector{T},z::Real=0;
        p::Integer=0,l::Integer=0,w0::Real=1,k::Real=1) where T<: Real

    lg(x::AbstractVector{T},y::AbstractVector{T},z::AbstractVector{T};
        p::Integer=0,l::Integer=0,w0::Real=1,k::Real=1) where T<: Real

Compute the Laguerre-Gaussian mode. 

For the first signature, the mode is calculated at point `(x,y,z)`

For the second signature, the mode is calculated over a grid defined by `x` and `y` at a distance `z` from the focus.

For the third signature, the mode is calculated over a grid defined by `x`, `y` and `z`.

The optional keyword arguments are:

`p`: radial index

`l`: topological charge

`w0`: beam's waist

`k`: wavenumber
"""
function lg(x::AbstractVector{T}, y::AbstractVector{T}, z::Real=0;
    p::Integer=0, l::Integer=0, w0::Real=1, k::Real=1) where {T<:Real}

    @assert p ≥ 0

    γ₀ = convert(float(T), w0 / √2)
    k = convert(float(T), k)

    α = 1 / (1 + im * z / (k * γ₀^2))
    prefactor = normalization_lg(p=p, l=l, γ₀=γ₀) * cis((2p + abs(l)) * angle(α))

    @tullio result[j, i] := prefactor * core_lg(x[i], y[j], α, γ₀, p, l)
end

function lg(x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T};
    p::Integer=0, l::Integer=0, w0::Real=1, k::Real=1) where {T<:Real}

    @assert p ≥ 0

    γ₀ = convert(float(T), w0 / √2)
    k = convert(float(T), k)

    function f(x, y, α)
        normalization_lg(p=p, l=l, γ₀=γ₀) * cis((2p + abs(l)) * angle(α)) * core_lg(x, y, α, γ₀, p, l)
    end

    @tullio result[j, i, l] := f(x[i], y[j], inv(1 + im * z[l] / (k * γ₀^2)))
end

function lg(x::Real, y::Real;
    p::Integer=0, l::Integer=0, w0::Real=1)

    @assert p ≥ 0

    T = typeof(float(x + y + w0))
    γ₀ = convert(T, w0 / √2)
    N = normalization_lg(; p, l, γ₀)

    L = abs(l)
    σ = sign(l)
    r2 = (x^2 + y^2) / γ₀^2

    #coeffs = laguerre_coefficients(p, convert(T, L))

    #N * (x + im * σ * y)^L * evalpoly(r2, coeffs) * exp(-r2 / 2)

    N * (x + im * σ * y)^L * laguerre(p, L, r2) * exp(-r2 / 2)
end


function hermite_coefficients(n)
    if iseven(n)
        ntuple(l -> -factorial(n) * (-1)^(n ÷ 2 - l) / (factorial(2l - 2) * factorial(n ÷ 2 - l + 1)) |> Integer, n ÷ 2 + 1)
    else
        ntuple(l -> -factorial(n) * (-1)^(n ÷ 2 - l) / (factorial(2l - 1) * factorial(n ÷ 2 - l + 1)) |> Integer, n ÷ 2 + 1)
    end
end

"""
    hermite(x,n,coeffs)

Evaluate the `n` th Hermite polynomial at `x`, given the coefficients `coeffs`.
"""
hermite(x, n, coeffs) = iseven(n) ? evalpoly(4x^2, coeffs) : 2x * evalpoly(4x^2, coeffs)

"""
    normalization_hg(;m,n,γ₀=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_hg(m, n, γ₀)
    oftype(float(γ₀), inv(γ₀ * √(π * 2^(m + n) * factorial(n) * factorial(m))))
end

function normalization_hg(n, γ₀)
    oftype(float(γ₀), inv((pi^(1 // 4) * √(γ₀ * 2^(n) * factorial(n)))))
end

function hermite(n::Int, x)
    T = typeof(x)
    p0, p1 = one(T), 2x
    n == 0 && return p0
    for k = 1:n-1
        p1, p0 = 2x * p1 - 2k * p0, p1
    end
    p1
end

function core_hg(x, y, α, γ₀, m, n, isdiagonal::Bool)
    ξ = isdiagonal ? (x + y) / (√2 * γ₀) : x / γ₀
    η = isdiagonal ? (x - y) / (√2 * γ₀) : y / γ₀
    α * exp(-α * (ξ^2 + η^2) / 2) * hermite(m, abs(α) * ξ) * hermite(n, abs(α) * η)
end

function core_hg(x, y, α, γ₀, m, n, x_coeffs, y_coeffs, isdiagonal)
    ξ = isdiagonal ? (x + y) / (√2 * γ₀) : x / γ₀
    η = isdiagonal ? (x - y) / (√2 * γ₀) : y / γ₀
    α * exp(-α * (ξ^2 + η^2) / 2) * hermite(abs(α) * ξ, m, x_coeffs) * hermite(abs(α) * η, n, y_coeffs)
end

function core_hg(x, α, γ₀, n, x_coeffs)
    ξ = x / γ₀
    α * exp(-α * ξ^2 / 2) * hermite(abs(α) * ξ, n, x_coeffs)
end

"""
    hg(x::Real,y::Real,z::Real=0;
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: Real

    function hg(x::AbstractVector{T},y::AbstractVector{T},z::Real=0;
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: Real

    hg(x::AbstractVector{T},y::AbstractVector{T},z::AbstractVector{T};
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: Real

Compute the Hermite-Gaussian mode.

For the first signature, the mode is calculated at point `(x,y,z)`

For the second signature, the mode is calculated over a grid defined by `x` and `y` at a distance `z` from the focus.

For the third signature, the mode is calculated over a grid defined by `x`, `y` and `z`.

The optional keyword arguments are:

`m`: horizontal index

`n`: vertical index

`w0`: beam's waist

`k`: wavenumber
"""
function hg(x::AbstractVector{T}, y::AbstractVector{T}, z::Real=0;
    m::Integer=0, n::Integer=0, w0::Real=1, k::Real=1) where {T<:Real}

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(float(T), w0 / √2)
    k = convert(float(T), k)

    x_coeffs = hermite_coefficients(m)
    y_coeffs = hermite_coefficients(n)

    α = inv(1 + im * z / (k * γ₀^2))
    prefactor = normalization_hg(m, n, γ₀) * cis((m + n) * angle(α))

    @tullio result[j, i] := prefactor * core_hg(x[i], y[j], α, γ₀, m, n, x_coeffs, y_coeffs, false)
end

function hg(x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T};
    m::Integer=0, n::Integer=0, w0::Real=1, k::Real=1) where {T<:Real}

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(float(T), w0 / √2)
    k = convert(float(T), k)

    x_coeffs = hermite_coefficients(m)
    y_coeffs = hermite_coefficients(n)

    function f(x, y, α)
        normalization_hg(m, n, γ₀) * cis((m + n) * angle(α)) * core_hg(x, y, α, γ₀, m, n, x_coeffs, y_coeffs, false)
    end

    @tullio result[j, i, l] := f(x[i], y[j], inv(1 + im * z[l] / (k * γ₀^2)))
end

function hg(x::Real, y::Real, z::Real=0;
    m::Integer=0, n::Integer=0, w0::Real=1, k::Real=1)

    first(hg([x], [y]; z, m, n, k, w0))
end

function hg(x::AbstractVector{T}, z::Real=0;
    n::Integer=0, w0::Real=1, k::Real=1) where {T<:Real}
    #TODO: The normalization only works for z=0

    @assert n ≥ 0

    γ₀ = convert(float(T), w0 / √2)
    k = convert(float(T), k)

    x_coeffs = hermite_coefficients(n)

    α = inv(1 + im * z / (k * γ₀^2))
    prefactor = normalization_hg(n, γ₀) * cis(n * angle(α))

    @tullio result[i] := prefactor * core_hg(x[i], α, γ₀, n, x_coeffs)
end

"""
    diagonal_hg(x::Real,y::Real,z::Real=0;
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: Real

    function diagonal_hg(x::AbstractVector{T},y::AbstractVector{T},z::Real=0;
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: Real

    diagonal_hg(x::AbstractVector{T},y::AbstractVector{T},z::AbstractVector{T};
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: Real

Compute the diagonal Hermite-Gaussian mode.

For the first signature, the mode is calculated at point `(x,y,z)`

For the second signature, the mode is calculated over a grid defined by `x` and `y` at a distance `z` from the focus.

For the third signature, the mode is calculated over a grid defined by `x`, `y` and `z`.

The optional keyword arguments are:

`m`: diagonal index

`n`: antidiagonal index

`w0`: beam's waist

`k`: wavenumber
"""
function diagonal_hg(x::AbstractVector{T}, y::AbstractVector{T}, z::Real=0;
    m::Integer=0, n::Integer=0, w0::Real=1, k::Real=1) where {T<:Real}

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(float(T), w0 / √2)
    k = convert(float(T), k)

    x_coeffs = hermite_coefficients(m)
    y_coeffs = hermite_coefficients(n)

    α = inv(1 + im * z / (k * γ₀^2))
    prefactor = normalization_hg(m, n, γ₀) * cis((m + n) * angle(α))

    @tullio result[j, i] := prefactor * core_hg(x[i], y[j], α, γ₀, m, n, x_coeffs, y_coeffs, true)
end

function diagonal_hg(x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T};
    m::Integer=0, n::Integer=0, w0::Real=1, k::Real=1) where {T<:Real}

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(float(T), w0 / √2)
    k = convert(float(T), k)

    x_coeffs = hermite_coefficients(m)
    y_coeffs = hermite_coefficients(n)

    function f(x, y, α)
        normalization_hg(m, n, γ₀) * cis((m + n) * angle(α)) * core_hg(x, y, α, γ₀, m, n, x_coeffs, y_coeffs, true)
    end

    @tullio result[j, i, l] := f(x[i], y[j], inv(1 + im * z[l] / (k * γ₀^2)))
end

function diagonal_hg(x::Real, y::Real, z::Real=0;
    m::Integer=0, n::Integer=0, w0::Real=1, k::Real=1)

    first(diagonal_hg([x], [y], z, m=m, n=n, k=k, w0=w0))
end

"""
    lens(x,y,fx,fy;k=1)

Output an array containing the phase shift introduced by a lens of focal lengths `fx` and `fy`.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* lens(x,y,fx,fy;k=k)`
"""
function lens(x, y, fx, fy; k=1)
    @tullio result[j, i] := cis(-k * (x[i]^2 / fx + y[j]^2 / fy) / 2)
end

"""
    tilted_lens(x,y,f,ϕ;k=1)

Output an array containing the phase shift introduced by a spherical lens of focal length `f` tilted by an angle `ϕ`.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* tilted_lens(x,y,f,ϕ;k=k)`
"""
function tilted_lens(x, y, f, ϕ; k=1)
    lens(x, y, sec(ϕ) * f, cos(ϕ) * f, k=k)
end