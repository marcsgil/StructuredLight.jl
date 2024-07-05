function assert_hermite_indices(idxs...)
    for idx ∈ idxs
        @assert idx ≥ 0 "All indices must be non negative"
    end
end

float_type(args...) = float(typeof(sum(first, args)))

"""
    normalization_hg(m,n,γ=1)

Compute the normalization constant for the Hermite-Gaussian modes.
"""
function normalization_hg(m, n, γ::T=1) where {T}
    convert(float(T), inv(γ * √(π * 2^(m + n))) * √(prod(inv, 1:m, init=1) * prod(inv, 1:n, init=1)))
end

function _hg(x, y, s, c; m::Integer=0, n::Integer=0, γ=one(eltype(x)))
    X = (x * c + y * s) / γ
    Y = (-x * s + y * c) / γ

    hermite(X, m) * hermite(Y, n) * exp((-X^2 - Y^2) / 2)
end

function _hg(x, y, α, s, c; m::Integer=0, n::Integer=0, γ=one(eltype(x)))
    X = (x * c + y * s) / γ
    Y = (-x * s + y * c) / γ

    α * exp(α * (-X^2 - Y^2) / 2) * hermite(abs(α) * X, m) * hermite(abs(α) * Y, n)
end

"""
    rotated_hg(x, y; θ, m::Integer=0, n::Integer=0, w=one(eltype(x)))
    rotated_hg(x, y, z; 
               θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a Hermite-Gaussian mode rotated by an angle `θ` (rad).

`x`, `y` and `z` can be numbers or vectors, but `x` and `y` must be always of the same kind.

# Other Arguments:

- `m`: x index when θ=0

- `n`: y index when θ=0

- `w0`: beam's waist

- `k`: wavenumber

# Examples
```jldoctest
julia> rs = LinRange(-5, 5, 256);
julia> ψ₁ = rotated_hg(rs, rs, θ=0.1, m=3, n=2);
julia> ψ₂ = [rotated_hg(x, y, θ=0.1, m=3, n=2) for x in rs, y in rs];
julia> ψ₁ ≈ ψ₂
true
julia> zs = LinRange(0, 1, 32);
julia> ψ₃ = rotated_hg(rs, rs, zs, θ=0.1, m=3, n=2);
julia> ψ₄ = [rotated_hg(x, y, z, θ=0.1, m=3, n=2) for x in rs, y in rs, z in zs];
julia> ψ₃ ≈ ψ₄
true
```

See also [`hg`](@ref), [`diagonal_hg`](@ref), [`lg`](@ref).
"""
function rotated_hg(x::Real, y::Real; θ, m::Integer=0, n::Integer=0, w=one(eltype(x)))
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

    @tullio _[j, k] := N * _hg(x[k], y[j], s, c; m, n, γ)
end

function rotated_hg(x::Real, y::Real, z::Real;
    θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))
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

    @tullio _[j, k] := prefactor * _hg(x[k], y[j], α, s, c; m, n, γ)
end

function rotated_hg(x, y, z;
    θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))
    assert_hermite_indices(m, n)

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)
    s, c = sincos(θ)

    α(z) = inv(1 + im * z / (k * γ^2))

    function prefactor(z)
        cis(-(m + n) * atan(z / (k * γ^2)))
    end
    N = normalization_hg(m, n, γ)

    @tullio _[j, k, l] := N * prefactor(z[l]) * _hg(x[k], y[j], α(z[l]), s, c; m, n, γ)
end

"""
    hg(x, y; m::Integer=0, n::Integer=0, w=one(eltype(x)))
    hg(x, y, z; m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a Hermite-Gaussian mode.

`x`, `y` and `z` can be numbers or vectors, but `x` and `y` must be always of the same kind.

# Other Arguments:

- `m`: x index

- `n`: y index

- `w0`: beam's waist

- `k`: wavenumber

# Examples
```jldoctest
julia> rs = LinRange(-5, 5, 256);
julia> ψ₁ = hg(rs, rs, m=3, n=2);
julia> ψ₂ = [hg(x, y, m=3, n=2) for x in rs, y in rs];
julia> ψ₁ ≈ ψ₂
true
julia> zs = LinRange(0, 1, 32);
julia> ψ₃ = hg(rs, rs, zs, θ=0.1, m=3, n=2);
julia> ψ₄ = [hg(x, y, z, m=3, n=2) for x in rs, y in rs, z in zs];
julia> ψ₃ ≈ ψ₄
true
```

See also [`rotated_hg`](@ref), [`diagonal_hg`](@ref), [`lg`](@ref).
"""
hg(x, y; kwargs...) = rotated_hg(x, y; θ=zero(first(x)), kwargs...)
hg(x, y, z; kwargs...) = rotated_hg(x, y, z; θ=zero(first(x)), kwargs...)

"""
    diagonal_hg(x, y; m::Integer=0, n::Integer=0, w=one(eltype(x)))
    diagonal_hg(x, y, z; m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a diagonal Hermite-Gaussian mode.

`x`, `y` and `z` can be numbers or vectors, but `x` and `y` must be always of the same kind.

# Other Arguments:

- `m`: diagonal index

- `n`: antidiagonal index

- `w0`: beam's waist

- `k`: wavenumber

# Examples
```jldoctest
julia> rs = LinRange(-5, 5, 256);

julia> ψ₁ = diagonal_hg(rs, rs, m=3, n=2);

julia> ψ₂ = [diagonal_hg(x, y, m=3, n=2) for x in rs, y in rs];

julia> ψ₁ ≈ ψ₂
true

julia> zs = LinRange(0, 1, 32);

julia> ψ₃ = diagonal_hg(rs, rs, zs, θ=0.1, m=3, n=2);

julia> ψ₄ = [diagonal_hg(x, y, z, m=3, n=2) for x in rs, y in rs, z in zs];

julia> ψ₃ ≈ ψ₄
true
```

See also [`rotated_hg`](@ref), [`hg`](@ref), [`lg`](@ref).
"""
diagonal_hg(x, y; kwargs...) = rotated_hg(x, y; θ=oftype(float(first(x)), π / 4), kwargs...)
diagonal_hg(x, y, z; kwargs...) = rotated_hg(x, y, z; θ=oftype(float(first(x)), π / 4), kwargs...)

function _lg(x, y; p::Integer, l::Integer, γ=one(eltype(x)))
    X = x / γ
    Y = y / γ
    r2 = X^2 + Y^2
    L = abs(l)
    exp(-r2 / 2) * (X + im * sign(l) * Y)^L * laguerre(r2, p, L)
end

function _lg(x, y, α; p::Integer, l::Integer, γ=one(eltype(x)))
    X = x / γ
    Y = y / γ
    r2 = X^2 + Y^2
    α * exp(-α * r2 / 2) * (abs(α) * (X + im * sign(l) * Y))^abs(l) * laguerre(abs2(α) * r2, p, abs(l))
end

"""
    normalization_lg(p,l,γ=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(p, l, γ=1)
    convert(float(eltype(γ)), √inv(prod(p+1:p+abs(l)) * π) / γ)
end

"""
    lg(x, y; m::Integer=0, n::Integer=0, w=one(eltype(x)))
    lg(x, y, z; m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a diagonal Hermite-Gaussian mode.

`x`, `y` and `z` can be numbers or vectors, but `x` and `y` must be always of the same kind.

# Other Arguments:

- `p`: radial index

- `l`: topological charge

- `w0`: beam's waist

- `k`: wavenumber

# Examples
```jldoctest
julia> rs = LinRange(-5, 5, 256);
julia> ψ₁ = lg(rs, rs, m=3, n=2);
julia> ψ₂ = [lg(x, y, m=3, n=2) for x in rs, y in rs];
julia> ψ₁ ≈ ψ₂
true
julia> zs = LinRange(0, 1, 32);
julia> ψ₃ = lg(rs, rs, zs, θ=0.1, m=3, n=2);
julia> ψ₄ = [lg(x, y, z, m=3, n=2) for x in rs, y in rs, z in zs];
julia> ψ₃ ≈ ψ₄
true
```

See also [`rotated_hg`](@ref), [`hg`](@ref), [`diagonal_hg`](@ref).
"""
function lg(x::Real, y::Real; p::Integer=0, l::Integer=0, w=one(eltype(x)))
    @assert p ≥ 0
    T = float_type(x, y, w)
    γ = convert(T, w / √2)

    normalization_lg(p, l, γ) * _lg(x, y; p, l, γ)
end

function lg(x, y; p::Integer=0, l::Integer=0, w=one(eltype(x)))
    @assert p ≥ 0
    T = float_type(x, y, w)
    γ = convert(T, w / √2)

    N = normalization_lg(p, l, γ)

    @tullio _[j, k] := N * _lg(x[k], y[j]; p, l, γ)
end

function lg(x::Real, y::Real, z::Real;
    p::Integer=0, l::Integer=0, w=one(eltype(x)), k=one(eltype(x)))
    @assert p ≥ 0

    T = float(typeof(sum((x, y, z, w, k))))
    γ = convert(T, w / √2)
    k = convert(T, k)

    α = inv(1 + im * z / (k * γ^2))
    prefactor = normalization_lg(p, l, γ) * cis((2p + abs(l)) * angle(α))
    prefactor * _lg(x, y, α; p, l, γ)
end

function lg(x, y, z::Real;
    p::Integer=0, l::Integer=0, w=one(eltype(x)), k=one(eltype(x)))
    @assert p ≥ 0

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)

    α = inv(1 + im * z / (k * γ^2))
    prefactor = normalization_lg(p, l, γ) * cis((2p + abs(l)) * angle(α))

    @tullio _[j, k] := prefactor * _lg(x[k], y[j], α; p, l, γ)
end

function lg(x, y, z;
    p::Integer=0, l::Integer=0, w=one(eltype(x)), k=one(eltype(x)))
    @assert p ≥ 0

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)

    α(z) = inv(1 + im * z / (k * γ^2))

    function prefactor(z)
        cis(-(2p + abs(l)) * atan(z / (k * γ^2)))
    end
    N = normalization_lg(p, l, γ)

    @tullio _[j, k, l] := N * prefactor(z[l]) * _lg(x[k], y[j], α(z[l]); p, l, γ)
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