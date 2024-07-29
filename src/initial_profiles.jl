function assert_hermite_indices(idxs...)
    for idx ∈ idxs
        @assert idx ≥ 0 "All indices must be non negative"
    end
end

float_type(args...) = promote_type((eltype(arg) for arg ∈ args)...) |> float

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

@kernel function hg_kernel!(dest, x, y, s, c, m, n, γ)
    j, k = @index(Global, NTuple)
    dest[j, k] = _hg(x[j], y[k], s, c; m, n, γ)
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

- `w`: beam's waist

- `k`: wavenumber

# Examples

```jldoctest
rs = LinRange(-5, 5, 256)
zs = LinRange(0, 1, 32)

# just x and y
ψ₁ = rotated_hg(rs, rs, θ = 0.1, m=3, n=2)
ψ₂ = [rotated_hg(x, y, θ = 0.1, m=3, n=2) for x in rs, y in rs]

# x, y and z
ψ₃ = rotated_hg(rs, rs, zs, θ = 0.1, m=3, n=2)
ψ₄ = [rotated_hg(x, y, z, θ = 0.1, m=3, n=2) for x in rs, y in rs, z ∈ zs]

ψ₁ ≈ ψ₂ && ψ₃ ≈ ψ₄

# output

true
```

See also [`hg`](@ref), [`diagonal_hg`](@ref), [`lg`](@ref).
"""
function rotated_hg(x::Real, y::Real; θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), include_normalization=true)
    assert_hermite_indices(m, n)
    T = float_type(x, y, w, θ)
    γ = convert(T, w / √2)
    s, c = sincos(θ)

    result = _hg(x, y, s, c; m, n, γ)

    if include_normalization
        result *= normalization_hg(m, n, γ)
    end

    result
end

function rotated_hg(x, y; θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), include_normalization=true)
    assert_hermite_indices(m, n)
    T = float_type(x, y, w, θ)
    γ = convert(T, w / √2)
    s, c = sincos(θ)

    result = similar(x, length(x), length(y))
    backend = get_backend(result)
    kernel! = hg_kernel!(backend, 256)
    kernel!(result, x, y, s, c, m, n, γ; ndrange=size(result))

    if include_normalization
        N = normalization_hg(m, n, γ)
        result *= N
    end

    result
end

function rotated_hg(x::Real, y::Real, z::Real;
    θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)), include_normalization=true)
    assert_hermite_indices(m, n)

    T = float_type(x, y, z, w, k, θ)
    γ = convert(T, w / √2)
    k = convert(T, k)
    s, c = sincos(θ)

    α = inv(1 + im * z / (k * γ^2))
    if include_normalization
        prefactor = normalization_hg(m, n, γ) * cis((m + n) * angle(α))
    else
        prefactor = cis((m + n) * angle(α))
    end
    prefactor * _hg(x, y, α, s, c; m, n, γ)
end

function rotated_hg(x, y, z::Real;
    θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)), include_normalization=true)
    assert_hermite_indices(m, n)

    T = float_type(x, y, z, w, k, θ)
    γ = convert(T, w / √2)
    k = convert(T, k)
    s, c = sincos(θ)

    α = inv(1 + im * z / (k * γ^2))
    if include_normalization
        prefactor = normalization_hg(m, n, γ) * cis((m + n) * angle(α))
    else
        prefactor = cis((m + n) * angle(α))
    end

    @tullio _[j, k] := prefactor * _hg(x[j], y[k], α, s, c; m, n, γ)
end

function rotated_hg(x, y, z;
    θ, m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)), include_normalization=true)
    assert_hermite_indices(m, n)

    T = float_type(x, y, z, w, k, θ)
    γ = convert(T, w / √2)
    k = convert(T, k)
    s, c = sincos(θ)

    α(z) = inv(1 + im * z / (k * γ^2))

    function prefactor(z)
        cis(-(m + n) * atan(z / (k * γ^2)))
    end

    @tullio result[j, k, l] := prefactor(z[l]) * _hg(x[j], y[k], α(z[l]), s, c; m, n, γ)

    if include_normalization
        N = normalization_hg(m, n, γ)
        result *= N
    end

    result
end

"""
    hg(x, y; m::Integer=0, n::Integer=0, w=one(eltype(x)))
    hg(x, y, z; m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a Hermite-Gaussian mode.

`x`, `y` and `z` can be numbers or vectors, but `x` and `y` must be always of the same kind.

# Other Arguments:

- `m`: x index

- `n`: y index

- `w`: beam's waist

- `k`: wavenumber

# Examples

```jldoctest
rs = LinRange(-5, 5, 256)
zs = LinRange(0, 1, 32)

# just x and y
ψ₁ = hg(rs, rs, m=3, n=2)
ψ₂ = [hg(x, y, m=3, n=2) for x in rs, y in rs]

# x, y and z
ψ₃ = hg(rs, rs, zs, m=3, n=2)
ψ₄ = [hg(x, y, z, m=3, n=2) for x in rs, y in rs, z ∈ zs]

ψ₁ ≈ ψ₂ && ψ₃ ≈ ψ₄

# output

true
```

See also [`rotated_hg`](@ref), [`diagonal_hg`](@ref), [`lg`](@ref).
"""
hg(x, y; kwargs...) = rotated_hg(x, y; θ=zero(eltype(x)), kwargs...)
hg(x, y, z; kwargs...) = rotated_hg(x, y, z; θ=zero(eltype(x)), kwargs...)

"""
    diagonal_hg(x, y; m::Integer=0, n::Integer=0, w=one(eltype(x)))
    diagonal_hg(x, y, z; m::Integer=0, n::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a diagonal Hermite-Gaussian mode.

`x`, `y` and `z` can be numbers or vectors, but `x` and `y` must be always of the same kind.

# Other Arguments:

- `m`: diagonal index

- `n`: antidiagonal index

- `w`: beam's waist

- `k`: wavenumber

# Examples

```jldoctest
rs = LinRange(-5, 5, 256)
zs = LinRange(0, 1, 32)

# just x and y
ψ₁ = diagonal_hg(rs, rs, m=3, n=2)
ψ₂ = [diagonal_hg(x, y, m=3, n=2) for x in rs, y in rs]

# x, y and z
ψ₃ = diagonal_hg(rs, rs, zs, m=3, n=2)
ψ₄ = [diagonal_hg(x, y, z, m=3, n=2) for x in rs, y in rs, z ∈ zs]

ψ₁ ≈ ψ₂ && ψ₃ ≈ ψ₄

# output

true
```

See also [`rotated_hg`](@ref), [`hg`](@ref), [`lg`](@ref).
"""
diagonal_hg(x, y; kwargs...) = rotated_hg(x, y; θ=convert(eltype(x), π / 4), kwargs...)
diagonal_hg(x, y, z; kwargs...) = rotated_hg(x, y, z; θ=convert(eltype(x), π / 4), kwargs...)

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
    L = abs(l)
    α * exp(-α * r2 / 2) * (abs(α) * (X + im * sign(l) * Y))^L * laguerre(abs2(α) * r2, p, L)
end

"""
    normalization_lg(p,l,γ=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(p, l, γ=1)
    convert(float(eltype(γ)), √inv(prod(p+1:p+abs(l)) * π) / γ)
end

@kernel function lg_kernel!(dest, x, y, p, l, γ)
    j, k = @index(Global, NTuple)
    dest[j, k] = _lg(x[j], y[k]; p, l, γ)
end
"""
    lg(x, y; p::Integer=0, l::Integer=0, w=one(eltype(x)))
    lg(x, y, z; p::Integer=0, l::Integer=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a diagonal Hermite-Gaussian mode.

`x`, `y` and `z` can be numbers or vectors, but `x` and `y` must be always of the same kind.

# Other Arguments:

- `p`: radial index

- `l`: topological charge

- `w`: beam's waist

- `k`: wavenumber

# Examples

```jldoctest
rs = LinRange(-5, 5, 256)
zs = LinRange(0, 1, 32)

# just x and y
ψ₁ = lg(rs, rs, p=1, l=2)
ψ₂ = [lg(x, y, p=1, l=2) for x in rs, y in rs]

# x, y and z
ψ₃ = lg(rs, rs, zs, p=1, l=2)
ψ₄ = [lg(x, y, z, p=1, l=2) for x in rs, y in rs, z ∈ zs]

ψ₃ ≈ ψ₄

# output

true
```

See also [`rotated_hg`](@ref), [`hg`](@ref), [`diagonal_hg`](@ref).
"""
function lg(x::Real, y::Real; p::Integer=0, l::Integer=0, w=one(eltype(x)), include_normalization=true)
    @assert p ≥ 0
    T = float_type(x, y, w)
    γ = convert(T, w / √2)

    result = _lg(x, y; p, l, γ)

    if include_normalization
        result *= normalization_lg(p, l, γ)
    end

    result
end

function lg(x, y; p::Integer=0, l::Integer=0, w=one(eltype(x)), include_normalization=true)
    @assert p ≥ 0
    T = float_type(x, y, w)
    γ = convert(T, w / √2)

    result = similar(x, complex(T), length(x), length(y))
    backend = get_backend(result)
    kernel! = lg_kernel!(backend, 256)
    kernel!(result, x, y, p, l, γ; ndrange=size(result))

    if include_normalization
        N = normalization_lg(p, l, γ)
        result *= N
    end

    result
end

function lg(x::Real, y::Real, z::Real;
    p::Integer=0, l::Integer=0, w=one(eltype(x)), k=one(eltype(x)), include_normalization=true)
    @assert p ≥ 0

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)

    α = inv(1 + im * z / (k * γ^2))
    prefactor = cis((2p + abs(l)) * angle(α))

    result = prefactor * _lg(x, y, α; p, l, γ)

    if include_normalization
        N = normalization_lg(p, l, γ)
        result *= N
    end

    result
end

function lg(x, y, z::Real;
    p::Integer=0, l::Integer=0, w=one(eltype(x)), k=one(eltype(x)), include_normalization=true)
    @assert p ≥ 0

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)

    α = inv(1 + im * z / (k * γ^2))

    if include_normalization
        prefactor = normalization_lg(p, l, γ) * cis((2p + abs(l)) * angle(α))
    else
        prefactor = cis((2p + abs(l)) * angle(α))
    end

    @tullio _[j, k] := prefactor * _lg(x[j], y[k], α; p, l, γ)
end

function lg(x, y, z;
    p::Integer=0, l::Integer=0, w=one(eltype(x)), k=one(eltype(x)), include_normalization=true)
    @assert p ≥ 0

    T = float_type(x, y, z, w, k)
    γ = convert(T, w / √2)
    k = convert(T, k)

    α(z) = inv(1 + im * z / (k * γ^2))

    function prefactor(z)
        cis(-(2p + abs(l)) * atan(z / (k * γ^2)))
    end


    @tullio result[j, k, m] := prefactor(z[m]) * _lg(x[j], y[k], α(z[m]); p, l, γ)

    if include_normalization
        N = normalization_lg(p, l, γ)
        result *= N
    end

    result
end

"""
    lens(x,y,fx,fy;k=1)

Output an array containing the phase shift introduced by a lens of focal lengths `fx` and `fy`.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* lens(x,y,fx,fy;k=k)`
"""
function lens(x, y, fx, fy; k=1)
    @tullio result[i, j] := cis(-k * (x[i]^2 / fx + y[j]^2 / fy) / 2)
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