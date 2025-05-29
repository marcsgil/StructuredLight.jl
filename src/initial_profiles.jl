complex_type(args...) = promote_type((eltype(arg) for arg ∈ args)...) |> complex

get_size(::Number) = ()
get_size(x) = (length(x),)
get_size(x, args...) = (get_size(x)..., get_size(args...)...)

get_α(z, w, k) = inv(1 + im * 2z / (k * w^2))

"""
    normalization_hg(m,n,w)

Compute the normalization constant for the Hermite-Gaussian modes.
"""
function normalization_hg(m, n, w::T) where {T}
    convert(float(T), √2 * inv(w * √(π * 2^(m + n))) * √(prod(inv, 1:m, init=1) * prod(inv, 1:n, init=1)))
end

"""
    hg(x, y, z=zero(eltype(x)); θ=zero(eltype(x)), m=0, n=0, w=one(eltype(x)), k=one(eltype(x)))

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

See also [`diagonal_hg`](@ref), [`lg`](@ref).
"""
function hg(x::Real, y::Real, z::Real=zero(x); θ=zero(x), m::Integer=0, n::Integer=0, w=one(x), k=one(x), N=normalization_hg(m, n, w))
    s, c = sincos(θ)
    γ = w / oftype(float(w), √2)
    X = (x * c + y * s) / γ
    Y = (-x * s + y * c) / γ
    α = get_α(z, w, k)

    N * α * exp(α * (-X^2 - Y^2) / 2 + im * (m + n) * angle(α)) * hermite(abs(α) * X, m) * hermite(abs(α) * Y, n)
end

@kernel function hg_kernel!(dest, x, y, z, θ, m, n, w, k, N)
    r, s, t = @index(Global, NTuple)
    dest[r, s, t] = hg(x[r], y[s], z[t]; θ, m, n, w, k, N)
end

function hg!(dest, x, y, z=zero(eltype(x)); θ=zero(eltype(x)), m=0, n=0, w=one(eltype(x)), k=one(eltype(x)), N=normalization_hg(m, n, w))
    backend = get_backend(dest)
    kernel! = hg_kernel!(backend)
    ndrange = (length(x), length(y), length(z))
    dest_3d = reshape(dest, ndrange...)
    kernel!(dest_3d, x, y, z, θ, m, n, w, k, N; ndrange)
end

function hg(x, y, z=zero(eltype(x)); θ=zero(eltype(x)), m=0, n=0, w=one(eltype(x)), k=one(eltype(x)), N=normalization_hg(m, n, w), array_constructor::Type{AC}=Array) where {AC<:AbstractArray}
    T = complex_type(x, y, z, θ, m, n, w, k, N)
    dims = get_size(x, y, z)
    N = length(dims)
    dest = AC{T,N}(undef, dims...)
    hg!(dest, x, y, z; θ, m, n, w, k, N)
    dest
end

#= function hg(x, y, z=zero(eltype(x)); θ=zero(eltype(x)), m=0, n=0, w=one(eltype(x)), k=one(eltype(x)), N=normalization_hg(m, n, w))
    T = complex_type(x, y, z, θ, m, n, w, k, N)
    dest = similar(x, T, get_size(x, y, z)...)
    hg!(dest, x, y, z; θ, m, n, w, k, N)
    dest
end =#

diagonal_hg!(dest, x, y, z=zero(eltype(x)); m=0, n=0, w=one(eltype(x)), k=one(eltype(x)), N=normalization_hg(m, n, w)) = hg!(dest, x, y, z; θ=π / 4, m, n, w, k, N)

"""
    diagonal_hg(x, y, z=zero(eltype(x)); m=0, n=0, w=one(eltype(x)), k=one(eltype(x)))

Compute a diagonal Hermite-Gaussian mode. It is calculated by setting `θ=π/4` in [`hg`](@ref).

See also [`lg`](@ref).
"""
diagonal_hg(x, y, z=zero(eltype(x)); m=0, n=0, w=one(eltype(x)), k=one(eltype(x)), N=normalization_hg(m, n, w)) = hg(x, y, z; θ=π / 4, m, n, w, k, N)

"""
    normalization_lg(p,l,w=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(p, l, w::T) where {T}
    convert(float(T), √inv(prod(p+1:p+abs(l)) * π) * √2 / w)
end

"""
    lg(x, y, z=zero(eltype(x)); p=0, l=0, w=one(eltype(x)), k=one(eltype(x)))

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

See also [`hg`](@ref), [`diagonal_hg`](@ref).
"""
function lg(x::Real, y::Real, z::Real=zero(x); p=0, l=0, w=one(x), k=one(x), N=normalization_lg(p, l, w))
    γ = w / oftype(float(w), √2)
    X = x / γ
    Y = y / γ
    r2 = X^2 + Y^2
    L = abs(l)
    α = get_α(z, w, k)

    N * α * exp(-α * r2 / 2 + im * (2p + abs(l)) * angle(α)) * (abs(α) * (X + im * sign(l) * Y))^L * laguerre(abs2(α) * r2, p, L)
end

@kernel function lg_kernel!(dest, x, y, z, p, l, w, k, N)
    r, s, t = @index(Global, NTuple)
    dest[r, s, t] = lg(x[r], y[s], z[t]; p, l, w, k, N)
end

function lg!(dest, x, y, z=zero(eltype(x)); p=0, l=0, w=one(eltype(x)), k=one(eltype(x)), N=normalization_lg(p, l, w))
    backend = get_backend(dest)
    kernel! = lg_kernel!(backend)
    ndrange = (length(x), length(y), length(z))
    dest_3d = reshape(dest, ndrange...)
    kernel!(dest_3d, x, y, z, p, l, w, k, N; ndrange)
end

function lg(x, y, z=zero(eltype(x)); p=0, l=0, w=one(eltype(x)), k=one(eltype(x)), N=normalization_lg(p, l, w))
    T = complex_type(x, y, z, p, l, w, k, N)
    dest = similar(x, T, get_size(x, y, z)...)
    lg!(dest, x, y, z; p, l, w, k, N)
    dest
end

"""
    rectangular_apperture(x, y, a, b)

Determine if points (x, y) are within a rectangular aperture centered at the origin.
The rectangle has width `a` and height `b`.

# Arguments
- `x`: x-coordinates of the points.
- `y`: y-coordinates of the points.
- `a`: Width of the rectangle.
- `b`: Height of the rectangle.

# Returns
- A boolean array indicating whether each point is within the rectangular aperture.
"""
function rectangular_apperture(x, y, a, b)
    @. abs(x) <= a / 2 && abs(y') <= b / 2
end

"""
    square(x, y, l)

Determine if points (x, y) are within a square aperture centered at the origin.
The square has side length `l`.

# Arguments
- `x`: x-coordinates of the points.
- `y`: y-coordinates of the points.
- `l`: Side length of the square.

# Returns
- A boolean array indicating whether each point is within the square aperture.
"""
function square(x, y, l)
    rectangular_apperture(x, y, l, l)
end

"""
    single_slit(x, y, a)

Determine if points (x, y) are within a single slit aperture centered at the origin.
The slit has width `a` and extends infinitely in the y-direction.

# Arguments
- `x`: x-coordinates of the points.
- `y`: y-coordinates of the points.
- `a`: Width of the slit.

# Returns
- A boolean array indicating whether each point is within the single slit aperture.
"""
function single_slit(x, y, a)
    rectangular_apperture(x, y, a, Inf)
end

"""
    double_slit(x, y, a, d)

Determine if points (x, y) are within a double slit aperture centered at the origin.
Each slit has width `a` and they are separated by distance `d`.

# Arguments
- `x`: x-coordinates of the points.
- `y`: y-coordinates of the points.
- `a`: Width of each slit.
- `d`: Distance between the centers of the slits.

# Returns
- A boolean array indicating whether each point is within the double slit aperture.
"""
function double_slit(x, y, a, d)
    @. (abs(x - d / 2) <= a / 2 && abs(y') <= Inf) || (abs(x + d / 2) <= a / 2 && abs(y') <= Inf)
end

"""
    pupil(x, y, r)

Determine if points (x, y) are within a circular aperture centered at the origin.

# Arguments
- `x`: x-coordinates of the points.
- `y`: y-coordinates of the points.
- `r`: Radius of the circle.

# Returns
- A boolean array indicating whether each point is within the circular aperture.
"""
function pupil(x, y, r)
    @. x^2 + y'^2 <= r^2
end

"""
    triangle(x, y, side_length)

Determine if points (x, y) are within an equilateral triangle aperture centered at the origin.
The triangle has side length `side_length`.

# Arguments
- `x`: x-coordinates of the points.
- `y`: y-coordinates of the points.
- `side_length`: Side length of the triangle.

# Returns
- A boolean array indicating whether each point is within the triangular aperture.
"""
function triangle(x, y, side_length)
    sqrt3 = √3
    @. y' > -side_length / 2 / sqrt3 && abs(x) < -y' / sqrt3 + side_length / 3
end

"""linear_combination(f, c) = (args...; kwargs...) -> map(c, f) do c, f
    c * f(args...; kwargs...)
end |> sum

function get_paralelized_function_2D(f)
    @kernel function kernel_2D!(dest, x, y)
        i, j = @index(Global, NTuple)
        dest[i, j] = f(x[i], y[j])
    end

    function f!(dest, x, y)
        backend = get_backend(dest)
        kernel! = kernel_2D!(backend)
        kernel!(dest, x, y; ndrange=size(dest))
    end

    function f_parallel(x, y)
        T = complex_type(x, y)
        dest = similar(x, T, size(x, 1), size(y, 1))
        f!(dest, x, y)
        dest
    end

    f_parallel, f!
end

function linear_combination_2D!(dest, fs, cs, x, y)
    f_parallel, f! = get_paralelized_function_2D(linear_combination(fs, cs))
    f_parallel, f!

end"""
