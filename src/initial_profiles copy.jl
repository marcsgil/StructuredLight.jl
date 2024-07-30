float_type(args...) = promote_type((eltype(arg) for arg ∈ args)...) |> float

γ(z, γ₀, k) = √(γ₀^2 + (z / (k * γ₀))^2)
K(z, γ₀, k) = k * z / (z^2 + (k * γ₀^2)^2)

"""
    normalization_lg(p,l,γ=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(p, l, γ=1)
    convert(float(eltype(γ)), √inv(prod(p+1:p+abs(l)) * π) / γ)
end

function lg(x::Number, y::Number; p=0, l=0, γ₀=one(x))
    X = x / γ₀
    Y = y / γ₀
    r2 = X^2 + Y^2
    L = abs(l)
    exp(-r2 / 2) * (X + im * sign(l) * Y)^L * laguerre(r2, p, L)
end

function lg(x::Number, y::Number, z::Number; p=0, l=0, γ₀=one(eltype(x)), k=one(eltype(x)))
    γz = γ(z, γ₀, k)
    Kz = K(z, γ₀, k)

    X = x / γz
    Y = y / γz
    r2 = X^2 + Y^2
    L = abs(l)
    order_p1 = 2p + abs(l) + 1

    exp(-r2 * (1 + im * Kz) / 2 - im * order_p1 * atan(z / k / γ₀^2)) * (X + im * sign(l) * Y)^L * laguerre(r2, p, L)
end

@kernel function lg_kernel!(dest, x, y, p, l, γ₀)
    r, s = @index(Global, NTuple)
    dest[r, s] = lg(x[r], y[s]; p, l, γ₀)
end

@kernel function lg_kernel!(dest, x, y, z::Number, p, l, γ₀, k)
    r, s = @index(Global, NTuple)
    dest[r, s] = lg(x[r], y[s], z; p, l, γ₀, k)
end

@kernel function lg_kernel!(dest, x, y, z, p, l, γ₀, k)
    r, s, t = @index(Global, NTuple)
    dest[r, s, t] = lg(x[r], y[s], z[t]; p, l, γ₀, k)
end

function lg!(dest, x, y; p=0, l=0, γ₀=one(eltype(x)))
    backend = get_backend(dest)
    kernel! = lg_kernel!(backend)
    kernel!(dest, x, y, p, l, γ₀, ndrange=size(dest))
end

function lg!(dest, x, y, z::Number; p=0, l=0, γ₀=one(eltype(x)), k=one(eltype(x)))
    backend = get_backend(dest)
    kernel! = lg_kernel!(backend)
    kernel!(dest, x, y, z, p, l, γ₀, k, ndrange=size(dest))
end

function lg!(dest, x, y, z; p=0, l=0, γ₀=one(eltype(x)), k=one(eltype(x)))
    backend = get_backend(dest)
    kernel! = lg_kernel!(backend)
    kernel!(dest, x, y, z, p, l, γ₀, k, ndrange=size(dest))
end

function lg(x, y; p=0, l=0, γ₀=one(eltype(x)))
    T = complex(float_type(x, y, γ₀))
    dest = similar(x, T, length(x), length(y))
    lg!(dest, x, y; p, l, γ₀)
    dest
end

function lg(x, y, z::Number; p=0, l=0, γ₀=one(eltype(x)), k=one(eltype(x)))
    T = complex(float_type(x, y, γ₀, k))
    dest = similar(x, T, length(x), length(y))
    lg!(dest, x, y, z; p, l, γ₀, k)
    dest
end

function lg(x, y, z; p=0, l=0, γ₀=one(eltype(x)), k=one(eltype(x)))
    T = complex(float_type(x, y, γ₀, k))
    dest = similar(x, T, length(x), length(y), length(z))
    lg!(dest, x, y, z; p, l, γ₀, k)
    dest
end