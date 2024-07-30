float_type(args...) = promote_type((eltype(arg) for arg ∈ args)...) |> float

"""
    normalization_lg(p,l,γ=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(p, l, γ::T) where {T}
    convert(float(T), √inv(prod(p+1:p+abs(l)) * π) / γ)
end

function lg(x::Real, y::Real, z::Real=zero(x); p, l, γ=one(x), k=one(x), N=normalization_lg(p, l, γ))
    X = x / γ
    Y = y / γ
    r2 = X^2 + Y^2
    L = abs(l)
    α = inv(1 + im * z / (k * γ^2))
    N * α * exp(-α * r2 / 2) * (abs(α) * (X + im * sign(l) * Y))^L * laguerre(abs2(α) * r2, p, L)
end

"""@kernel function lg_kernel!(x, y, z::Number, p, l, γ, k, N)
    r, s = @index(Global, NTuple)
    lg(x[r], y[s], z; p, l, γ, k, N)
end"""

@kernel function lg_kernel!(x, y, z, p, l, γ, k, N)
    r, s, t = @index(Global, NTuple)
    lg(x[r], y[s], z[t]; p, l, γ, k, N)
end

function lg!(dest, x, y, z=zero(eltype(x)); p, l, γ=one(eltype(x)), k=one(eltype(x)), N=normalization_lg(p, l, γ))
    backend = get_backend(dest)
    kernel! = lg_kernel!(backend)
    kernel!(dest, x, y, z, p, l, γ, k, N; ndrange=size(dest))
end

function lg(x, y, z=zero(eltype(x)); p, l, γ=one(eltype(x)), k=one(eltype(x)), N=normalization_lg(p, l, γ))
    T = float_type(x, y, z, p, l, γ, k, N)
    dest = similar(x, complex(T), length(x), length(y), length(z))
    @show size(dest)
    lg!(dest, x, y, z; p, l, γ, k, N)
    dest
end