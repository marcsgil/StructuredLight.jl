complex_type(args...) = promote_type((eltype(arg) for arg ∈ args)...) |> complex

function reshape_3d(A, xs, ys, zs)
    reshape(A, length(xs), length(ys), length(zs))
end

get_size(::Number) = ()
get_size(x) = (length(x))
get_size(x, args...) = (get_size(x)..., get_size(args...)...)


"""
    normalization_lg(p,l,γ=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
function normalization_lg(p, l, γ::T) where {T}
    convert(float(T), √inv(prod(p+1:p+abs(l)) * π) / γ)
end

function lg(x::Real, y::Real, z::Real=zero(x); p=0, l=0, γ=one(x), k=one(x), N=normalization_lg(p, l, γ))
    X = x / γ
    Y = y / γ
    r2 = X^2 + Y^2
    L = abs(l)
    α = inv(1 + im * z / (k * γ^2))
    N * α * exp(-α * r2 / 2 + im * (2p + abs(l)) * angle(α)) * (abs(α) * (X + im * sign(l) * Y))^L * laguerre(abs2(α) * r2, p, L)
end

@kernel function lg_kernel!(dest, x, y, z, p, l, γ, k, N)
    r, s, t = @index(Global, NTuple)
    dest[r, s, t] = lg(x[r], y[s], z[t]; p, l, γ, k, N)
end

function lg!(dest, x, y, z=zero(eltype(x)); p, l, γ=one(eltype(x)), k=one(eltype(x)), N=normalization_lg(p, l, γ))
    backend = get_backend(dest)
    kernel! = lg_kernel!(backend)
    kernel!(dest, x, y, z, p, l, γ, k, N; ndrange=(length(x), length(y), length(z)))
end

function lg(x, y, z=zero(eltype(x)); p, l, γ=one(eltype(x)), k=one(eltype(x)), N=normalization_lg(p, l, γ))
    T = complex_type(x, y, z, p, l, γ, k, N)
    dest = similar(x, T, get_size(x, y, z)...)
    lg!(dest, x, y, z; p, l, γ, k, N)
    dest
end