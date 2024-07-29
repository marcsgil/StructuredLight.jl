function ψ!(dest, mode_kernel!, x, y; kwargs...)
    backend = get_backend(dest)
    kernel! = mode_kernel!(backend)
    ψ!(dest, kernel!, x, y; kwargs..., ndrange=size(dest))
end

function ψ!(dest, mode_kernel!, normalization_func, x, y; kwargs...)
    ψ!(dest, mode_kernel!, x, y; kwargs...)
    N = normalization_func(kwargs...)
    dest ./= N
end

function ψ(::Type{T}, args...; kwargs...) where {T}
    dest = Matrix{T}(undef, length(x), length(y))
    ψ!(dest, args...; kwargs...)
end

function ψ(args...; kwargs...)
    
end

#Laguerre-Gaussian modes

function _lg(x, y; p::Integer, l::Integer, γ=one(eltype(x)))
    X = x / γ
    Y = y / γ
    r2 = X^2 + Y^2
    L = abs(l)
    exp(-r2 / 2) * (X + im * sign(l) * Y)^L * laguerre(r2, p, L)
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

function lg(x, y)

end