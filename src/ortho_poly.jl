"""
    laguerre(x, n::Integer, α::Real=0)

Evaluate the `n` th Laguerre polynomial at `x`.
"""
function laguerre(x, n::Integer, α::Real=0)
    T = typeof(x)
    α = convert(T, α)
    p0, p1 = one(T), -x + (α + 1)
    n == 0 && return p0
    for k = 1:n-1
        p1, p0 = ((2k + α + 1) / (k + 1) - x / (k + 1)) * p1 - (k + α) / (k + 1) * p0, p1
    end
    p1
end

"""
    hermite(x, n::Integer)

Evaluate the `n` th Hermite polynomial at `x`.
"""
function hermite(x, n::Integer)
    T = typeof(x)
    p0, p1 = one(T), 2x
    n == 0 && return p0
    for k = 1:n-1
        p1, p0 = 2x * p1 - 2k * p0, p1
    end
    p1
end