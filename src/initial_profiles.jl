"""
    binomial(x::Number,y::Number)

Compute the binomial coefficient for noninteger `x` and `y`.
"""
Base.binomial(x::Number, y::Number) = inv((x+1) * beta(x-y+1, y+1))

"""
    laguerre_coefficients(n::Integer,α=0)

Compute the coefficients of the nth generalized Laguerre Polynomial.
"""
function laguerre_coefficients(n,α=0)
    ntuple(i->-(-1)^i*binomial(n+α,n-i+1)/factorial(i-1),n+1)
end

"""
    normalization_lg(;p,l,γ₀=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
normalization_lg(;p,l,γ₀=1) = 1/(γ₀*√( oftype(float(γ₀),π)*prod(p+1:p+abs(l))))


function core_lg(x,y,α,γ₀,l,coefs)
    r2 = (x^2 + y^2)/γ₀^2
    α*exp(-α*r2/2)*(abs(α)*(x+im*sign(l)*y)/γ₀)^abs(l)*evalpoly(abs2(α)*r2,coefs)
end

"""
    lg(xs::AbstractVector{T},ys::AbstractVector{T},z::Real=0;
        p::Integer=0,l::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

Compute the Laguerre-Gaussian mode over a cartesian grid defined by `xs` and `ys`. One may give a distance `z` away from the focus.

The optional keyword arguments are:

`p`: radial index

`l`: topological charge

`w0`: beam's waist

`k`: wavenumber
"""
function lg(xs::AbstractVector{T},ys::AbstractVector{T},z::Real=0;
    p::Integer=0,l::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

    @assert p ≥ 0

    γ₀ = convert(T,w0/√2)
    k = convert(T,k)

    coefs = laguerre_coefficients(p,convert(T,abs(l)))

    α = 1/(1+im*z/(k*γ₀^2))
    prefactor = normalization_lg(p=p,l=l,γ₀=γ₀) * cis((2p+abs(l))*angle(α))

    ThreadsX.map(r->prefactor*core_lg(r...,α,γ₀,l,coefs), Iterators.product(xs,ys))
end

"""
    lg(xs::AbstractVector{T},ys::AbstractVector{T},zs::AbstractVector{T};
        p::Integer=0,l::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

Compute the Laguerre-Gaussian mode over a cartesian grid defined by `xs`, `ys` and `zs`.

The optional keyword arguments are:

`p`: radial index

`l`: topological charge

`w0`: beam's waist

`k`: wavenumber
"""
function lg(xs::AbstractVector{T},ys::AbstractVector{T},zs::AbstractVector{T};
    p::Integer=0,l::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

    @assert p ≥ 0

    γ₀ = convert(T,w0/√2)
    k = convert(T,k)

    coefs = laguerre_coefficients(p,convert(T,abs(l)))

    function f(x,y,α)
        normalization_lg(p=p,l=l,γ₀=γ₀) * cis((2p+abs(l))*angle(α)) * core_lg(x,y,α,γ₀,l,coefs)
    end

    ThreadsX.map(z -> map( rs -> f(rs...,inv(1+im*z/(k*γ₀^2))), Iterators.product(xs,ys) ), zs) |> stack
end


function hermite_coefficients(n)
    if iseven(n)
        ntuple(l -> - factorial(n) * (-1) ^ (n ÷ 2 - l) / ( factorial(2l-2) * factorial( n÷2 - l + 1 ) ) |> Integer, n÷2+1)
    else
        ntuple(l -> - factorial(n) * (-1) ^ (n ÷ 2 - l) / ( factorial(2l-1) * factorial( n÷2 - l + 1 ) ) |> Integer, n÷2+1)
    end
end

"""
    hermite(x,n,coefs)

Evaluate the `n` th Hermite polynomial at `x`, given the coefficients `coefs`.
"""
hermite(x,n,coefs) = iseven(n) ? evalpoly(4x^2,coefs) : 2x*evalpoly(4x^2,coefs)

"""
    normalization_hg(;m,n,γ₀=1)

Compute the normalization constant for the Laguerre-Gaussian modes.
"""
normalization_hg(;m,n,γ₀=1) = 1/(γ₀*√( oftype(float(γ₀),π)*2^(m+n)*factorial(n)*factorial(m)))

function core_hg(x,y,α,γ₀,m,n,x_coefs,y_coefs,isdiagonal)
    ξ = isdiagonal ? (x+y)/( √2 * γ₀) : x/γ₀
    η = isdiagonal ? (x-y)/( √2 * γ₀) : y/γ₀
    α*exp(-α*(ξ^2+η^2)/2)*hermite(abs(α)*ξ,m,x_coefs)*hermite(abs(α)*η,n,y_coefs)
end

"""
    function hg(xs::AbstractVector{T},ys::AbstractVector{T},z::Real=0;
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

Compute the Hermite-Gaussian mode over a cartesian grid defined by `xs` and `ys`. One may give a distance `z` away from the focus.

The optional keyword arguments are:

`m`: horizontal index

`n`: vertical index

`w0`: beam's waist

`k`: wavenumber
"""
function hg(xs::AbstractVector{T},ys::AbstractVector{T},z::Real=0;
    m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(T,w0/√2)
    k = convert(T,k)

    x_coefs = hermite_coefficients(m)
    y_coefs = hermite_coefficients(n)

    α = inv(1+im*z/(k*γ₀^2))
    prefactor = normalization_hg(m=m,n=n,γ₀=γ₀) * cis((m+n)*angle(α))

    ThreadsX.map(r->prefactor*core_hg(r...,α,γ₀,m,n,x_coefs,y_coefs,false), Iterators.product(xs,ys))
end

"""
    hg(xs::AbstractVector{T},ys::AbstractVector{T},zs::AbstractVector{T};
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

Compute the Hermite-Gaussian mode over a cartesian grid defined by `xs`, `ys` and `zs`.

The optional keyword arguments are:

`m`: horizontal index

`n`: vertical index

`w0`: beam's waist

`k`: wavenumber
"""
function hg(xs::AbstractVector{T},ys::AbstractVector{T},zs::AbstractVector{T};
    m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(T,w0/√2)
    k = convert(T,k)

    x_coefs = hermite_coefficients(m)
    y_coefs = hermite_coefficients(n)

    function f(x,y,α)
        normalization_hg(m=m,n=n,γ₀=γ₀) * cis((m+n)*angle(α)) * core_hg(x,y,α,γ₀,m,n,x_coefs,y_coefs,false)
    end

    ThreadsX.map(z -> map( rs -> f(rs...,inv(1+im*z/(k*γ₀^2))), Iterators.product(xs,ys) ), zs) |> stack
end

"""
    function diagonal_hg(xs::AbstractVector{T},ys::AbstractVector{T},z::Real=0;
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

Compute the diagonal Hermite-Gaussian mode over a cartesian grid defined by `xs` and `ys`. One may give a distance `z` away from the focus.

The optional keyword arguments are:

`m`: diagonal index

`n`: antidiagonal index

`w0`: beam's waist

`k`: wavenumber
"""
function diagonal_hg(xs::AbstractVector{T},ys::AbstractVector{T},z::Real=0;
    m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(T,w0/√2)
    k = convert(T,k)

    x_coefs = hermite_coefficients(m)
    y_coefs = hermite_coefficients(n)

    α = inv(1+im*z/(k*γ₀^2))
    prefactor = normalization_hg(m=m,n=n,γ₀=γ₀) * cis((m+n)*angle(α))

    ThreadsX.map(r->prefactor*core_hg(r...,α,γ₀,m,n,x_coefs,y_coefs,true), Iterators.product(xs,ys))
end

"""
    diagonal_hg(xs::AbstractVector{T},ys::AbstractVector{T},zs::AbstractVector{T};
        m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

Compute the Hermite-Gaussian mode over a cartesian grid defined by `xs`, `ys` and `zs`.

The optional keyword arguments are:

`m`: diagonal index

`n`: antidiagonal index

`w0`: beam's waist

`k`: wavenumber
"""
function diagonal_hg(xs::AbstractVector{T},ys::AbstractVector{T},zs::AbstractVector{T};
    m::Integer=0,n::Integer=0,w0::Real=1,k::Real=1) where T<: AbstractFloat

    @assert m ≥ 0
    @assert n ≥ 0

    γ₀ = convert(T,w0/√2)
    k = convert(T,k)

    x_coefs = hermite_coefficients(m)
    y_coefs = hermite_coefficients(n)

    function f(x,y,α)
        normalization_hg(m=m,n=n,γ₀=γ₀) * cis((m+n)*angle(α)) * core_hg(x,y,α,γ₀,m,n,x_coefs,y_coefs,true)
    end

    ThreadsX.map(z -> map( rs -> f(rs...,inv(1+im*z/(k*γ₀^2))), Iterators.product(xs,ys) ), zs) |> stack
end

function lens(xs,ys,fx,fy;k=1)
    ThreadsX.map(r -> cis( -k * (r[1]^2/fx + r[2]^2/fy) / 2 ), Iterators.product(xs,ys))
end

function tilted_lens(xs,ys,f,ϕ;k=1)
    lens(xs,ys,sec(ϕ)*f,cos(ϕ)*f,k=k)
end