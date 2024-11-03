"""
    find_zero(f, x₁, x₂; tol=1e-6, maxiter=10^3)

Finds a zero of the function `f` within the interval `[x₁, x₂]` using the bisection method.

`f(x₁)` and `f(x₂)` should have opposite signs, so that the method is guaranteed to converge if `f` is continuous.

# Examples
```jldoctest
using StructuredLight: find_zero
f(x) = x^2 - 4
zero = find_zero(f, 0, 5)
isapprox(zero, 2; atol=1e-6)

# output

true
```
"""
function find_zero(f, x₁, x₂; tol=1e-6, maxiter=10^3)
    @assert f(x₁) * f(x₂) ≤ 0 "Error: f(x₁) and f(x₂) must have opposite signs"

    isapprox(f(x₁), 0; atol=tol) && return 2x₁ / 2
    isapprox(f(x₂), 0; atol=tol) && return 2x₂ / 2

    for _ ∈ 1:maxiter
        x = (x₁ + x₂) / 2
        y = f(x)

        if abs(y) < tol
            return x
        end

        if f(x) * f(x₁) < 0
            x₂ = x
        else
            x₁ = x
        end
    end
end

"""
    inverse(f, y, x₁, x₂; tol=1e-6, maxiter=10^3)

Finds the inverse of the function `f` at the value `y` within the interval `[x₁, x₂]` using the bisection method.

`f(x₁)` and `f(x₂)` should have opposite signs, so that the method is guaranteed to converge if `f` is continuous.

# Examples
```jldoctest
using StructuredLight: inverse
f(x) = 2x + 1
inv = inverse(f, 3, 0, 5)
isapprox(inv, 1; atol=1e-6)

# output

true
```
"""
function inverse(f, y, x₁, x₂)
    find_zero(x -> f(x) - y, x₁, x₂)
end

"""
    closest_index(number, grid)

Finds the index of the closest value to `number` in the uniform `grid`.

# Examples
```jldoctest
using StructuredLight: closest_index
grid = [1, 3, 5, 7, 9]
idx = closest_index(3.2, grid)
idx == 2

# output

true
```
"""
function closest_index(number, grid)
    step = grid[2] - grid[1]
    idx = round(Int, (number - grid[1]) / step) + 1
    clamp(idx, 1, length(grid))
end

"""
    zero_order_interpolation(x, xs, ys)

Performs zero-order interpolation to find the value in `ys` corresponding to the closest value to `x` in `xs`.

# Examples
```jldoctest
using StructuredLight: zero_order_interpolation
xs = [1, 2, 3, 4, 5]
ys = [10, 20, 30, 40, 50]
value = zero_order_interpolation(2.7, xs, ys)
value == 30

# output

true
```
"""
function zero_order_interpolation(x, xs, ys)
    i = closest_index(x, xs)
    ys[i]
end

abs_div(x, y) = abs(x / y)


"""
    max_abs_div(x, y [, base_size])

Finds the maximum absolute value of the division of the elements of `x` and `y`.

If is provided, the function will run a divide and conquer type algorithm to parallelize the computation, using `base_size` as the base case.
Otherwise, the function will use a simple loop to find the maximum value.

# Examples
```jldoctest
using StructuredLight: max_abs_div
x = [-1, 2, -3]
y = [4, 5, 6]
result = max_abs_div(x, y)
result == 0.5

# output

true
```

```jldoctest
using StructuredLight: max_abs_div
x = rand(10^3)
y = rand(10^3) .+ 1
result = max_abs_div(x, y, 10^2)
result ≤ 1

# output

true
```
"""
function max_abs_div(x, y)
    M = real(zero(eltype(x)))

    for (x, y) ∈ zip(x, y)
        M = max(M, abs_div(x, y))
    end

    M
end

function max_abs_div(x::Union{Array,SubArray}, y::Union{Array,SubArray}, base_size)
    if length(x) < base_size || length(y) < base_size
        return max_abs_div(x, y)
    end

    x₁ = view(x, 1:length(x)÷2)
    x₂ = view(x, length(x)÷2+1:length(x))
    y₁ = view(y, 1:length(y)÷2)
    y₂ = view(y, length(y)÷2+1:length(y))

    M₁ = Threads.@spawn max_abs_div(x₁, y₁, base_size)
    M₂ = Threads.@spawn max_abs_div(x₂, y₂, base_size)

    max(fetch(M₁), fetch(M₂))
end

function max_abs_div(x, y, base_size)
    mapreduce(abs_div, max, x, y)
end

normalize_angle(angle) = mod2pi(angle + π) - π

abstract type HologramMethod end

"""
    BesselJ1 <: HologramMethod

A hologram generation method based on the inverse of the BesselJ1 function. 
Corresponds to method of 3 of [1] and F of [2].
It is not compatible with GPUs.

See also [`Simple`](@ref).

# References

[1] Victor Arrizón, Ulises Ruiz, Rosibel Carrada, and Luis A. González,
"Pixelated phase computer holograms for the accurate encoding of scalar complex fields,"
J. Opt. Soc. Am. A 24, 3500-3507 (2007)

[2] Thomas W. Clark, Rachel F. Offer, Sonja Franke-Arnold, Aidan S. Arnold, and Neal Radwell, 
"Comparison of beam generation techniques using a phase only spatial light modulator," 
Opt. Express 24, 6249-6264 (2016)
"""
struct BesselJ1 <: HologramMethod end

"""
    Simple <: HologramMethod

A hologram generation with a naive amplitude modulation.
Method of 3 of [1] and F of [2].
It is compatible with GPUs.

See also [`BesselJ1`](@ref).

# References

[1] Victor Arrizón, Ulises Ruiz, Rosibel Carrada, and Luis A. González,
"Pixelated phase computer holograms for the accurate encoding of scalar complex fields,"
J. Opt. Soc. Am. A 24, 3500-3507 (2007)

[2] Thomas W. Clark, Rachel F. Offer, Sonja Franke-Arnold, Aidan S. Arnold, and Neal Radwell, 
"Comparison of beam generation techniques using a phase only spatial light modulator," 
Opt. Express 24, 6249-6264 (2016)
"""
struct Simple <: HologramMethod end

const x_min_besselj1 = 0
const x_max_besselj1 = 0.5818
const y_min_besselj1 = 0
const y_max_besselj1 = 1.82337
const xs_besselj1 = LinRange(x_min_besselj1, x_max_besselj1, 1024)
const ys_besselj1 = inverse.(x -> besselj1(x), xs_besselj1, y_min_besselj1, y_max_besselj1)

inverse_besselj1(x) = zero_order_interpolation(x, xs_besselj1, ys_besselj1)

@kernel function hologram_kernel!(dest, desired, incoming, two_pi_modulation, M, x_period, y_period, ::Type{BesselJ1})
    i, j = @index(Global, NTuple)

    relative = desired[i, j] / incoming[i, j] / M
    ψ = (inverse_besselj1(x_max_besselj1 * abs(relative))
         *
         sin(2π * (i / x_period + j / y_period) + angle(relative)))

    dest[i, j] = round(two_pi_modulation * 0.586 * (ψ / y_max_besselj1 + 1) / 2)
end

@kernel function hologram_kernel!(dest, desired, incoming, two_pi_modulation, M, x_period, y_period, ::Type{Simple})
    i, j = @index(Global, NTuple)

    relative = desired[i, j] / incoming[i, j] / M
    ψ = abs(relative) * normalize_angle(angle(relative) + 2π * (i / x_period + j / y_period))

    dest[i, j] = round(two_pi_modulation * (ψ + π) / 2π)
end

"""
    generate_hologram!(dest, desired, incoming, two_pi_modulation, 
        x_period, y_period, method::Type{T}=BesselJ1) where {T<:HologramMethod}

Same as [`generate_hologram`](@ref), but writes the result to `dest`.
"""
function generate_hologram!(dest, desired, incoming, two_pi_modulation, x_period, y_period, method::Type{T}=BesselJ1) where {T<:HologramMethod}
    M = max_abs_div(desired, incoming, length(desired) ÷ (2 * Threads.nthreads()))

    backend = get_backend(dest)
    kernel! = hologram_kernel!(backend)
    kernel!(dest, desired, incoming, two_pi_modulation, M, x_period, y_period, method; ndrange=size(dest))
end


"""
    generate_hologram(desired, incoming, two_pi_modulation, 
        x_period, y_period, method::Type{T}=BesselJ1) where {T<:HologramMethod}

Returns a hologram built from the desired and incoming images, using the specified `method`.

`two_pi_modulation` is the maximum value of the hologram (should be an integer between 0 and 255). 
`x_period` and `y_period` are the periods, in units of pixels, of the hologram in the x and y directions.
The avaliable methods are [`BesselJ1`](@ref) and [`Simple`](@ref).

Se also [`generate_hologram!`](@ref).

# Examples

```julia
rs = LinRange(-4.0f0, 4, 1024)

desired = lg(rs, rs, l=1, p=1)
incoming = fill(one(eltype(desired)), size(desired)...) # plane wave

# Simple method
dest = generate_hologram(desired, incoming, rs, rs, 255, 50, -50, Simple) 

# BesselJ1 method
dest = generate_hologram(desired, incoming, rs, rs, 255, 50, -50, BesselJ1)
```

# References

[1] Victor Arrizón, Ulises Ruiz, Rosibel Carrada, and Luis A. González,
"Pixelated phase computer holograms for the accurate encoding of scalar complex fields,"
J. Opt. Soc. Am. A 24, 3500-3507 (2007)

[2] Thomas W. Clark, Rachel F. Offer, Sonja Franke-Arnold, Aidan S. Arnold, and Neal Radwell, 
"Comparison of beam generation techniques using a phase only spatial light modulator," 
Opt. Express 24, 6249-6264 (2016)
"""

function generate_hologram(desired, incoming, two_pi_modulation, x_period, y_period, method::Type{T}=BesselJ1) where {T<:HologramMethod}
    dest = similar(desired, UInt8)
    generate_hologram!(dest, desired, incoming, two_pi_modulation, x_period, y_period, method)
    dest
end