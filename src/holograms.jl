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

    isapprox(f(x₁), 0; atol=tol) && return float(x₁)
    isapprox(f(x₂), 0; atol=tol) && return float(x₂)

    x = (x₁ + x₂) / 2

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

    x
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

@kernel function hologram_kernel!(dest, relative, M, two_pi_modulation, x_period, y_period, ::BesselJ1)
    i, j = @index(Global, NTuple)

    ψ = (inverse_besselj1(x_max_besselj1 * abs(relative[i, j] / M))
         *
         sin(2π * (i / x_period + j / y_period) + angle(relative[i, j])))

    dest[i, j] = round(two_pi_modulation * 0.586 * (ψ / y_max_besselj1 + 1) / 2)
end

@kernel function hologram_kernel!(dest, relative, M, two_pi_modulation, x_period, y_period, ::Simple)
    i, j = @index(Global, NTuple)

    ψ = abs(relative[i, j] / M) * normalize_angle(angle(relative[i, j]) + 2π * (i / x_period + j / y_period))

    dest[i, j] = round(two_pi_modulation * (ψ + π) / 2π)
end

"""
    generate_hologram!(dest, desired, incoming, two_pi_modulation, 
        x_period, y_period, method::Type{T}=BesselJ1) where {T<:HologramMethod}

Same as [`generate_hologram`](@ref), but writes the result to `dest`.
"""
function generate_hologram!(dest, relative, two_pi_modulation, x_period, y_period, method::T=BesselJ1()) where {T<:HologramMethod}
    M = maximum(abs, relative)
    backend = get_backend(dest)
    kernel! = hologram_kernel!(backend)
    kernel!(dest, relative, M, two_pi_modulation, x_period, y_period, method; ndrange=size(dest))
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

```jldoctest
x = -5:5
y = -5:5
relative = lg(x, y, l=1)
holo = generate_hologram(relative, 255, 2, 2)

# output

11×11 Matrix{UInt8}:
 0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4b  0x4c  0x4b  0x4a  0x4b  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4c  0x38  0x4b  0x5d  0x49  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x46  0x95  0x4b  0x00  0x4f  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4c  0x38  0x4b  0x5d  0x49  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4b  0x4c  0x4b  0x4a  0x4b  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b
 0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b  0x4b
```

# References

[1] Victor Arrizón, Ulises Ruiz, Rosibel Carrada, and Luis A. González,
"Pixelated phase computer holograms for the accurate encoding of scalar complex fields,"
J. Opt. Soc. Am. A 24, 3500-3507 (2007)

[2] Thomas W. Clark, Rachel F. Offer, Sonja Franke-Arnold, Aidan S. Arnold, and Neal Radwell, 
"Comparison of beam generation techniques using a phase only spatial light modulator," 
Opt. Express 24, 6249-6264 (2016)
"""
function generate_hologram(relative, two_pi_modulation, x_period, y_period, method::T=BesselJ1()) where {T<:HologramMethod}
    dest = similar(relative, UInt8)
    generate_hologram!(dest, relative, two_pi_modulation, x_period, y_period, method)
    dest
end