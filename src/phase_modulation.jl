@kernel function zernike_kernel!(dest, x, y, m, n)
    i, j = @index(Global, NTuple)
    dest[i, j] = zernike_polynomial(x[i], y[j], m, n)
end

"""
    zernike_polynomial!(dest, x, y, m, n)

Same as [`zernike_polynomial`](@ref), but writes the result to `dest`.
"""
function zernike_polynomial!(dest, x, y, m, n)
    backend = get_backend(dest)
    kernel! = zernike_kernel!(backend)
    kernel!(dest, x, y, m, n; ndrange=size(dest))
end

"""
    zernike_polynomial(x, y, m, n; backend=CPU())

    Evaluate the `n` th Zernike polynomial at `(x, y)` with azimuthal order `m`.
"""
function zernike_polynomial(x, y, m, n; backend=CPU())
    dest = allocate(backend, promote_type(eltype(x), eltype(y)), length(x), length(y))
    zernike_polynomial!(dest, x, y, m, n)
    dest
end

function lens(x::Number, y::Number, fx, fy; k=1)
    cis(-k * (x^2 / fx + y^2 / fy) / 2)
end

@kernel function lens_kernel!(dest, x, y, fx, fy, k)
    i, j = @index(Global, NTuple)
    dest[i, j] = lens(x[i], y[j], fx, fy; k)
end

"""
    lens!(dest, x, y, fx, fy; k=1)

Same as [`lens`](@ref), but writes the result to `dest`.
"""
function lens!(dest, x, y, fx, fy; k=1)
    backend = get_backend(dest)
    kernel! = lens_kernel!(backend)
    kernel!(dest, x, y, fx, fy, k; ndrange=size(dest))
end

@doc raw"""
    lens(x, y, fx, fy; k=1, backend=CPU())

Output an array containing the phase shift introduced by a lens of focal lengths `fx` and `fy`:

```math
\exp\left[-i k \left( \frac{x^2}{2 f_x} + \frac{y^2}{2 f_y} \right)\right]
```

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* lens(x,y,fx,fy;k=k)`
"""
function lens(x, y, fx, fy; k=1, backend=CPU())
    dest = allocate(backend, complex_type(x, y, fx, fy, k), length(x), length(y))
    lens!(dest, x, y, fx, fy; k)
    dest
end

"""
    tilted_lens!(dest, x, y, f, ϕ; k=1)

Same as [`tilted_lens`](@ref), but writes the result to `dest`.
"""
function tilted_lens!(dest, x, y, f, ϕ; k=1)
    lens!(dest, x, y, sec(ϕ) * f, cos(ϕ) * f; k)
end

@doc raw"""
    tilted_lens(x,y,f,ϕ;k=1)

Output an array containing the phase shift introduced by a spherical lens of focal length `f` tilted by an angle `ϕ`.

This is a special case of the lens function, where the focal lengths are adjusted based on the tilt angle: ``f_x = \sec(ϕ) * f`` and ``f_y = \cos(ϕ) * f``.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* tilted_lens(x,y,f,ϕ;k=k)`
"""
function tilted_lens(x, y, f, ϕ; k=1, backend=CPU())
    lens(x, y, sec(ϕ) * f, cos(ϕ) * f; k, backend)
end

@kernel function apply_phase_kernel!(dest, funcs, coeffs, grid)
    J = @index(Global, NTuple)
    dest[J...] *= cis(linear_combination(funcs, coeffs, ntuple(j -> grid[j][J[j]], length(grid))))
end

"""
    apply_phase!(dest, funcs, coeffs, grid)

Apply a phase modulation defined by `funcs` and `coeffs` to the array `dest`.

The phase modulation consists of multiplying `dest` by `exp(i * φ)`, where `φ` is a linear combination of functions `funcs` with coefficients `coeffs` evaluated at coordinates defined by `grid`.

`grid` is a tuple of collections, each representing a coordinate axis.

This function has many applications including aberration correction, lens simulation, phase masks, and general wavefront shaping. The `funcs` can include [`zernike_polynomial`](@ref), [`lens`](@ref), or any custom phase functions.

```jldoctest
rs = LinRange(-3, 3, 3)
u = hg(rs, rs)

f1(args) = zernike_polynomial(args..., 1, 1)
f2(args) = zernike_polynomial(args..., 2, 2)

coeffs = (0.1, 0.2)

apply_phase!(u, (f1, f2), coeffs, (rs, rs))

u

# output

3×3 Matrix{ComplexF64}:
   1.1609e-8-3.59109e-9im   6.96526e-6+9.82201e-5im    1.1609e-8-3.59109e-9im
 -2.23719e-5-9.58916e-5im     0.797885+0.0im         -2.23719e-5-9.58916e-5im
   1.1609e-8+3.59109e-9im  -4.97106e-5+8.49974e-5im    1.1609e-8+3.59109e-9im
```
"""
function apply_phase!(dest, funcs, coeffs, grid)
    backend = get_backend(dest)
    kernel! = apply_phase_kernel!(backend)
    kernel!(dest, funcs, coeffs, grid; ndrange=size(dest))
end