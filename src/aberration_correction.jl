@kernel function zernike_kernel!(dest, x, y, m, n)
    i, j = @index(Global, NTuple)
    dest[i, j] = zernike_polynomial(x[i], y[j], m, n)
end

function zernike_polynomial!(dest, x, y, m, n)
    backend = get_backend(dest)
    kernel! = zernike_kernel!(backend)
    kernel!(dest, x, y, m, n; ndrange=size(dest))
end

function zernike_polynomial(x, y, m, n)
    dest = similar(x, promote_type(eltype(x), eltype(y)), size(x, 1), size(y, 1))
    zernike_polynomial!(dest, x, y, m, n)
    dest
end

@kernel function lens_kernel!(dest, x, y, fx, fy, k)
    i, j = @index(Global, NTuple)
    dest[i, j] *= cis(-k * (x[i]^2 / fx + y[j]^2 / fy) / 2)
end

function lens!(dest, x, y, fx, fy; k=1)
    backend = get_backend(dest)
    kernel! = lens_kernel!(backend)
    kernel!(dest, x, y, fx, fy, k; ndrange=size(dest))
end

"""
    lens(x,y,fx,fy;k=1)

Output an array containing the phase shift introduced by a lens of focal lengths `fx` and `fy`.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* lens(x,y,fx,fy;k=k)`
"""
function lens(x, y, fx, fy; k=1)
    dest = fill(one(complex_type(x, y, fx, fy, k)), length(x), length(y))
    lens!(dest, x, y, fx, fy; k)
    dest
end

function tilted_lens!(dest, x, y, f, ϕ, ; k=1)
    lens!(dest, x, y, sec(ϕ) * f, cos(ϕ) * f; k)
end

"""
    tilted_lens(x,y,f,ϕ;k=1)

Output an array containing the phase shift introduced by a spherical lens of focal length `f` tilted by an angle `ϕ`.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* tilted_lens(x,y,f,ϕ;k=k)`
"""
function tilted_lens(x, y, f, ϕ; k=1)
    lens(x, y, sec(ϕ) * f, cos(ϕ) * f; k)
end

@kernel function aberration_correction_kernel!(dest, funcs, coeffs, grid)
    J = @index(Global, NTuple)
    dest[J...] *= cis(linear_combination(funcs, coeffs, ntuple(j -> grid[j][J[j]], length(grid))))
end

"""
    aberration_correction!(dest, funcs, coeffs, grid)

Apply the aberration correction defined by `funcs` and `coeffs` to the array `dest`.

The correction consists of applying a phase defined by linear combination of functions `funcs` with coefficients `coeffs` to the values at the coordinates defined by `grid`.

`grid` is a tuple of collections, each representing a coordinate axis.

This is used to correct for aberrations in optical systems, where `funcs` can include [`zernike_polynomials`](@ref) or other phase functions.

```jldoctest
rs = LinRange(-3, 3, 3)
u = hg(rs, rs)

f1(args) = zernike_polynomial(args..., 1, 1)
f2(args) = zernike_polynomial(args..., 2, 2)

coeffs = (0.1, 0.2)

aberration_correction!(u, (f1, f2), coeffs, (rs, rs))

u

# output

3×3 Matrix{ComplexF64}:
   1.1609e-8-3.59109e-9im   6.96526e-6+9.82201e-5im    1.1609e-8-3.59109e-9im
 -2.23719e-5-9.58916e-5im     0.797885+0.0im         -2.23719e-5-9.58916e-5im
   1.1609e-8+3.59109e-9im  -4.97106e-5+8.49974e-5im    1.1609e-8+3.59109e-9im
```
"""
function aberration_correction!(dest, funcs, coeffs, grid)
    backend = get_backend(dest)
    kernel! = aberration_correction_kernel!(backend)
    kernel!(dest, funcs, coeffs, grid; ndrange=size(dest))
end