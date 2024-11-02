# TODO
function linear_combination(fs, coeffs)
    (args...; kwargs...) -> sum(f(args...; kwargs) * c for (f, c) in zip(fs, coeffs))
end

@kernel function aberration_correction_kernel!(dest, correction_func, x, y)
    i, j = @index(Global, NTuple)
    dest[i, j] *= cis(correction_func(x[i], y[j]))
end

function aberration_correction!(dest, x, y, scale, idxs, coeffs)
    fs = ((x, y) -> zernike_polynomial(x / scale, y / scale, idxs...) for idx ∈ idxs)
    correction_func = linear_combination(fs, coeffs)
end
# End TODO

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