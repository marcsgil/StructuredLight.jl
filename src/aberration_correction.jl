function linear_combination(fs, coeffs)
    (args...; kwargs...) -> sum(f(args...; kwargs) * c for (f, c) in zip(fs, coeffs))
end

@kernel function aberration_correction_kernel!(dest, correction_func, x, y)
    i, j = @index(Gloabal, NTuple)
    dest[i, j] *= cis(correction_func(x[i], y[j]))
end

function aberration_correction!(dest, x, y, scale, idxs, coeffs)
    fs = ((x, y) -> zernike_polynomial(x / scale, y / scale, idxs...) for idx ∈ idxs)
    correction_func = linear_combination(fs, coeffs)


end

"""
    lens(x,y,fx,fy;k=1)

Output an array containing the phase shift introduced by a lens of focal lengths `fx` and `fy`.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* lens(x,y,fx,fy;k=k)`
"""
function lens(x, y, fx, fy; k=1)
    @tullio result[i, j] := cis(-k * (x[i]^2 / fx + y[j]^2 / fy) / 2)
end

"""
    tilted_lens(x,y,f,ϕ;k=1)

Output an array containing the phase shift introduced by a spherical lens of focal length `f` tilted by an angle `ϕ`.

The calculation is done over a grid defined by `x` and `y`.

`k` is the incident wavenumber.

To apply the lens at a beam `ψ₀`, just calculate `ψ = ψ₀ .* tilted_lens(x,y,f,ϕ;k=k)`
"""
function tilted_lens(x, y, f, ϕ; k=1)
    lens(x, y, sec(ϕ) * f, cos(ϕ) * f, k=k)
end