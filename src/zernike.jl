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