function get_shifted_ψ(ψ, x, y, z)
    shifted_ψ = similar(ψ, complex(eltype(ψ)), get_size(x, y, z)...)
    for slice ∈ eachslice(shifted_ψ, dims=3)
        fftshift!(slice, ψ, (1, 2))
    end
    shifted_ψ
end

function get_shifted_ψ(ψ, x, y, z::Number)
    shifted_ψ = similar(ψ, complex(eltype(ψ)), get_size(x, y, z)...)
    fftshift!(shifted_ψ, ψ)
    shifted_ψ
end

#Without scaling
@kernel function fresnel_kernel!(ψ, x, y, z, k)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(-z[l] * (x[i]^2 + y[j]^2) / 2k)
end

#With scaling
@kernel function pre_kernel!(ψ, x, y, z, k, scaling)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(k * (x[i]^2 + y[j]^2) * (1 - scaling[l]) / 2z[l]) / scaling[l]
end

@kernel function fresnel_kernel!(ψ, x, y, z, k, scaling)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(-z[l] / scaling[l] * (x[i]^2 + y[j]^2) / 2k)
end

@kernel function post_kernel!(ψ, x, y, z, k, scaling)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(k * (x[i]^2 + y[j]^2) * (scaling[l] - 1) * scaling[l] / 2z[l])
end

three_d_size(x::AbstractMatrix) = (size(x)..., 1)
three_d_size(x::AbstractArray{T,3}) where {T} = size(x)


"""
    free_propagation(ψ, x, y, z [, scaling]; k=1)
Propagate an inital profile `ψ`.

The propagation is the solution of `∇² ψ + 2ik ∂_z ψ = 0` at distance `z` under the initial condition `ψ`.

`x` and `y` are the grids over which `ψ` is calculated.

If `z` is an `AbstractArray`, the output is a 3D array representing the solution at every element of `z`.

The output at a distance `z[n]` is calculated on a scalled grid defined by `scaling[n] * x` and `scaling[n] * y`.

`k` is the wavenumber.

# Example

```jldoctest
x = LinRange(-10, 10, 256)
y = LinRange(-10, 10, 512)
z = LinRange(0.1, 1, 10)

ψ = hg(x, y; m=3, n=2)
ψ′ = hg(2x, 2y; m=3, n=2)

(
    free_propagation(ψ, x, y, z) ≈ stack(free_propagation(ψ, x, y, z) for z ∈ z)
    &&
    free_propagation(ψ, x, y, z, fill(2, length(z))) ≈ stack(free_propagation(ψ, x, y, z, 2) for z ∈ z)
    &&
    free_propagation(ψ, x, y, 0.5, 2) ≈ free_propagation(ψ′, 2x, 2y, 0.5)
)

# output
true
```
"""
function free_propagation(ψ, x, y, z; k=1)
    qx = fftfreq(length(x), 2π / step(x))
    qy = fftfreq(length(y), 2π / step(y))
    result = stack(ψ for _ in z)

    backend = get_backend(result)
    _fresnel_kernel! = fresnel_kernel!(backend)
    ndrange = three_d_size(result)

    fft!(result, (1, 2))
    _fresnel_kernel!(result, qx, qy, z, k; ndrange)
    ifft!(result, (1, 2))

    result
end

function free_propagation(ψ, x, y, z, scaling; k=1)
    @assert length(z) == length(scaling) "`z` and `scaling` should have the same length"
    @assert 0 ∉ z "This method does not support `z` containing `0`"

    _x = to_device(ψ, ifftshift(x))
    _y = to_device(ψ, ifftshift(y))
    _z = to_device(ψ, z)
    qx = fftfreq(length(x), 2π / step(x))
    qy = fftfreq(length(y), 2π / step(y))

    shifted_ψ = get_shifted_ψ(ψ, x, y, z)

    shifted_ψ_3d = reshape(shifted_ψ, length(x), length(y), length(z))

    backend = get_backend(shifted_ψ_3d)
    _fresnel_kernel! = fresnel_kernel!(backend)
    _pre_kernel! = pre_kernel!(backend)
    _post_kernel! = post_kernel!(backend)
    ndrange = size(shifted_ψ_3d)

    _pre_kernel!(shifted_ψ_3d, _x, _y, _z, k, scaling; ndrange)
    fft!(shifted_ψ_3d, (1, 2))
    _fresnel_kernel!(shifted_ψ_3d, qx, qy, z ./ scaling, k; ndrange)
    ifft!(shifted_ψ_3d, (1, 2))
    _post_kernel!(shifted_ψ_3d, _x, _y, _z, k, scaling; ndrange)

    fftshift(shifted_ψ, (1, 2))
end