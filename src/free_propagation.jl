to_device(::Any, y) = y

function create_shifted_ψ(ψ, ::Number)
    shifted_ψ = similar(ψ, complex(eltype(ψ)))
    ifftshift!(shifted_ψ, ψ)
end

function create_shifted_ψ(ψ, z)
    shifted_ψ = similar(ψ, complex(eltype(ψ)), (size(ψ)..., length(z)))
    for dest ∈ eachslice(shifted_ψ, dims=3)
        ifftshift!(dest, ψ)
    end
    shifted_ψ
end

get_ndrange(ψ, ::Number) = size(ψ)
get_ndrange(ψ, z) = (size(ψ)..., length(z))

@kernel function fresnel_kernel!(ψ, x, y, z::Number, k)
    i, j = @index(Global, NTuple)
    ψ[i, j] *= cis(-z * (x[j]^2 + y[i]^2) / 2k)
end

@kernel function fresnel_kernel!(ψ, x, y, z::Number, k, scaling)
    i, j = @index(Global, NTuple)
    ψ[i, j] *= cis(-z / scaling * (x[j]^2 + y[i]^2) / 2k)
end

@kernel function pre_kernel!(ψ, x, y, z::Number, k, scaling)
    i, j = @index(Global, NTuple)
    ψ[i, j] *= cis(k * (x[j]^2 + y[i]^2) * (1 - scaling) / 2z) / scaling
end

@kernel function post_kernel!(ψ, x, y, z::Number, k, scaling)
    i, j = @index(Global, NTuple)
    ψ[i, j] *= cis(k * (x[j]^2 + y[i]^2) * (scaling - 1) * scaling / 2z)
end

@kernel function fresnel_kernel!(ψ, x, y, z, k)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(-z[l] * (x[j]^2 + y[i]^2) / 2k)
end

@kernel function fresnel_kernel!(ψ, x, y, z, k, scaling)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(-z[l] / scaling[l] * (x[j]^2 + y[i]^2) / 2k)
end

@kernel function pre_kernel!(ψ, x, y, z, k, scaling)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(k * (x[j]^2 + y[i]^2) * (1 - scaling[l]) / 2z[l]) / scaling[l]
end

@kernel function post_kernel!(ψ, x, y, z, k, scaling)
    i, j, l = @index(Global, NTuple)
    ψ[i, j, l] *= cis(k * (x[j]^2 + y[i]^2) * (scaling[l] - 1) * scaling[l] / 2z[l])
end


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
    qx = to_device(ψ, reciprocal_grid(x, shift=true))
    qy = to_device(ψ, reciprocal_grid(y, shift=true))

    shifted_ψ = create_shifted_ψ(ψ, z)

    backend = get_backend(ψ)
    _fresnel_kernel! = fresnel_kernel!(backend, 256)
    ndrange = get_ndrange(ψ, z)

    fft!(shifted_ψ, (1, 2))
    _fresnel_kernel!(shifted_ψ, qx, qy, z, k; ndrange)
    ifft!(shifted_ψ, (1, 2))

    fftshift(shifted_ψ, (1, 2))
end

function free_propagation(ψ, x, y, z, scaling; k=1)
    @assert length(z) == length(scaling) "`z` and `scaling` should have the same length"
    @assert 0 ∉ z "This method does not support `z` containing `0`"

    _x = to_device(ψ, ifftshift(x))
    _y = to_device(ψ, ifftshift(y))
    _z = to_device(ψ, z)
    qx = to_device(ψ, reciprocal_grid(x, shift=true))
    qy = to_device(ψ, reciprocal_grid(y, shift=true))

    shifted_ψ = create_shifted_ψ(ψ, z)

    backend = get_backend(ψ)
    _fresnel_kernel! = fresnel_kernel!(backend, 256)
    _pre_kernel! = pre_kernel!(backend, 256)
    _post_kernel! = post_kernel!(backend, 256)
    ndrange = get_ndrange(ψ, z)

    _pre_kernel!(shifted_ψ, _x, _y, _z, k, scaling; ndrange)
    fft!(shifted_ψ, (1, 2))
    _fresnel_kernel!(shifted_ψ, qx, qy, z ./ scaling, k; ndrange)
    ifft!(shifted_ψ, (1, 2))
    _post_kernel!(shifted_ψ, _x, _y, _z, k, scaling; ndrange)

    fftshift(shifted_ψ, (1, 2))
end