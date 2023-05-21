"""
    symetric_integer_range(N,shift=false)

Return a range of integers of length N that is symetric with respect to 0 when N is odd, and has an extra negative number when N is even.

Examples:
    symetric_integer_range(4) = -2:1
    symetric_integer_range(5) = -2:2

If `shift = true`, then we we return an array containing the nonegative and then the negative values. This is equivalent to calling `ifftshift( symetric_integer_range(N,false))`.
"""
symetric_integer_range(N,shift=false) = shift ? vcat(0:N - N÷2 - 1,-N÷2:-1) : (-N÷2:N - N÷2 - 1)

interval(grid) = CUDA.@allowscalar grid[2] - grid[1]

reciprocal_interval(grid) = CUDA.@allowscalar convert(eltype(grid),2π) / length(grid) / interval(grid)

function direct_grid(max,N,shift,use_cuda)
    if use_cuda
        (max / (N ÷ 2)) * symetric_integer_range(length(grid),shift) |> CuArray
    else
        (max / (N ÷ 2)) * symetric_integer_range(length(grid),shift)
    end
end
direct_grid(grid,shift=false) = interval(grid) * symetric_integer_range(length(grid),shift)
direct_grid(grid::CuArray,shift=false) = interval(grid) * symetric_integer_range(length(grid),shift) |> CuArray

reciprocal_grid(grid, shift = false) = reciprocal_interval(grid) * symetric_integer_range(length(grid),shift)
reciprocal_grid(grid::CuArray, shift = false) = reciprocal_interval(grid) * symetric_integer_range(length(grid),shift) |> CuArray

function my_fftshift(x::CuArray)
    fftshift(x)
end

function my_fftshift(x)
    fftshift_view(x)
end

function my_fftshift(x::CuArray,dims)
    fftshift(x,dims)
end

function my_fftshift(x,dims)
    fftshift_view(x,dims)
end

function my_ifftshift(x::CuArray)
    ifftshift(x)
end

function my_ifftshift(x)
    ifftshift_view(x)
end

function my_ifftshift(x::CuArray,dims)
    ifftshift(x,dims)
end

function my_ifftshift(x,dims)
    ifftshift_view(x,dims)
end