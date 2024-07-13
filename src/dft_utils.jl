"""
    symetric_integer_range(N,shift=false)

Return a range of integers of length N that is symetric with respect to 0 when N is odd, and has an extra negative number when N is even.

Examples:
    symetric_integer_range(4) = -2:1
    symetric_integer_range(5) = -2:2

If `shift = true`, then we we return an array containing the nonegative and then the negative values. This is equivalent to calling `ifftshift( symetric_integer_range(N,false))`.
"""
symetric_integer_range(N; shift=false) = shift ? vcat(0:N-N÷2-1, -N÷2:-1) : (-N÷2:N-N÷2-1)

interval(grid) = grid[2] - grid[1]

reciprocal_interval(grid) = convert(eltype(grid), 2π) / length(grid) / interval(grid)

function direct_grid(max, N; shift=false)
    (max / (N ÷ 2)) * symetric_integer_range(N; shift)
end

direct_grid(grid; shift=false) = interval(grid) * symetric_integer_range(length(grid); shift)

reciprocal_grid(grid; shift=false) = reciprocal_interval(grid) * symetric_integer_range(length(grid); shift)