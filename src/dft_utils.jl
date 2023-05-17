function symetric_integer_range(N)
    #=Returns a range of integers of length N that is symetric with respect to 0 when N is odd, and has an extra negative number when N is even.

    Examples:
    symetric_integer_range(4) = -2:1
    symetric_integer_range(5) = -2:2
    =#

    -N÷2:N - N÷2 - 1
end

interval(grid) = grid[2] - grid[1]

reciprocal_interval(grid) = π/last(grid)

direct_grid(grid) = interval(grid)*symetric_integer_range(length(grid))
direct_grid(grids...) = Iterators.product((direct_grid(grid) for grid in grids)...)

reciprocal_grid(grid) = reciprocal_interval(grid)*symetric_integer_range(length(grid))
reciprocal_grid(grids...) = Iterators.product((reciprocal_grid(grid) for grid in grids)...)