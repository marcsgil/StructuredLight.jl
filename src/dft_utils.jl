"""
    symetric_integer_range(N,shift=false)

Return a range of integers of length N that is symetric with respect to 0 when N is odd, and has an extra negative number when N is even.

Examples:
    symetric_integer_range(4) = -2:1
    symetric_integer_range(5) = -2:2

If `shift = true`, then we we return an array containing the nonegative and then the negative values. This is equivalent to calling `ifftshift( symetric_integer_range(N,false))`.
"""
symetric_integer_range(N,shift=false) = shift ? vcat(0:N - N÷2 - 1,-N÷2:-1) : (-N÷2:N - N÷2 - 1)

struct DFTGrid{use_gpu}
    max::AbstractFloat
    N::Integer
end

DFTGrid(max,N,use_gpu=false) = DFTGrid{use_gpu}(max,N)
DFTGrid(grid::AbstractArray) = DFTGrid(-minimum(grid),length(grid),grid isa CuArray)

interval(grid::DFTGrid) = grid.max / (grid.N ÷ 2)
interval(grid) = CUDA.@allowscalar grid[2] - grid[1]
reciprocal_interval(grid::DFTGrid) = oftype(grid.max,2π) / grid.N / interval(grid)

direct_grid(grid::DFTGrid{false},shift=false) = interval(grid) .* symetric_integer_range(grid.N,shift)
direct_grid(grid::DFTGrid{true},shift=false) = interval(grid) .* symetric_integer_range(grid.N,shift) |> CuArray
direct_grid(grid::AbstractArray,shift=false) = direct_grid(DFTGrid(grid,false),shift)
direct_grid(grid::CuArray,shift=false) = direct_grid(DFTGrid(grid,true),shift)

reciprocal_grid(grid::DFTGrid{false},shift=false) = reciprocal_interval(grid) .* symetric_integer_range(grid.N,shift)
reciprocal_grid(grid::DFTGrid{true},shift=false) = reciprocal_interval(grid) .* symetric_integer_range(grid.N,shift) |> CuArray
reciprocal_grid(grid::AbstractArray,shift=false) = reciprocal_grid(DFTGrid(grid,false),shift)
reciprocal_grid(grid::CuArray,shift=false) = reciprocal_grid(DFTGrid(grid,true),shift)