module StructuredLight

using FourierTools, LinearAlgebra, Tullio

include("dft_utils.jl")
export direct_grid, reciprocal_grid

include("ortho_poly.jl")
export laguerre, hermite

include("initial_profiles.jl")
export _hg, rotated_hg, hg, diagonal_hg, _lg, lg, lens, tilted_lens

include("free_propagation.jl")
export free_propagation

include("kerr_propagation.jl")
export kerr_propagation

include("misc.jl")
export overlap

end