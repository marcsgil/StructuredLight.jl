module StructuredLight

using FFTW, LinearAlgebra, KernelAbstractions

FFTW.set_num_threads(8)

function visualize end
function save_animation end
export visualize, save_animation

include("utils.jl")

include("dft_utils.jl")
export direct_grid, reciprocal_grid

include("ortho_poly.jl")
export laguerre, hermite

include("initial_profiles.jl")
export _hg, rotated_hg, hg, diagonal_hg, _lg, lg, lens, tilted_lens
export lg!

include("free_propagation.jl")
export free_propagation

include("kerr_propagation.jl")
export kerr_propagation

include("misc.jl")
export overlap

end