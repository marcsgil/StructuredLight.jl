module StructuredLight

using FFTW, LinearAlgebra, KernelAbstractions, Bessels

FFTW.set_num_threads(8)

function visualize end
function save_animation end
export visualize, save_animation

include("dft_utils.jl")
export direct_grid, reciprocal_grid

include("ortho_poly.jl")
export laguerre, hermite

include("initial_profiles.jl")
export hg!, hg, diagonal_hg, diagonal_hg!, lg!, lg

include("free_propagation.jl")
export free_propagation

include("kerr_propagation.jl")
export kerr_propagation

include("misc.jl")
export overlap

include("holograms.jl")
export BesselJ!, generate_hologram, generate_hologram!

end