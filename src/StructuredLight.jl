module StructuredLight

using Reexport

using FourierTools, LinearAlgebra

using Images, VideoIO
@reexport using ColorSchemes
using Tullio

include("dft_utils.jl")
export direct_grid, reciprocal_grid

include("ortho_poly.jl")
export laguerre, hermite

include("initial_profiles.jl")
export rotated_hg, hg, diagonal_hg, lg, lens, tilted_lens

include("visualization.jl")
export visualize, show_animation, save_animation

include("free_propagation.jl")
export free_propagation

include("kerr_propagation.jl")
export kerr_propagation

include("misc.jl")
export overlap

end