module StructuredLight

using Reexport

using FFTW, CUDA, CUDA.CUFFT

include("dft_utils.jl")
export symetric_integer_range,interval,reciprocal_interval,direct_grid,reciprocal_grid

include("steps.jl")
export dispersion_step!,get_dispersion_phases,evolve_phase

using Images, VideoIO
@reexport using ColorSchemes
using ThreadsX
using SpecialFunctions

include("initial_profiles.jl")
export lg,hg,diagonal_hg,lens,tilted_lens

include("visualization.jl")
export visualize,show_animation,save_animation

include("free_propagation.jl")
export free_propagation

include("kerr_propagation.jl")
export kerr_propagation

end
