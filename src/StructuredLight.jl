module StructuredLight

using Reexport

using FourierTools

include("dft_utils.jl")

using Images, VideoIO
@reexport using ColorSchemes
using Tullio, CUDAKernels, KernelAbstractions
import SpecialFunctions: beta 

include("initial_profiles.jl")
export lg,hg,diagonal_hg,lens,tilted_lens

include("visualization.jl")
export visualize,show_animation,save_animation

include("propagation_functions.jl")

include("free_propagation.jl")
export free_propagation

include("kerr_propagation.jl")
export kerr_propagation,kerr_propagation_old

include("misc.jl")
export overlap

end
