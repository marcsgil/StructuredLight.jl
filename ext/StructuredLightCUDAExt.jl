module StructuredLightCUDAExt

using StructuredLight, CUDA, CUDA.CUFFT
using Tullio, KernelAbstractions

include("free_propagation.jl")
include("kerr_propagation.jl")
include("visualization.jl")

end