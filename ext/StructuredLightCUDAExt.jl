module StructuredLightCUDAExt

using StructuredLight, CUDA, CUDA.CUFFT
using Tullio, CUDAKernels, KernelAbstractions

include("free_propagation.jl")
include("kerr_propagation.jl")

    
end