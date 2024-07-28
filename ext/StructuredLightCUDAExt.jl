module StructuredLightCUDAExt

using StructuredLight, CUDA

StructuredLight.to_device(::CuArray, y::AbstractArray) = CuArray(y)

using Tullio, KernelAbstractions, CUDA.CUFFT
include("kerr_propagation.jl")

end