module StructuredLightCUDAExt

using StructuredLight, CUDA

StructuredLight.to_device(::CuArray, y::AbstractArray) = CuArray(y)

end