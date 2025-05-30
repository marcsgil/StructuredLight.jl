module StructuredLightCUDAExt
using CUDA
import StructuredLight: xs_besselj1, ys_besselj1, select_itp, BesselJ1, HologramMethod, generate_hologram!, hologram_kernel!

select_itp(::CuArray, ::BesselJ1) = (xs_besselj1, CuArray(ys_besselj1))

end