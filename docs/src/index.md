# Introduction

This package is a comprehensive toolkit for working with spatially structured light beams in Julia. It offers a range of functionalities for simulating, propagating, and visualizing structured light beams, such as Laguerre-Gauss and Hermite-Gauss modes. It also implements the computation of holograms used for producing such beams in Spatial Light Modulators (SLMs). The package supports both CPU and GPU computations, allowing for efficient simulations on GPUs.

Here is a brief overview of the main features:
- **[Beam Profiles](@ref)**: Calculate a large variety of beam profiles, such as ([Hermite-Gauss](@ref)) and [Laguerre-Gauss](@ref) modes.
- **[Visualization](@ref)**: Tools for visualizing the beam profiles and their propagation.
- **[Propagation](@ref)**: Simulate the propagation of beams in free space and in Kerr media.
- **[Holograms](@ref)**: Compute holograms for generating structured light beams on SLMs.
- **[Phase Modulation](phase_modulation.md)**: Tools for applying phase modulation to structured light beams, including aberration correction, lens simulation, and custom phase masks using Zernike polynomials.
- **[GPU Support](@ref)**: Efficient computation on Nvidia GPUs using CUDA.

Check [Quick Start](@ref) for a quick introduction to the package.