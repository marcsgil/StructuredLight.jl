using Preferences: set_preferences!
set_preferences!("StructuredLight", "dispatch_doctor_mode" => "error")

using StructuredLight, LinearAlgebra, CUDA
using Test, Documenter

DocMeta.setdocmeta!(StructuredLight, :DocTestSetup, :(using StructuredLight); recursive=true)

doctest(StructuredLight)

include("propagation_test.jl")

include("kerr_test.jl")

include("hologram_tests.jl")

include("normalization_tests.jl")

include("shapes_test.jl")

include("aberration_correction_tests.jl")