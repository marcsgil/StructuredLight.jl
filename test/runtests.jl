using StructuredLight, LinearAlgebra, CUDA
using Test, Documenter
DocMeta.setdocmeta!(StructuredLight, :DocTestSetup, :(using StructuredLight); recursive=true)

doctest(StructuredLight)

include("propagation_test.jl")

include("kerr_test.jl")