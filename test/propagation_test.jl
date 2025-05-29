function propagation_test(mode, xs, ys, zs; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps(eltype(xs)), kwargs...)
    ψ₀ = mode(xs, ys; kwargs...)
    for z ∈ zs
        ψ₁ = mode(xs, ys, z; kwargs...)
        ψ₂ = free_propagation(ψ₀, xs, ys, z)
        @test isapprox(overlap(ψ₁, ψ₂, xs, ys), 1; rtol, atol)
    end
end

function propagation_test(mode, xs, ys, zs, scaling; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps(eltype(xs)), kwargs...)
    ψ₀ = mode(xs, ys; kwargs...)
    for z ∈ 2 .* zs
        ψ₁ = mode(2xs, 2ys, z; kwargs...)
        ψ₂ = free_propagation(ψ₀, xs, ys, z, scaling)
        @test isapprox(overlap(ψ₁, ψ₂, xs, ys), 1 / scaling^2; rtol, atol)
    end
end

function run_propagation_tests(xs, ys, identifier=""; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps(eltype(xs)))
    zs = (0.1, 0.5, 1)
    scaling = 2

    @testset "Hermite Gauss$identifier" begin
        for m ∈ 0:3, n ∈ 0:3
            propagation_test(diagonal_hg, xs, ys, zs; m, n, rtol)
            propagation_test(diagonal_hg, xs, ys, zs, scaling; m, n, rtol)
            for θ ∈ LinRange(0, π, 5)
                propagation_test(hg, xs, ys, zs; m, n, θ, rtol)
                propagation_test(hg, xs, ys, zs, scaling; m, n, rtol)
            end
        end
    end

    @testset "Laguerre Gauss$identifier" begin
        for p ∈ 0:3, l ∈ 0:3
            propagation_test(lg, xs, ys, zs; p, l, rtol)
            propagation_test(lg, xs, ys, zs, scaling; p, l, rtol)
        end
    end
end

xs = LinRange(-20, 20, 1024)
ys = LinRange(-10, 10, 512)
rtol = 0.03

run_propagation_tests(xs, ys; rtol)

#= if CUDA.functional()
    CUDA.allowscalar(false)

    run_propagation_tests(xs |> collect |> CuArray, ys |> collect |> CuArray, " (CUDA)"; rtol=rtol)
end =#