using StructuredLight, LinearAlgebra, CUDA
using Test, Documenter

DocMeta.setdocmeta!(StructuredLight, :DocTestSetup, :(using StructuredLight); recursive=true)
doctest(StructuredLight)

@testset "Hermite Gauss" begin
    xs = LinRange(-10, 10, 1024)
    ys = LinRange(-10, 10, 512)

    for m in 0:3, n in 0:3
        ψ0 = hg(xs, ys, m=m, n=n)
        for z in (0.1, 0.5, 1)
            ψ1 = hg(xs, ys, z, m=m, n=n)
            ψ2 = free_propagation(ψ0, xs, ys, z)
            @test isapprox(overlap(ψ1, ψ2, xs, ys), 1, atol=0.03)
        end
    end

end

@testset "Diagonal Hermite Gauss" begin
    xs = LinRange(-10, 10, 1024)
    ys = LinRange(-10, 10, 512)

    for m in 0:3, n in 0:3, z in (0.1, 0.5, 1)
        ψ1 = diagonal_hg(xs, ys, z, m=m, n=n)
        ψ2 = free_propagation(diagonal_hg(xs, ys, m=m, n=n), xs, ys, z)
        @test isapprox(overlap(ψ1, ψ2, xs, ys), 1, atol=0.03)
    end

end

@testset "Laguerre Gauss" begin
    xs = LinRange(-10, 10, 1024)
    ys = LinRange(-10, 10, 512)

    for p in 0:3, l in 0:3, z in (0.1, 0.5, 1)
        ψ1 = lg(xs, ys, z, p=p, l=l)
        ψ2 = free_propagation(lg(xs, ys, p=p, l=l), xs, ys, z)
        @test isapprox(overlap(ψ1, ψ2, xs, ys), 1, atol=0.03)
    end

end

@testset "Scalled Hermite Gauss" begin
    xs = LinRange(-10, 10, 1024)
    ys = LinRange(-10, 10, 512)

    for m in 0:3, n in 0:3, z in (0.2, 1, 2)
        ψ1 = hg(2xs, 2ys, z, m=m, n=n)
        ψ2 = free_propagation(hg(xs, ys, m=m, n=n), xs, ys, z, 2)
        @test isapprox(overlap(ψ1, ψ2, xs, ys), 1 / 4, atol=0.03)
    end

end

@testset "Scalled Diagonal Hermite Gauss" begin
    xs = LinRange(-10, 10, 1024)
    ys = LinRange(-10, 10, 512)

    for m in 0:3, n in 0:3, z in (0.2, 1, 2)
        ψ1 = diagonal_hg(2xs, 2ys, z, m=m, n=n)
        ψ2 = free_propagation(diagonal_hg(xs, ys, m=m, n=n), xs, ys, z, 2)
        @test isapprox(overlap(ψ1, ψ2, xs, ys), 1 / 4, atol=0.03)
    end
end

@testset "Scalled Laguerre Gauss" begin
    xs = LinRange(-10, 10, 1024)
    ys = LinRange(-10, 10, 512)

    for p in 0:3, l in 0:3, z in (0.2, 1, 2)
        ψ1 = lg(2xs, 2ys, z, p=p, l=l)
        ψ2 = free_propagation(lg(xs, ys, p=p, l=l), xs, ys, z, 2)
        @test isapprox(overlap(ψ1, ψ2, xs, ys), 1 / 4, atol=0.03)
    end

end

if CUDA.functional()
    CUDA.allowscalar(false)

    @testset "Hermite Gauss (CUDA)" begin
        xs = LinRange(-10, 10, 1024)
        ys = LinRange(-10, 10, 512)

        for m in 0:3, n in 0:3, z in (0.1, 0.5, 1)
            ψ1 = hg(xs, ys, z, m=m, n=n) |> CuArray
            ψ2 = free_propagation(hg(xs, ys, m=m, n=n) |> CuArray, xs, ys, z)
            @test isapprox(overlap(ψ1, ψ2, xs, ys), 1, atol=0.03)
        end

    end

    @testset "Diagonal Hermite Gauss (CUDA)" begin
        xs = LinRange(-10, 10, 1024)
        ys = LinRange(-10, 10, 512)

        for m in 0:3, n in 0:3, z in (0.1, 0.5, 1)
            ψ1 = diagonal_hg(xs, ys, z, m=m, n=n) |> CuArray
            ψ2 = free_propagation(diagonal_hg(xs, ys, m=m, n=n) |> CuArray, xs, ys, z)
            @test isapprox(overlap(ψ1, ψ2, xs, ys), 1, atol=0.03)
        end

    end

    @testset "Scalled Laguerre Gauss (CUDA)" begin
        xs = LinRange(-10, 10, 1024)
        ys = LinRange(-10, 10, 512)

        for p in 0:3, l in 0:3, z in (0.2, 1, 2)
            ψ1 = lg(2xs, 2ys, z, p=p, l=l) |> CuArray
            ψ2 = free_propagation(lg(xs, ys, p=p, l=l) |> CuArray, xs, ys, z, 2)
            @test isapprox(overlap(ψ1, ψ2, xs, ys), 1 / 4, atol=0.03)
        end

    end

    @testset "Scalled Hermite Gauss (CUDA)" begin
        xs = LinRange(-10, 10, 1024)
        ys = LinRange(-10, 10, 512)

        for m in 0:3, n in 0:3, z in (0.2, 1, 2)
            ψ1 = hg(2xs, 2ys, z, m=m, n=n) |> CuArray
            ψ2 = free_propagation(hg(xs, ys, m=m, n=n) |> CuArray, xs, ys, z, 2)
            @test isapprox(overlap(ψ1, ψ2, xs, ys), 1 / 4, atol=0.03)
        end

    end

    @testset "Scalled Diagonal Hermite Gauss (CUDA)" begin
        xs = LinRange(-10, 10, 1024)
        ys = LinRange(-10, 10, 512)

        for m in 0:3, n in 0:3, z in (0.2, 1, 2)
            ψ1 = diagonal_hg(2xs, 2ys, z, m=m, n=n) |> CuArray
            ψ2 = free_propagation(diagonal_hg(xs, ys, m=m, n=n) |> CuArray, xs, ys, z, 2)
            @test isapprox(overlap(ψ1, ψ2, xs, ys), 1 / 4, atol=0.03)
        end
    end

    @testset "Laguerre Gauss (CUDA)" begin
        xs = LinRange(-10, 10, 1024)
        ys = LinRange(-10, 10, 512)

        for p in 0:3, l in 0:3, z in (0.1, 0.5, 1)
            ψ1 = lg(xs, ys, z, p=p, l=l) |> CuArray
            ψ2 = free_propagation(lg(xs, ys, p=p, l=l) |> CuArray, xs, ys, z)
            @test isapprox(overlap(ψ1, ψ2, xs, ys), 1, atol=0.03)
        end
    end
end