function test_overlaps(xs, ys, device="")
    @testset "Normalization HG $device" begin
        for m ∈ 0:10, n ∈ 0:10
            ψ = hg(xs, ys; m, n)
            @test overlap(ψ, ψ, xs, ys) ≈ 1
        end
    end

    @testset "Normalization LG $device" begin
        for p ∈ 0:5, l ∈ 0:10
            ψ = lg(xs, ys; p, l)
            @test overlap(ψ, ψ, xs, ys) ≈ 1
        end
    end
end

rs = LinRange(-10, 10, 1024)

test_overlaps(rs, rs)

if CUDA.functional()
    test_overlaps(rs |> collect |> cu, rs |> collect |> cu, "(CUDA)")
end