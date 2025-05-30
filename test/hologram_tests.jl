if CUDA.functional()
    @testset "Holograms (CUDA)" begin
        x = -5:5
        y = -5:5
        relative = lg(x, y, l=1)
        relative_cuda = CuArray(relative)

        for method in (BesselJ1(), Simple())
            holo = generate_hologram(relative, 255, 2, 2, method)
            holo_cuda = generate_hologram(relative_cuda, 255, 2, 2, method)
            @test holo ≈ Array(holo_cuda)
        end
    end
end