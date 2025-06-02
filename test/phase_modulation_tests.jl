if CUDA.functional()
    @testset "Phase Modulation (CUDA)" begin
        rs = LinRange(-3, 3, 3)
        u = hg(rs, rs)

        f1(args) = zernike_polynomial(args..., 1, 1)
        f2(args) = zernike_polynomial(args..., 2, 2)

        coeffs = (0.1, 0.2)

        apply_phase!(u, (f1, f2), coeffs, (rs, rs))

        u_cuda = hg(rs, rs, backend=CUDABackend())
        apply_phase!(u_cuda, (f1, f2), coeffs, (rs, rs))
        @test Array(u_cuda) ≈ Array(u)
    end
end