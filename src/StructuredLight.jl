module StructuredLight

using FFTW, LinearAlgebra, KernelAbstractions, Bessels, DispatchDoctor

using DispatchDoctor
@stable default_mode = "disable" begin

    function relu(x)
        if x > 0
            return x
        else
            return 0.0
        end
    end

    function visualize end
    function save_animation end
    export visualize, save_animation

    include("misc.jl")
    export overlap

    include("ortho_poly.jl")
    export laguerre, hermite

    include("initial_profiles.jl")
    export hg!, hg, diagonal_hg, diagonal_hg!, lg!, lg,
        rectangular_apperture, square, single_slit, double_slit, pupil, triangle

    include("free_propagation.jl")
    export free_propagation

    include("kerr_propagation.jl")
    export kerr_propagation

    include("holograms.jl")
    export BesselJ1, Simple, generate_hologram, generate_hologram!

    include("aberration_correction.jl")
    export lens!, lens, tilted_lens!, tilted_lens, zernike_polynomial!, zernike_polynomial

    include("precompile.jl")

end
end