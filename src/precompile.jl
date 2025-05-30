using PrecompileTools

@setup_workload begin
    rs_list = [LinRange(-4.0f0, 4, 32), LinRange(-4.0, 4, 32)]


    @compile_workload begin
        for rs ∈ rs_list
            relative = lg(rs, rs, l=1, p=1)

            for z in (1, [1])
                free_propagation(relative, rs, rs, z)
                free_propagation(relative, rs, rs, 2, z)
                kerr_propagation(relative, rs, rs, z, 2)
            end

            holo = generate_hologram(relative, 255, 50, -50, Simple())
            holo = generate_hologram(relative, 255, 50, -50, BesselJ1())
        end
    end

end