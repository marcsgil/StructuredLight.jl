using PrecompileTools

@setup_workload begin
    rs_list = [LinRange(-4.0f0, 4, 32), LinRange(-4.0, 4, 32)]


    @compile_workload begin
        for rs ∈ rs_list
            desired = lg(rs, rs, l=1, p=1)

            for z in (1, [1])
                free_propagation(desired, rs, rs, z)
                free_propagation(desired, rs, rs, 2, z)
                kerr_propagation(desired, rs, rs, z, 2)
            end

            input = hg(rs, rs; w=10)

            holo = generate_hologram(desired, input, rs, rs, 255, 50, -50, Simple)
            holo = generate_hologram(desired, input, rs, rs, 255, 50, -50, BesselJ1)
        end
    end

end