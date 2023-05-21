

function fourier_propagation_kernel!(dest,phases,z)
    @tullio dest[i,j] = cis( z * phases[i,j] )
end

function dispersion_step!(ψ,kernel,plan,iplan)
    plan*ψ
    @tullio ψ[i,j] *= kernel[i,j]
    iplan*ψ
end

function apply_kerr_phase!(ψ, phase_factor)
    @tullio ψ[i,j] *= cis( phase_factor * abs2(ψ[i,j]) )
end