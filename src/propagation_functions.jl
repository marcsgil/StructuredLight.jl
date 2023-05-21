function reciprocal_quadratic_phase(qxs,qys,k)
    @tullio result[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
end

function direct_quadratic_phase(xs,ys,k)
    @tullio result[i,j] := k * ( xs[i]^2 + ys[j]^2 ) / 2
end

function apply_reciprocal_phases(ψ₀,phases,zs)
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( phases[i,j] * zs[l] )
end

function apply_reciprocal_phases!(ψ,phases,z)
    @tullio ψ[i,j] *= cis( phases[i,j] * z )
end

function apply_reciprocal_phases!(ψ,phases,zs,scaling)
    @tullio ψ[i,j,l] *= cis( phases[i,j] * zs[l] / scaling[l] )
end

function apply_direct_phases(ψ₀,phases,zs,scaling)
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( phases[i,j] * ( 1 - scaling[l] ) / zs[l] ) / scaling[l]
end

function apply_direct_phases!(ψ,phases,zs,scaling)
    @tullio ψ[i,j,l] *= cis( - phases[i,j] * ( 1 - scaling[l] ) * scaling[l] / zs[l])
end

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