function StructuredLight.reciprocal_quadratic_phase(qxs::CuArray,qys::CuArray,k)
    @tullio result[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
end

function StructuredLight.direct_quadratic_phase(xs::CuArray,ys::CuArray,k)
    @tullio result[i,j] := k * ( xs[i]^2 + ys[j]^2 ) / 2
end

function StructuredLight.apply_reciprocal_phases(ψ₀::CuArray,phases,zs)
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( phases[i,j] * zs[l] )
end

function StructuredLight.apply_reciprocal_phases!(ψ::CuArray,phases,zs,scaling)
    @tullio ψ[i,j,l] *= cis( phases[i,j] * zs[l] / scaling[l] )
end

function StructuredLight.apply_direct_phases(ψ₀::CuArray,phases,zs,scaling)
    @tullio ψ[i,j,l] := ψ₀[i,j] * cis( phases[i,j] * ( 1 - scaling[l] ) / zs[l] ) / scaling[l]
end

function StructuredLight.apply_direct_phases!(ψ::CuArray,phases,zs,scaling)
    @tullio ψ[i,j,l] *= cis( - phases[i,j] * ( 1 - scaling[l] ) * scaling[l] / zs[l])
end

function StructuredLight.fourier_propagation_kernel!(dest::CuArray,phases::CuArray,z)
    @tullio dest[i,j] = cis( z * phases[i,j] )
end

function StructuredLight.dispersion_step!(ψ::CuArray,kernel::CuArray,plan,iplan)
    plan*ψ
    ψ .*= kernel
    iplan*ψ
end

function StructuredLight.apply_kerr_phase!(ψ::CuArray, phase_factor)
    @tullio ψ[i,j] *= cis( phase_factor * abs2(ψ[i,j]) )
end