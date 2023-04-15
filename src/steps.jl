function dispersion_step!(ψ,phases,plan,iplan)
    plan*ψ
    map!(*,ψ,ψ,phases)
    iplan*ψ
end

evolve_phase(ψ,phase) = cis(phase) * ψ