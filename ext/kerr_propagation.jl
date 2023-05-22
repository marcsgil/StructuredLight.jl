function StructuredLight.dispersion_step!(ψ::CuArray,kernel::CuArray,plan,iplan)
    plan*ψ
    ψ .*= kernel
    iplan*ψ
end

function StructuredLight.type_2A_step!(ψ::CuArray,kernel::CuArray,plan,iplan,factor,repetitions)
    if !iszero(repetitions)
        @tullio ψ[i,j] *= cis( factor * abs2(ψ[i,j]) / 2)
        for _ in 1:repetitions-1
            StructuredLight.dispersion_step!(ψ,kernel,plan,iplan)
            @tullio ψ[i,j] *= cis( factor * abs2(ψ[i,j]) )
        end
        StructuredLight.dispersion_step!(ψ,kernel,plan,iplan)
        @tullio ψ[i,j] *= cis( factor * abs2(ψ[i,j]) / 2)
    end
end

function StructuredLight.kerr_propagation_loop!(dest::CuArray,ψ₀::CuArray,kernel::CuArray,phases::CuArray,z,divisions,g,k,plan,iplan)
    Δz = z / divisions
    phase_evolution_factor = g * Δz / 2k

    @tullio kernel[i,j] = cis( Δz * phases[i,j] )

    StructuredLight.type_2A_step!(ψ₀,kernel,plan,iplan,phase_evolution_factor,divisions)

    fftshift!(dest, ψ₀)
end

function StructuredLight.kerr_propagation(ψ₀::CuArray,xs,ys,zs,total_steps;k=1,g=1)
    Zs = vcat(0,Array(zs))

    steps = StructuredLight.distribute(total_steps,Zs)

    result = similar(ψ₀,size(ψ₀)...,length(zs))

    plan = plan_fft!(ψ₀)
    iplan = plan_ifft!(ψ₀)

    ψ = ifftshift(ψ₀)

    qxs = StructuredLight.reciprocal_grid(xs,shift=true) |> CuArray
    qys = StructuredLight.reciprocal_grid(ys,shift=true) |> CuArray

    @tullio phases[i,j] := - ( qxs[i]^2 + qys[j]^2 ) / 2k
    kernel = similar(ψ₀)

    for (i,divisions) in enumerate(steps)
        StructuredLight.kerr_propagation_loop!(view(result,:,:,i),ψ,kernel,phases,Zs[i+1] - Zs[i],divisions,g,k,plan,iplan)
    end

    result
end