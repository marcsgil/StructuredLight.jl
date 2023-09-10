using GLMakie, StructuredLight, FourierTools, Tullio
include("src/dft_utils.jl")

round2(x) = round(x,sigdigits=2)

function interactive_visualization(ψ₀, rs, zs)
    qs = reciprocal_grid(rs, shift=true)
    @tullio phases[i, j] := -(qs[j]^2 + qs[i]^2) / 2
    fig = Figure()
    ax = Axis(fig[1, 1], aspect=1, titlesize=32)
    I = Observable(permutedims(abs2.(ψ₀),(2,1)))
    hidedecorations!(ax)
    image!(ax, I, colormap=:hot)
    plan, iplan = plan_fft!(ψ₀), plan_ifft!(ψ₀)

    sl = Slider(fig[2, 1], range = zs, startvalue = first(zs))

    lift(sl.value) do z
        permutedims!(I[], step(ψ₀, z, phases, plan, iplan), (2, 1))
        I[] = I[]
        ax.title = L"z = %$(round2(z))"
    end

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            #=
            if event.key == Keyboard.right
                if sl.selected_index[] < length(zs)
                    sl.selected_index[] += 1
                end
            elseif event.key == Keyboard.left
                if sl.selected_index[] > 1
                    sl.selected_index[] -= 1
                end
            end
            =#

            @async while Keyboard.right ∈ events(fig).keyboardstate
                if sl.selected_index[] < length(zs)
                    sl.selected_index[] += 1
                end
                sleep(0.1)
            end

            @async while Keyboard.left ∈ events(fig).keyboardstate
                if sl.selected_index[] > 1
                    sl.selected_index[] -= 1
                end
                sleep(0.1)
            end
        end
    end

    fig
end

function step(ψ₀, z, phases, plan, iplan)
    cache = ifftshift(ψ₀)
    plan * cache
    @tullio cache[i, j] *= cis(phases[i, j] * z)
    iplan * cache
    abs2.(fftshift(cache))
end
##
rs = LinRange(-5, 5, 512)
zs = LinRange(0, 2, 128)
ψ₀ = hg(rs, rs, m=1)

interactive_visualization(ψ₀, rs, zs)
##
fig = Figure()
label = Slider(fig[2, 1], range = zs, startvalue = first(zs))
propertynames(label)
label.text=L"2"