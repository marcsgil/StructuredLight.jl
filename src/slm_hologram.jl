mymod(x) = mod(x + π,2π) - π

function generate_hologram(desired,input,xs,Λ)
    relative = desired ./ input
    M = maximum(abs,relative)

    @tullio result[i,j] := abs(relative[i,j]) / M * mymod(angle(relative[i,j]) + 2π*xs[j]/Λ)
end