function generate_hologram(desired,input,xs,Λ)
    relative = desired 
    M = maximum(abs,relative)

    #@tullio result[i,j] := abs(relative[i,j]) / M *( mod( (angle(relative[i,j]) + π + ),2π)) - π

    @tullio result[i,j] := mod(angle(relative[i,j]) + 2π*xs[j]/Λ,2π)
end