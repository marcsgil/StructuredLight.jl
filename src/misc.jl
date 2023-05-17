"""
    overlap(ψ₁,ψ₂,xs,ys)

Calculate ∫ `ψ₁` `conj(ψ₂)` dx dy.

`xs` and `ys` are the grids over which `ψ₁` and `ψ₂` are calculated.
"""
overlap(ψ₁,ψ₂,xs,ys) = interval(xs) * interval(ys) * ψ₁ ⋅ ψ₂

"""
    overlap(ψ₁::AbstractArray{T1,3},ψ₂::AbstractArray{T2,3},xs,ys) where {T1,T2}

Calculate ∫ `ψ₁[:,:,j]` `conj(ψ₂[:,:,j])` dx dy, where j runs over all indices of the third dimension. 

The output is a vector.

`xs` and `ys` are the grids over which `ψ₁` and `ψ₂` are calculated.
"""
function overlap(ψ₁::AbstractArray{T1,3},ψ₂::AbstractArray{T2,3},xs,ys) where {T1,T2}
    map( (a,b) -> overlap(a,b,xs,ys), eachslice(ψ₁,dims=3), eachslice(ψ₂,dims=3) )
end