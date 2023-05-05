# Miscellany

Here we show a few small functionalities that didn't fit anywhere else.

```@docs
overlap
```

```julia
using StructuredLight

rs = LinRange(-5,5,256) 
ψ₁ = hg(rs,rs)
ψ₂ = hg(rs,rs,m=1)

overlap(ψ₁,ψ₁,rs,rs) #The modes are normalized
overlap(ψ₁,ψ₂,rs,rs) #Modes with different indices are orthogonal

#Free propagation is an unitary transformation
zs = LinRange(0,.5,8)
ψ₁s = free_propagation(ψ₁,rs,rs,zs)
ψ₂s = free_propagation(ψ₂,rs,rs,zs)

overlap(ψ₁s,ψ₁s,rs,rs) 
overlap(ψ₁s,ψ₂s,rs,rs) 
```