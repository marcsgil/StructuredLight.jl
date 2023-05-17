using StructuredLight,LinearAlgebra
using Test

@testset "Laguerre Gauss" begin
    rs = LinRange(-10,10,1024)

    for p in 0:3, l in 0:3, z in (.1,.5,1)
        ψ1 = lg(rs,rs,z,p=p,l=l)
        ψ2 = free_propagation(lg(rs,rs,p = p,l = l), rs,rs,z)
        @test isapprox(ψ1 ⋅ ψ2 * (rs[2]-rs[1])^2, 1, atol=0.03)
    end
    
end

@testset "Hermite Gauss" begin
    rs = LinRange(-10,10,1024)

    for m in 0:3, n in 0:3, z in (.1,.5,1)
        ψ1 = hg(rs,rs,z,m=m,n=n)
        ψ2 = free_propagation(hg(rs,rs,m=m,n=n), rs,rs,z)
        @test isapprox(ψ1 ⋅ ψ2 * (rs[2]-rs[1])^2, 1, atol=0.03)
    end
    
end

@testset "Diagonal Hermite Gauss" begin
    rs = LinRange(-10,10,1024)

    for m in 0:3, n in 0:3, z in (.1,.5,1)
        ψ1 = diagonal_hg(rs,rs,z,m=m,n=n)
        ψ2 = free_propagation(diagonal_hg(rs,rs,m=m,n=n), rs,rs,z)
        @test isapprox(ψ1 ⋅ ψ2 * (rs[2]-rs[1])^2, 1, atol=0.03)
    end
    
end

@testset "Scalled Laguerre Gauss" begin
    rs = LinRange(-10,10,1024)

    for p in 0:3, l in 0:3, z in (.2,1,2)
        ψ1 = lg(2rs,2rs,z,p=p,l=l)
        ψ2 = free_propagation(lg(rs,rs,p = p,l = l), rs,rs,z,2)
        @test isapprox(ψ1 ⋅ ψ2 * 4 * (rs[2]-rs[1])^2, 1, atol=0.03)
    end
    
end

@testset "Scalled Hermite Gauss" begin
    rs = LinRange(-10,10,1024)

    for m in 0:3, n in 0:3, z in (.2,1,2)
        ψ1 = hg(2rs,2rs,z,m=m,n=n)
        ψ2 = free_propagation(hg(rs,rs,m=m,n=n), rs,rs,z,2)
        @test isapprox(ψ1 ⋅ ψ2 * 4 * (rs[2]-rs[1])^2, 1, atol=0.03)
    end
    
end

@testset "Scalled Diagonal Hermite Gauss" begin
    rs = LinRange(-10,10,1024)

    for m in 0:3, n in 0:3, z in (.2,1,2)
        ψ1 = diagonal_hg(2rs,2rs,z,m=m,n=n)
        ψ2 = free_propagation(diagonal_hg(rs,rs,m=m,n=n), rs,rs,z,2)
        @test isapprox(ψ1 ⋅ ψ2 * 4 * (rs[2]-rs[1])^2, 1, atol=0.03)
    end
    
end
