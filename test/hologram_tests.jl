using StructuredLight: find_zero, inverse, closest_index, zero_order_interpolation

# Test for find_zero function
@testset "find_zero tests" begin
    f(x) = x - 3
    @test isapprox(find_zero(f, 0, 5), 3, rtol=1e-6)

    g(x) = x^2 - 4
    @test isapprox(find_zero(g, 0, 5), 2, rtol=1e-6)
    @test isapprox(find_zero(g, -5, 0), -2, rtol=1e-6)
end

# Test for inverse function
@testset "inverse tests" begin
    f(x) = 2x + 1
    @test isapprox(inverse(f, 3, 0, 5), 1, rtol=1e-6)

    g(x) = x^2
    for i in 1:10
        @test isapprox(inverse(g, i, 0, 10), √i, rtol=1e-6)
    end
end


# Test for closest_index function
@testset "closest_index tests" begin
    grid = 1:5
    @test closest_index(2.4, grid) == 2
    @test closest_index(4.6, grid) == 5
end

# Test for zero_order_interpolation function
@testset "zero_order_interpolation tests" begin
    xs = 1:5
    ys = 10:10:50
    @test zero_order_interpolation(2.4, xs, ys) == 20
    @test zero_order_interpolation(4.6, xs, ys) == 50
end