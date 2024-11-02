using Test

# Test rectangular_apperture
@testset "rectangular_apperture" begin
    x_inside = [0.0, 0.5]
    y_inside = [0.0, 0.5]
    x_outside = [1.0, -1.0]
    y_outside = [1.0, -1.0]
    a, b = 1.0, 1.0

    @test all(rectangular_apperture(x_inside, y_inside, a, b))
    @test !any(rectangular_apperture(x_outside, y_outside, a, b))

    x = LinRange(-1.0, 1.0, 10)
    y = LinRange(-1.0, 1.0, 20)
    result = rectangular_apperture(x, y, a, b)
    @test size(result) == (10, 20)
end

# Test square
@testset "square" begin
    x_inside = [0.0, 0.5]
    y_inside = [0.0, 0.5]
    x_outside = [1.0, -1.0]
    y_outside = [1.0, -1.0]
    l = 1.0

    @test all(square(x_inside, y_inside, l))
    @test !any(square(x_outside, y_outside, l))

    x = LinRange(-1.0, 1.0, 10)
    y = LinRange(-1.0, 1.0, 20)
    result = square(x, y, l)
    @test size(result) == (10, 20)
end

# Test single_slit
@testset "single_slit" begin
    x_inside = [0.0, 0.5]
    y_inside = [0.0, 0.5]
    x_outside = [1.0, -1.0]
    y_outside = [1.0, -1.0]
    a = 1.0

    @test all(single_slit(x_inside, y_inside, a))
    @test !any(single_slit(x_outside, y_outside, a))

    x = LinRange(-1.0, 1.0, 10)
    y = LinRange(-1.0, 1.0, 20)
    result = single_slit(x, y, a)
    @test size(result) == (10, 20)
end

# Test double_slit
@testset "double_slit" begin
    x_inside = [-0.5, 0.5]
    y_inside = [0.0, 0.5]
    x_outside = [0, 3]
    y_outside = [0, 3]
    a = 1.0
    d = 2.0

    @test all(double_slit(x_inside, y_inside, a, d))
    @test !any(double_slit(x_outside, y_outside, a, d))

    x = LinRange(-2.0, 2.0, 10)
    y = LinRange(-1.0, 1.0, 20)
    result = double_slit(x, y, a, d)
    @test size(result) == (10, 20)
end

# Test pupil
@testset "pupil" begin
    x_inside = [0.0, 0.5]
    y_inside = [0.0, 0.5]
    x_outside = [1.0, -1.0]
    y_outside = [1.0, -1.0]
    r = 1.0

    @test all(pupil(x_inside, y_inside, r))
    @test !any(pupil(x_outside, y_outside, r))

    x = LinRange(-1.0, 1.0, 10)
    y = LinRange(-1.0, 1.0, 20)
    result = pupil(x, y, r)
    @test size(result) == (10, 20)
end

# Test triangle
@testset "triangle" begin
    x_inside = [0.0, -0.3]
    y_inside = [0.0, -0.2]
    x_outside = [1.0, -1.0]
    y_outside = [1.0, -1.0]
    side_length = 1.0

    @test all(triangle(x_inside, y_inside, side_length))
    @test !any(triangle(x_outside, y_outside, side_length))

    x = LinRange(-1.0, 1.0, 10)
    y = LinRange(-1.0, 1.0, 20)
    result = triangle(x, y, side_length)
    @test size(result) == (10, 20)
end