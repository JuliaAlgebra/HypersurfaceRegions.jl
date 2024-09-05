using Chambers
using Test

@testset "2 concentric circles" begin
    @var x y
    f_1 = x^2 + y^2 - 1
    f_2 = x^2 + y^2 - 4
    f = System([f_1; f_2])
    R = chambers(f)
    R = chambers(f; show_progress = false)


    # functions for R
    @test ncritical_complex(R) == 9
    @test nchambers(R) == 3
    @test nbounded(R) == 2
    @test nunbounded(R) == 1

    @test length(bounded(R)) == 2
    @test length(unbounded(R)) == 1

    @test count(is_bounded, chambers(R)) == 2

    # affine chambers
    A = affine_chambers(f)
    @test nchambers(A) == 3

    # seed
    Rseed = chambers(f, seed = 0x801124df)
    @test nchambers(Rseed) == 3

    # setting the exponents
    Rs = chambers(f; s = [1; 1; -20])
    @test nchambers(Rs) == 3

    # membership
    R0 = membership(R, [0; 0])
    @test isa(R0, Chamber)
    @test Base.sign(R0) == [-1; -1]

    R1 = membership(R, [3 / 2; 0])
    @test Base.sign(R1) == [1; -1]

    # parameters 
    @var p
    f_1 = x^2 + y^2 - p
    f_2 = x^2 + y^2 - 4
    f = System([f_1; f_2], parameters = [p])
    R = chambers(f, target_parameters = [0.0])
    @test isa(R, ChambersResult)
end

@testset "a circle and a line" begin
    @var x y
    R = chambers([x^2 + y^2 - 1; x])
    p = projective_chambers(R)

    @test length(p) == 3

    D = filter(pi -> length(pi) == 2, p)
    @test !isempty(D)
    @test !isempty(D[1])
    @test all(Ci -> !is_bounded(Ci), D[1])
end

@testset "a circle and a parabola" begin
    @var x y
    f_1 = y^2 + x^2 - 1
    f_2 = y - x^2
    f = System([f_1; f_2])
    R = chambers(f)
    p = projective_chambers(R)

    @test nbounded(R) == 2
    @test nunbounded(R) == 1
    @test nundecided(R) == 1

    @test length(p) == 4
end

@testset "elliptope" begin
    @var x y z
    h = 2 * x * y * z - x^2 - y^2 - z^2 + 1
    f = System([h])
    R = chambers(f)

    @test nchambers(R) == 6

    E = [h; 1 - x; 1 + x; 1 - y; 1 + y; 1 - z; 1 + z]
    RE = chambers(E)
    @test nchambers(RE) == 43
end

@testset "3 ellipsoids in R^3" begin
    @var x y z

    f = [x^2 + y^2 + z^2 - 1, x^2 + y^2 + z^2 - 4, 100 * x^2 + 100 * y^2 + z^2 - 9]
    R = chambers(f)

    @test nchambers(R) == 8
end

@testset "6 lines in R^2" begin
    @var x y
    coeff_matrix = randn(3, 6)
    f_list = map(v -> v[1] * x + v[2] * y + v[3], eachcol(coeff_matrix))
    R = chambers(f_list)

    @test isa(R, ChambersResult)
end


@testset "a net plane of sextics" begin
    @var a b
    f = [a, a + 1, 3a - 1, 3a + b + 6, 3a + b - 3,
            9a^3 - 3a^2*b + a*b^2 - 3a - b + 2]
    C = chambers(f)
    P = projective_chambers(C)
    K = length.(P)

    @test count(k -> k == 1, K) == 13
    @test count(k -> k == 2, K) == 2
end