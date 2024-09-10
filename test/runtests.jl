using HypersurfaceRegions
using Test

@testset "2 concentric circles" begin
    @var x y
    f_1 = x^2 + y^2 - 1
    f_2 = x^2 + y^2 - 4
    f = System([f_1; f_2])
    R = regions(f)
    R = regions(f; show_progress = false)


    # functions for R
    @test ncritical_complex(R) == 9
    @test nregions(R) == 3
    @test nbounded(R) == 2
    @test nunbounded(R) == 1

    @test length(bounded(R)) == 2
    @test length(unbounded(R)) == 1

    @test count(is_bounded, regions(R)) == 2

    # affine regions
    A = affine_regions(f)
    @test nregions(A) == 3

    # seed
    Rseed = regions(f, seed = 0x801124df)
    @test nregions(Rseed) == 3

    # setting the exponents
    Rs = regions(f; s = [1; 1; -20])
    @test nregions(Rs) == 3

    # membership
    R0 = membership(R, [0; 0])
    @test isa(R0, Region)
    @test Base.sign(R0) == [-1; -1]

    R1 = membership(R, [3 / 2; 0])
    @test Base.sign(R1) == [1; -1]

    # parameters 
    @var p
    f_1 = x^2 + y^2 - p
    f_2 = x^2 + y^2 - 4
    f = System([f_1; f_2], parameters = [p])
    R = regions(f, target_parameters = [0.0])
    @test isa(R, RegionsResult)
end

@testset "a circle and a line" begin
    @var x y
    R = regions([x^2 + y^2 - 1; x])
    p = projective_regions(R)

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
    R = regions(f)
    p = projective_regions(R)

    @test nbounded(R) == 2
    @test nunbounded(R) == 1
    @test nundecided(R) == 1

    @test length(p) == 4
end

@testset "elliptope" begin
    @var x y z
    h = 2 * x * y * z - x^2 - y^2 - z^2 + 1
    f = System([h])
    R = regions(f)

    @test nregions(R) == 6

    E = [h; 1 - x; 1 + x; 1 - y; 1 + y; 1 - z; 1 + z]
    RE = regions(E)
    @test nregions(RE) == 43
end

@testset "3 ellipsoids in R^3" begin
    @var x y z

    f = [x^2 + y^2 + z^2 - 1, x^2 + y^2 + z^2 - 4, 100 * x^2 + 100 * y^2 + z^2 - 9]
    R = regions(f)

    @test nregions(R) == 8
end

@testset "6 lines in R^2" begin
    @var x y
    coeff_matrix = randn(3, 6)
    f_list = map(v -> v[1] * x + v[2] * y + v[3], eachcol(coeff_matrix))
    R = regions(f_list)

    @test isa(R, RegionsResult)
end


@testset "a net plane of sextics" begin
    @var a b
    f = [a, a + 1, 3a - 1, 3a + b + 6, 3a + b - 3,
            9a^3 - 3a^2*b + a*b^2 - 3a - b + 2]
    C = regions(f)
    P = projective_regions(C)
    K = length.(P)

    @test count(k -> k == 1, K) == 13
    @test count(k -> k == 2, K) == 2
end

@testset "one paraboloid" begin
    @var  x y z
    f = [3 + x + 3*y - z + (1 + 2*x + 4*y - 4*z)^2 + (2 + 3*x + 2*y + 3*z)^2]
    R = regions(f)
    @test nundecided(R) == 1
end

@testset "two paraboloids" begin
    @var  x y z
    f = [3 + x + 3*y - z + (1 + 2*x + 4*y - 4*z)^2 + (2 + 3*x + 2*y + 3*z)^2,
        3 + x + 3*z + (3 - 3*x - 2*z)^2 + (3 + 3*x + 3*y + 4*z)^2
        ]
    R = regions(f)
    @test nundecided(R) == 2
end

@testset "three paraboloids" begin
    @var  x y z
    f = [3 + x + 3*y - z + (1 + 2*x + 4*y - 4*z)^2 + (2 + 3*x + 2*y + 3*z)^2,
        3 + x + 3*z + (3 - 3*x - 2*z)^2 + (3 + 3*x + 3*y + 4*z)^2,
        2 - 2*x - 2*y - 3*z + (2 - x + 4*z)^2 + (2 + 3*x + y + 2*z)^2
        ]
    R = regions(f)
    @test nundecided(R) == 3
end

@testset "four paraboloids" begin
    @var  x y z
    f = [3 + x + 3*y - z + (1 + 2*x + 4*y - 4*z)^2 + (2 + 3*x + 2*y + 3*z)^2,
        3 + x + 3*z + (3 - 3*x - 2*z)^2 + (3 + 3*x + 3*y + 4*z)^2,
        2 - 2*x - 2*y - 3*z + (2 - x + 4*z)^2 + (2 + 3*x + y + 2*z)^2,
        1 - 3*x + 3*y - 3*z + (1 - 2*y + 2*z)^2 + (2 + x + 4*y)^2
        ]
    R = regions(f)
    @test undecided(R) == 4
end