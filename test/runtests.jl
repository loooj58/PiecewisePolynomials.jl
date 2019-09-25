using Test, Polynomials

# include("../src/PiecewisePolynomials.jl")
include("../src/utils.jl")
include("../src/implementation.jl")
# using .PiecewisePolynomials

polys1 = [Poly([1, 2, 3]), Poly([3, 2, 1]), Poly([2])]
knots1 = [-1, 4, 5, 11.5]

polys2 = [Poly([5, 6, 7, 2, 3]), Poly([2, -6, 1]), Poly([2, -10, 2]), Poly([2, -10, 2])]
knots2 = [1, 3, 5, 10, 15]

@testset "External" begin
    @test @returntrue my_poly1 = PiecewisePoly(polys1, knots1)
    @test @returntrue my_poly2 = PiecewisePoly(polys2, knots2)
    @test @returntrue derivative(my_poly1, order=3)(2)
    @test @returntrue derivative(my_poly1, order=0)
    @test @returntrue antiderivative(my_poly2)
    @test @returntrue antiderivative(my_poly2)(5)
    @test @returntrue integral(my_poly1, 2, 4)
    @test @returntrue show(my_poly2)
    @test @returntrue my_poly1 + my_poly2
    @test @returntrue my_poly1 - my_poly2
    @test @returntrue my_poly1 * my_poly2
end

@testset "Internal" begin
    @test @returntrue find_poly(my_poly1, 2)
    @test @returntrue merge_knots(knots1, knots2)
end
