include("../src/PiecewisePolynomials.jl")
using .PiecewisePolynomials
using Polynomials


polys = [Poly([1, 2, 3]), Poly([3, 2, 1]), Poly([2])]
knots = [1, 3, 5, 10]
my_poly = PiecewisePoly(polys, knots)

polys1 = [Poly([5, 6, 7, 2, 3]), Poly([2, -6, 1]), Poly([2, -10, 2]), Poly([2, -10, 2])]
knots1 = [1, 3, 5, 10, 15]
my_poly1 = PiecewisePoly(polys1, knots1)


println(derivative(my_poly, order=3)(2))
println(antiderivative(my_poly1)(5))
println(integral(my_poly, 2, 4))
Base.show(my_poly + my_poly1)

my_poly + my_poly1

my_poly - my_poly1

my_poly * my_poly1
