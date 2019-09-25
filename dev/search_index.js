var documenterSearchIndex = {"docs":
[{"location":"api/#API-1","page":"API","title":"API","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"PiecewisePoly","category":"page"},{"location":"api/#Main.PiecewisePolynomials.PiecewisePoly","page":"API","title":"Main.PiecewisePolynomials.PiecewisePoly","text":"Type for function piecewise polynomial.\n\nPiecewisePoly(\n    polys::Vector{<:Poly},\n    knots::Vector{<:Real};\n    issafe::Bool = false\n    )\n\nArguments\n\npolys – polynomials\nknots – polynomials segments\nissafe – if true, correctness of agruments is not cheched\n\nFields\n\npolys::Vector{<:Poly} – function (type depends on the basis)\nknots::Vector{<:Real} – support of the function\n\n\n\n\n\n","category":"type"},{"location":"api/#","page":"API","title":"API","text":"derivative(poly::PiecewisePoly, x::Real)","category":"page"},{"location":"api/#Main.PiecewisePolynomials.derivative-Tuple{PiecewisePoly,Real}","page":"API","title":"Main.PiecewisePolynomials.derivative","text":"derivative(poly::PiecewisePoly, x::Real)\n\nReturns: the value of the first order derivative of polynomial at point x.\n\n\n\n\n\n","category":"method"},{"location":"api/#","page":"API","title":"API","text":"derivative(poly::PiecewisePoly; order::Int=1)","category":"page"},{"location":"api/#Main.PiecewisePolynomials.derivative-Tuple{PiecewisePoly}","page":"API","title":"Main.PiecewisePolynomials.derivative","text":"derivative(poly::PiecewisePoly; order::Int=1)::PiecewisePoly\n\nReturns: derivative of order order for given polynomial.\n\n\n\n\n\n","category":"method"},{"location":"api/#","page":"API","title":"API","text":"antiderivative(poly::PiecewisePoly)","category":"page"},{"location":"api/#Main.PiecewisePolynomials.antiderivative-Tuple{PiecewisePoly}","page":"API","title":"Main.PiecewisePolynomials.antiderivative","text":"antiderivative(poly::PiecewisePoly)::PiecewisePoly\n\nReturns: antiderivative for given polynomial.\n\n\n\n\n\n","category":"method"},{"location":"api/#","page":"API","title":"API","text":"integral(poly::PiecewisePoly, a::Real, b::Real)","category":"page"},{"location":"api/#Main.PiecewisePolynomials.integral-Tuple{PiecewisePoly,Real,Real}","page":"API","title":"Main.PiecewisePolynomials.integral","text":"integral(poly::PiecewisePoly, a::Real, b::Real)\n\nReturns: integral of given polynomial from point a to point b.\n\n\n\n\n\n","category":"method"},{"location":"#PiecewisePolynomials.jl-1","page":"Home","title":"PiecewisePolynomials.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This is documentation for PiecewisePolynomials.jl – a Julia package that allows to apply piecewise polynomials, integrate, differentiate them and do simple arithmetical operations.","category":"page"}]
}