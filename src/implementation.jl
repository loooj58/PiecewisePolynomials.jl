using Polynomials


"""
Type for function piecewise polynomial.

```julia
PiecewisePoly(
    polys::Vector{<:Poly},
    knots::Vector{<:Real};
    issafe::Bool = false
    )
```

**Arguments**

* `polys` -- polynomials
* `knots` -- polynomials segments
* `issafe` -- if true, correctness of agruments is not cheched

**Fields**

* `polys::Vector{<:Poly}` -- function (type depends on the basis)
* `knots::Vector{<:Real}` -- support of the function
"""
struct PiecewisePoly
    polys::Vector{<:Poly}
    knots::Vector{<:Real}

    function PiecewisePoly(
        polys::Vector{<:Poly},
        knots::Vector{<:Real};
        issafe::Bool = false
        )
        if issafe
            return new(polys, knots)
        else
            if length(polys) + 1 != length(knots)
                error(
                "length(poly) + 1 != length(knots), $(length(polys)), $(length(knots))"
                )
            end
            if !issorted(knots)
                error("Knots should be sorted")
            end
            return new(polys, knots)
        end
    end

end


function find_poly(
    poly::PiecewisePoly, x::Real;
    segmetize::Bool=false, atol=1e-7, out_of_segment="error"
    )

    if poly.knots[1] > x
        if out_of_segment == "error"
            error("x should be inside the segment")
        end
        return Poly([0])
    end
    for (index, value) in enumerate(poly.knots)
        if value > x
            return poly.polys[index - 1]
        end
    end

    if segmetize == true
        if isapprox(x, poly.knots[end], atol=atol) && x <= poly.knots[end]
            return poly.polys[end]
        end
    end
    if out_of_segment == "error"
        error("x should be inside the segment")
    end
    return Poly([0])
end


"""
Computes the value of function in point `x`
"""
(poly::PiecewisePoly)(x::Real) = find_poly(
    poly, x, out_of_segment="zero", segmetize=true
    )(x)


"""
```julia
derivative(poly::PiecewisePoly, x::Real)
```

**Returns:** the value of the first order derivative of polynomial at point `x`.
"""
derivative(poly::PiecewisePoly, x::Real) = polyder(
    find_poly(poly, x, out_of_segment="zero", segmetize=true)
    )(x)


"""
```julia
derivative(poly::PiecewisePoly; order::Int=1)::PiecewisePoly
```

**Returns:** derivative of order `order` for given polynomial.
"""
function derivative(poly::PiecewisePoly; order::Int=1)::PiecewisePoly
    if order < 0
        error("The order of derivative shiuld be non-negative.")
    end
    if order == 0
        return poly
    end
    if order == 1
        return PiecewisePoly(map(polyder, poly.polys), poly.knots, issafe=true)
    end
    return derivative(derivative(poly, order=order-1))
end


"""
```julia
antiderivative(poly::PiecewisePoly)::PiecewisePoly
```

**Returns:** antiderivative for given polynomial.
"""
function antiderivative(poly::PiecewisePoly)::PiecewisePoly
    return PiecewisePoly(map(polyint, poly.polys), poly.knots, issafe=true)
end


"""
```julia
integral(poly::PiecewisePoly, a::Real, b::Real)
```

**Returns:** integral of given polynomial from point `a` to point `b`.
"""
function integral(poly::PiecewisePoly, a::Real, b::Real)
    if !(a >= poly.knots[1] && a <= b && b <= poly.knots[end])
        Base.error("a or b not in the segment")
    end
    integral = 0.
    for (index, value) in enumerate(poly.polys)
        if a < poly.knots[index + 1] && b > poly.knots[index]
            s = polyint(value, poly.knots[index], poly.knots[index + 1])
            if a > poly.knots[index]
                s -= polyint(value, poly.knots[index], a)
            end
            if b < poly.knots[index + 1]
                s -= polyint(value, b, poly.knots[index + 1])
            end
            integral += s
        end
    end
    return integral
end

function merge_knots(a::Vector{<:Real}, b::Vector{<:Real}, atol::Real=1e-7)
    i = 1
    j = 1
    result = []
    while (i != length(a) + 1 && j != length(b) + 1)
        if isapprox(a[i], b[j], atol=atol)
            append!(result, a[i])
            i += 1
            j += 1
        elseif a[i] < b[j]
            append!(result, a[i])
            i += 1
        else
            append!(result, b[j])
            j += 1
        end
    end
    while i != length(a) + 1
        append!(result, a[i])
        i += 1
    end
    while j != length(b) + 1
        append!(result, b[j])
        j += 1
    end
    return result
end


import Base.+
function (+)(poly1::PiecewisePoly, poly2::PiecewisePoly)::PiecewisePoly
    new_knots = merge_knots(poly1.knots, poly2.knots)
    new_polys = []
    for (l, r) in collect(zip(new_knots[1:end-1], new_knots[2:end]))
        mid = (l + r) / 2
        append!(new_polys, [
        find_poly(poly1, mid, out_of_segment="zero", segmetize=true) +
        find_poly(poly2, mid, out_of_segment="zero", segmetize=true)
        ])
    end
    return PiecewisePoly(
        convert(Array{Poly}, new_polys), convert(Array{Float64}, new_knots),
        issafe=true
        )
end

import Base.-
function (-)(poly1::PiecewisePoly, poly2::PiecewisePoly)::PiecewisePoly
    return poly1 + PiecewisePoly(-poly2.polys, poly2.knots, issafe=true)
end

import Base.*
function (*)(poly1::PiecewisePoly, poly2::PiecewisePoly)::PiecewisePoly
    new_knots = merge_knots(poly1.knots, poly2.knots)
    new_polys = []
    for (l, r) in collect(zip(new_knots[1:end-1], new_knots[2:end]))
        mid = (l + r) / 2
        append!(new_polys, [
        find_poly(poly1, mid, out_of_segment="zero", segmetize=true) *
        find_poly(poly2, mid, out_of_segment="zero", segmetize=true)
        ])
    end
    return PiecewisePoly(
        convert(Array{Poly}, new_polys), convert(Array{Float64}, new_knots),
        issafe=true
        )
end


function Base.show(io::IO, poly::PiecewisePoly)
    function make_polystring(poly)
        a = poly.a
        second = ""
        if length(a) == 1
            return "$(a[1])"
        elseif length(a) == 2
            return "$(a[1]) + $(a[2])x"
        end
        return string(
            a[1],
            " + $(a[2])x",
            [" + $(value)x^$(key+1)" for (key, value) in enumerate(a[3:end])]...
            )
    end

    poly_string = [make_polystring(poly_) for poly_ in poly.polys]
    knots_string = [" from $(knot_begin) to $(knot_end)\n"
        for (knot_begin, knot_end)
        in collect(zip(poly.knots[1:end-1], poly.knots[2:end]))
        ]
    list_of_strings = [string(poly_string_, knots_string_)
        for (poly_string_, knots_string_)
        in zip(poly_string, knots_string)
        ]
    println(io, list_of_strings...)
end
