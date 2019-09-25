macro returntrue(f::Expr)
    try
        eval(f)
    catch e
        error(e)
        return false
    end
    return true
end
