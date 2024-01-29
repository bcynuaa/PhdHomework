# author: bcynuaa
# date: 2023/01/12

###################################################################################################

struct Case
    alpha::Float64
    b::Float64
    a::Float64
    boundary_value::Float64
end

function Case(
    angle_type::String="deg",
    alpha::Float64=30.,
    b::Float64=1.,
    a::Float64=1.,
    boundary_value::Float64=0.
)::Case
    if angle_type == "deg"
        alpha *= pi / 180.
    end
    return Case(
        alpha,
        b,
        a,
        boundary_value
    );
end

const error::Float64 = 1e-4;
function whether_on_boundary(case::Case, xy::Vector{Float64})::Bool
    x::Float64, y::Float64 = xy;
    if isapprox(abs(x - case.b), 0.; atol=error) == true
        return true;
    elseif isapprox(abs(y - tan(case.alpha) * x), 0.; atol=error) == true
        return true;
    elseif isapprox(abs(y + tan(case.alpha) * x), 0.; atol=error) == true
        return true;
    else
        return false;
    end
end