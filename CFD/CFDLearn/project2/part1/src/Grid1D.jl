# author: bcynuaa
# date: 2022/11/06

module Grid1D # module begin

###################################################################################################

struct UniformGrid1D
    dx::Float64
    x_begin::Float64
    x_end::Float64
    N::Int64
    x::Vector{Float64}
    dt::Float64
    t_end::Float64
    Step::Int64
    t::Vector{Float64}
end

function UniformGrid1D(
    x_begin::Float64,
    x_end::Float64,
    dx::Float64,
    t_end::Float64,
    dt::Float64
)::UniformGrid1D
    x::Vector{Float64} = Vector(x_begin: dx: x_end);
    t::Vector{Float64} = Vector(0.: dt: t_end);
    N::Int64 = length(x);
    Step::Int64 = length(t) - 1;
    return UniformGrid1D(
        dx, x_begin, x_end, N, x,
        dt, t_end, Step, t
    );
end

function UniformGrid1D(
    x_begin::Float64,
    x_end::Float64,
    N::Int64,
    t_end::Float64,
    Step::Int64
)::UniformGrid1D
    dx::Float64 = (x_end - x_begin) / N;
    N += 1;
    x::Vector{Float64} = Vector(x_begin: dx: x_end);
    dt::Float64 = t_end / Step;
    t::Vector{Float64} = Vector(0.: dt: t_end);
    return UniformGrid1D(
        dx, x_begin, x_end, N, x,
        dt, t_end, Step, t
    );
end

###################################################################################################

export UniformGrid1D;

end # module end