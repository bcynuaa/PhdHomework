# author: bcynuaa
# date: 2023/01/17

struct Case
    L::Float64
    Es::Float64
    d1::Float64
    d2::Float64
    J1::Float64
    J2::Float64
    q::Float64
    F::Float64
    M::Float64
    devide::Int64
    h::Float64
    N::Int64
    x::Vector{Float64}
    middle_ID::Int64
end

function moment_of_inertia(d::Float64)::Float64
    return pi * d^4 / 64.;
end

function Case(
    L::Float64 = 0.12,
    Es::Float64 = 200e9,
    d1::Float64 = 0.03,
    d2::Float64 = 0.02,
    q::Float64 = -200.,
    F::Float64 = -1000.,
    M::Float64 = 2000.,
    devide::Int64 = 10
)::Case
    N::Int64 = 2*devide + 1;
    h::Float64 = L / devide;
    middle_ID::Int64 = devide + 1;
    x::Vector{Float64} = Vector(LinRange(0., 2*L, N));
    return Case(
        L,
        Es,
        d1,
        d2,
        moment_of_inertia(d1),
        moment_of_inertia(d2),
        q,
        F,
        M,
        devide,
        h,
        N,
        x,
        middle_ID
    );
end