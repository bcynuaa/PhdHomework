"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-12 19:50:43
 # @ license: MIT
 # @ description: contains the advection operation
 """

include("Case.jl");

function binaryInterpolation(
    p::Float64, q::Float64,
    N1::Int64, M1::Int64,
    u_list::Vector{Float64},
)::Float64
    ul::Float64 = (N1+1-p) * u_list[1] + (p-N1) * u_list[2];
    uu::Float64 = (N1+1-p) * u_list[4] + (p-N1) * u_list[3];
    return (M1+1-q) * ul + (q-M1) * uu;
end

function advection!(
    case::Case,
    w1u::Matrix{Float64},
    w1v::Matrix{Float64},
)::Int64
    w2u::Matrix{Float64} = zeros(Float64, case.N_, case.N_);
    w2v::Matrix{Float64} = zeros(Float64, case.N_, case.N_);
    u::Float64 = 0.;
    v::Float64 = 0.;
    for j = 2: case.N_ - 1
        for i = 2: case.N_ - 1
            p::Float64 = i - w1u[j, i] * case.dt_ / case.dx_;
            q::Float64 = j - w1v[j, i] * case.dt_ / case.dx_;
            N1::Int64 = floor(p);
            M1::Int64 = floor(q);
            if (M1 + 1) > case.N_
                u = case.u_upper_;
                v = 0.;
            elseif N1 > 0 && N1 < case.N_ && M1 > 0 && M1 < case.N_
                u_list::Vector{Float64} = [w1u[M1, N1], w1u[M1, N1+1], w1u[M1+1, N1+1], w1u[M1+1, N1]];
                u = binaryInterpolation(p, q, N1, M1, u_list);
                v_list::Vector{Float64} = [w1v[M1, N1], w1v[M1, N1+1], w1v[M1+1, N1+1], w1v[M1+1, N1]];
                v = binaryInterpolation(p, q, N1, M1, v_list);
            else
                u = 0.;
                v = 0.;
            end
            w2u[j, i] = u;
            w2v[j, i] = v;
        end
    end
    w1u[2: case.N_-1, 2: case.N_-1] .= w2u[2: case.N_-1, 2: case.N_-1];
    w1v[2: case.N_-1, 2: case.N_-1] .= w2v[2: case.N_-1, 2: case.N_-1];
    return 0;
end