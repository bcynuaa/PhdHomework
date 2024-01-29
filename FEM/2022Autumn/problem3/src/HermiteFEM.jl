# author: bcynuaa
# date: 2023/01/17

using LinearAlgebra;
using SparseArrays;
include("Case.jl");

function generate_K_b(case::Case)::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    K::SparseMatrixCSC{Float64, Int64} = spzeros(Float64, 2*case.N, 2*case.N);
    b::Vector{Float64} = zeros(Float64, 2*case.N);
    h::Float64 = case.h;
    k_part::Matrix{Float64} = [
        12/h^3 6/h^2 -12/h^3 6/h^2;
        6/h^2 4/h -6/h^2 2/h;
        -12/h^3 -6/h^2 12/h^3 -6/h^2;
        6/h^2 2/h -6/h^2 4/h
    ];
    qb_part::Vector{Float64} = case.q * [
        h/2, h^2/12, h/2, -h^2/12
    ];
    for k = 1: case.N-1
        if k < case.middle_ID
            b[2*k-1: 2*k+2] += qb_part;
            K[2*k-1: 2*k+2, 2*k-1: 2*k+2] .+= case.Es * case.J1 * k_part;
        else
            K[2*k-1: 2*k+2, 2*k-1: 2*k+2] .+= case.Es * case.J2 * k_part;
        end
    end
    b[2*case.middle_ID-1] += case.F;
    b[end] += -case.M;
    return K, b;
end

function solve(case::Case)::Vector{Float64}
    K::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64} = generate_K_b(case);
    index::Vector{Int64} = push!(Vector(3: 2*case.N-2), 2*case.N);
    bm::Vector{Float64} = b[index];
    Km::SparseMatrixCSC{Float64, Int64} = K[index, index];
    um::Vector{Float64} = Km \ bm;
    u::Vector{Float64} = zeros(Float64, 2*case.N);
    u[index] = um;
    return u;
end

function hermite_interpolation(case::Case, u::Vector{Float64})::Function
    w::Vector{Float64} = u[1: 2: end];
    theta::Vector{Float64} = u[2: 2: end];
    h::Float64 = case.h;
    k::Int64 = 0;
    function her_in(x::Float64)::Vector{Float64}
        res::Vector{Float64} = [0., 0.];
        if x < minimum(case.x)
            nothing;
        elseif x > maximum(case.x)
            nothing;
        elseif x in case.x
            k = findfirst(item->item==x, case.x);
            res = [w[k], theta[k]];
        else
            k = findfirst(item->item>x, case.x) - 1;
            x1::Float64 = case.x[k];
            t::Float64 = x - x1;
            wtheta::Vector{Float64} = [w[k], theta[k], w[k+1], theta[k+1]];
            coeff1::Vector{Float64} = [
                (h-t)^2 * (h + 2*t), t * (h-t)^2 * h, t^2 * (3*h - 2*t), t^2 * (t-h) * h
            ] ./ h^3;
            coeff2::Vector{Float64} = [
                6 * t * (t-h), (h - 3*t) * (h - t) * h, 6 * t * (h-t), t * (3*t - 2*h) * h
            ] ./ h^3;
            res = [coeff1' * wtheta, coeff2' * wtheta];
        end
        return res;
    end
    return her_in;
end

const default_minor::Int64 = 4;
function interpolation_minor(case::Case, u::Vector{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    f::Function = hermite_interpolation(case, u);
    N_minor::Int64 = default_minor * (case.N-1) + 1;
    x_minor::Vector{Float64} = Vector( LinRange(case.x[1], case.x[end], N_minor) );
    w_minor::Vector{Float64} = zeros(Float64, N_minor);
    theta_minor::Vector{Float64} = zeros(Float64, N_minor);
    for k = 1: N_minor
        w_minor[k], theta_minor[k] = f(x_minor[k]);
    end
    return x_minor, w_minor, theta_minor;
end