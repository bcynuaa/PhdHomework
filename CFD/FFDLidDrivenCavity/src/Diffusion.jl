"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-12 20:13:15
 # @ license: MIT
 # @ description: contains the diffusion operation
 """

using SparseArrays;

include("Case.jl");

function diffusionU!(
    case::Case, w2u::Matrix{Float64}, matrixA::SparseMatrixCSC{Float64, Int64}
)::Int64
    u_cur::Vector{Float64} = reshape(
        transpose(w2u[2: case.N_ - 1, 2: case.N_ - 1]), (case.N_ - 2)^2);
    u_cur[end-N+1: end] .+= case.lambda_ * case.u_upper_;
    u_next::Vector{Float64} = matrixA \ u_cur;
    w2u[2: case.N_ - 1, 2: case.N_ - 1] .= transpose(
        reshape(u_next, case.N_ - 2, case.N_ - 2));
    return 0;
end

function diffusionV!(
    case::Case, w2v::Matrix{Float64}, matrixA::SparseMatrixCSC{Float64, Int64}
)::Int64
    v_cur::Vector{Float64} = reshape(
        transpose(w2v[2: case.N_ - 1, 2: case.N_ - 1]), (case.N_ - 2)^2);
    v_next::Vector{Float64} = matrixA \ v_cur;
    w2v[2: case.N_ - 1, 2: case.N_ - 1] .= transpose(
        reshape(v_next, case.N_ - 2, case.N_ - 2));
    return 0;
end

function diffusion!(
    case::Case, w2u::Matrix{Float64}, w2v::Matrix{Float64},
    matrixA::SparseMatrixCSC{Float64, Int64}
)::Int64
    diffusionU!(case, w2u, matrixA);
    diffusionV!(case, w2v, matrixA);
    return 0;
end