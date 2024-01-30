"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-12 19:27:35
 # @ license: MIT
 # @ description: contains the sparse matrix and related function
 """

using SparseArrays;
include("Case.jl");

function addToSparseMatrix!(
    J::Vector{Int64},
    I::Vector{Int64},
    V::Vector{Float64},
    index::Int64,
    m::Int64,
    n::Int64,
    v::Float64,
)::Int64
    J[index] = m;
    I[index] = n;
    V[index] = v;
    return index + 1;
end

function addByPoint!(
    J::Vector{Int64},
    I::Vector{Int64},
    V::Vector{Float64},
    index::Int64,
    i::Int64,
    diagv::Float64,
    subv::Float64,
    point::Vector{Int64},
)::Int64
    index = addToSparseMatrix!(J, I, V, index, i, i, diagv);
    for j in point
        index = addToSparseMatrix!(J, I, V, index, i, i+j, subv);
    end
    return index;
end

function makeSparseMatrix(N::Int64, diagv::Float64, subv::Float64)::SparseMatrixCSC{Float64, Int64}
    NA::Int64 = 5*(N-2)^2 + 4*4*(N-2) + 3*4;
    index::Int64 = 1;
    J::Vector{Int64} = zeros(Int64, NA);
    I::Vector{Int64} = zeros(Int64, NA);
    V::Vector{Float64} = zeros(Float64, NA);
    for i = 1:N^2
        if i == 1
            point = [1, N];
        elseif i == N
            point = [-1, N];
        elseif i == (N-1)*N+1
            point = [1, -N];
        elseif i == N^2
            point = [-1, -N]
        elseif 1 < i < N
            point = [-1, 1, N]
        elseif (N-1)*N < i < N^2
            point = [-1, 1, -N]
        elseif i % N == 1
            point = [1, -N, N]
        elseif i % N == 0
            point = [-1, -N, N]
        else
            point = [-1, 1, -N, N]
        end
        index = addByPoint!(J, I, V, index, i, diagv, subv, point);
    end
    return sparse(J, I, V, N^2, N^2);
end

function makeMatrixA(case::Case)::SparseMatrixCSC{Float64, Int64}
    N::Int64 = case.N_;
    lambda::Float64 = case.lambda_;
    diagv::Float64 = 1 + 4 * lambda;
    subv::Float64 = -lambda;
    return makeSparseMatrix(N-2, diagv, subv);
end

function makeMatrixAProjection(case::Case)::SparseMatrixCSC{Float64, Int64}
    N::Int64 = case.N_;
    dx::Float64 = case.dx_;
    diagv::Float64 = -4. / dx^2;
    subv::Float64 = 1. / dx^2;
    return makeSparseMatrix(N, diagv, subv);
end