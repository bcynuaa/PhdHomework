"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-12 20:30:54
 # @ license: MIT
 # @ description: contains the projection operation
 """

using SparseArrays;

include("Case.jl");

"""
flag:
- flag=1: u(x+dx)-u(x)
- flag=2: u(x)-u(x-dx)
- flag=3: (u(x+dx) - u(x-dx))/2
"""
function partialX(j::Int64, i::Int64, case::Case, u::Matrix{Float64}, flag::Int64)::Float64
    if flag==1
        return (u[j, i+1] - u[j, i])/case.dx_;
    elseif flag==2
        return (u[j, i] - u[j, i-1])/case.dx_;
    else
        return (u[j, i+1] - u[j, i-1])/2/case.dx_;
    end
end

"""
flag:
- flag=1: v(y+dy)-v(y)
- flag=2: v(y)-v(y-dy)
- flag=3: (v(y+dy) - v(y-dy))/2
"""
function partialY(j::Int64, i::Int64, case::Case, v::Matrix{Float64}, flag::Int64)::Float64
    if flag==1
        return (v[j+1, i] - v[j, i])/case.dx_;
    elseif flag==2
        return (v[j, i] - v[j-1, i])/case.dx_;
    else
        return (v[j+1, i] - v[j-1, i])/2/case.dx_;
    end
end

function generateNabla(case::Case, u::Matrix{Float64}, v::Matrix{Float64})::Vector{Float64}
    N::Int64 = case.N_;
    dx::Float64 = case.dx_;
    nabla::Matrix{Float64} = zeros(Float64, N, N);
    for j=1:N
        for i=1:N
            if j==1 && i==1
                nabla[j, i] = partialX(j, i, case, u, 1) + partialY(j, i, case, v, 1);
            elseif j==1 && i==N
                nabla[j, i] = partialX(j, i, case, u, 2) + partialY(j, i, case, v, 1);
            elseif j==N && i==1
                nabla[j, i] = partialX(j, i, case, u, 1) + partialY(j, i, case, v, 2);
            elseif j==N && i==N
                nabla[j, i] = partialX(j, i, case, u, 2) + partialY(j, i, case, v, 2);
            elseif j==1 && 1<i<N
                nabla[j, i] = partialX(j, i, case, u, 3) + partialY(j, i, case, v, 1);
            elseif j==N && 1<i<N
                nabla[j, i] = partialX(j, i, case, u, 3) + partialY(j, i, case, v, 2);
            elseif 1<j<N && i==1
                nabla[j, i] = partialX(j, i, case, u, 1) + partialY(j, i, case, v, 3);
            elseif 1<j<N && i==N
                nabla[j, i] = partialX(j, i, case, u, 2) + partialY(j, i, case, v, 3);
            else
                nabla[j, i] = partialX(j, i, case, u, 3) + partialY(j, i, case, v, 3);
            end
        end
    end
    return reshape(transpose(nabla), N^2);
end

function generateNablaDot(case::Case, q::Matrix{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}}
    u::Matrix{Float64} = zeros(Float64, case.N_-2, case.N_-2);
    v::Matrix{Float64} = zeros(Float64, case.N_-2, case.N_-2);
    for j = 1: case.N_-2
        for i = 1: case.N_ - 2
            u[j, i] = (q[j+1, i+2] - q[j+1, i]) / 2 / case.dx_;
            v[j, i] = (q[j+2, i+1] - q[j, i+1]) / 2 / case.dx_;
        end
    end
    return u, v;
end

function projection!(
    case::Case,
    w3u::Matrix{Float64}, w3v::Matrix{Float64},
    matrixA_projection::SparseMatrixCSC{Float64, Int64},
)::Int64
    nabla_w2::Vector{Float64} = generateNabla(case, w3u, w3v);
    q::Vector{Float64} = matrixA_projection \ nabla_w2;
    nabla_dot_q::Tuple{Matrix{Float64}, Matrix{Float64}} = generateNablaDot(case, Matrix(
        transpose(reshape(q, case.N_, case.N_))
    ));
    w3u[2: case.N_ - 1, 2: case.N_ - 1] .-= nabla_dot_q[1];
    w3v[2: case.N_ - 1, 2: case.N_ - 1] .-= nabla_dot_q[2];
    return 0;
end