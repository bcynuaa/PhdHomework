# author: bcynuaa
# date: 2023/01/13

include("Case.jl");
include("Mesh.jl");
using SparseArrays;

function generate_K_b(case::Case, mesh::Mesh)::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}}
    K::SparseMatrixCSC{Float64, Int64} = spzeros(Float64, mesh.nnodes, mesh.nnodes);
    b::Vector{Float64} = zeros(Float64, mesh.nnodes);

    ID1::Int64 = 0;
    ID2::Int64 = 0;
    ID3::Int64 = 0;
    x1::Float64 = 0.;
    y1::Float64 = 0.;
    x2::Float64 = 0.;
    y2::Float64 = 0.;
    x3::Float64 = 0.;
    y3::Float64 = 0.;
    A2::Float64 = 0.;
    k11::Float64 = 0.;
    k12::Float64 = 0.;
    k13::Float64 = 0.;
    k22::Float64 = 0.;
    k23::Float64 = 0.;
    k33::Float64 = 0.;
    b1::Float64 = 0.;
    b2::Float64 = 0.;
    b3::Float64 = 0.;
    
    for k = 1: mesh.ncells
        ID1, ID2, ID3 = mesh.icell[k, :];
        x1, y1 = mesh.xy[ID1, :];
        x2, y2 = mesh.xy[ID2, :];
        x3, y3 = mesh.xy[ID3, :];
        A2 = 2* det([
            1 x1 y1;
            1 x2 y2;
            1 x3 y3
        ]);
        k11 = (x2 - x3)^2 + (y2 - y3)^2;
        k12 = (x1 - x3) * (-x2 + x3) + (y1 - y3) * (-y2 + y3);
        k13 = (x1 - x2) * (x2 - x3) + (y1 - y2) * (y2 - y3);
        k22 = (x1 - x3)^2 + (y1 - y3)^2;
        k23 = -(x1^2 + x2 * x3 - x1 * (x2 + x3) + (y1 - y2) * (y1 - y3));
        k33 = (x1 - x2)^2 + (y1 - y2)^2;
        b1 = (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3)) / 6;
        b2 = (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1+y3)) / 6;
        b3 = (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1+y3)) / 6;

        K[ID1, ID1] += k11 / A2;
        K[ID1, ID2] += k12 / A2;
        K[ID1, ID3] += k13 / A2;

        K[ID2, ID1] += k12 / A2;
        K[ID2, ID2] += k22 / A2;
        K[ID2, ID3] += k23 / A2;

        K[ID3, ID1] += k13 / A2;
        K[ID3, ID2] += k23 / A2;
        K[ID3, ID3] += k33 / A2;

        b[ID1] += b1;
        b[ID2] += b2;
        b[ID3] += b3;
    end

    return K, case.a * b;
end

function boundary_node_ID(case::Case, mesh::Mesh)::Vector{Int64}
    node_ID_boundary::Vector{Int64} = [];
    for j = 1: mesh.nnodes
        if whether_on_boundary(case, mesh.xy[j, :]) == true
            push!(node_ID_boundary, j);
        end
    end
    return node_ID_boundary;
end

function solve(case::Case, mesh::Mesh)::Vector{Float64}
    K::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64} = generate_K_b(case, mesh);
    node_ID_boundary::Vector{Int64} = boundary_node_ID(case, mesh);
    node_ID_inner::Vector{Int64} = setdiff(Vector(1: mesh.nnodes), node_ID_boundary);
    n_inner::Int64 = length(node_ID_inner);
    b_inner::Vector{Float64} = b[node_ID_inner];
    K_inner::SparseMatrixCSC{Float64, Int64} = spzeros(Float64, n_inner, n_inner);
    for j = 1: n_inner
        for i = 1: n_inner
            K_inner[j, i] = K[node_ID_inner[j], node_ID_inner[i]];
            b_inner[j] -= case.boundary_value * K_inner[j, i];
        end
    end
    u_inner::Vector{Float64} = K_inner \ b_inner;
    u::Vector{Float64} = zeros(Float64, mesh.nnodes);
    u[node_ID_boundary] .= case.boundary_value;
    u[node_ID_inner] .= u_inner;
    return u;
end