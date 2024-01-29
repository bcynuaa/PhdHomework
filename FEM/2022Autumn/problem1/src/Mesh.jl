# author: bcynuaa
# date: 2023/01/12

###################################################################################################

using Parsers;
using LinearAlgebra;

function readfile_Vector_String(filename::String, skiprow::Int64=0, row::Int64=1)::Vector{String}
    file::IOStream = open(filename, "r");
    data::Vector{String} = readlines(file)[skiprow+1: skiprow+row];
    close(file);
    return data;
end

function readMatrix(filename::String, T::DataType=Float64, skiprow::Int64=0, row::Int64=1)::Matrix{T}
    data::Vector{String} = readfile_Vector_String(filename, skiprow, row);
    col::Int64 = length(parse.(T, split(data[1])));
    mat::Matrix{T} = zeros(T, row, col);
    for j=1:row
        mat[j, :] = parse.(T, split(data[j]));
    end
    return mat;
end

function read_mesh_ply(file_name::String)::Tuple{Int64, Int64, Matrix{Float64}, Matrix{Int64}}
    file::IOStream = open(file_name, "r");
    info::Vector{String} = readlines(file)[[4, 8]];
    nnodes::Int64 = parse(Int64, match(r"[0-9]+", info[1]).match);
    ncells::Int64 = parse(Int64, match(r"[0-9]+", info[2]).match);
    xy::Matrix{Float64} = readMatrix(file_name, Float64, 10, nnodes)[:, 1: 2];
    icell::Matrix{Int64} = readMatrix(file_name, Int64, 10+nnodes, ncells)[:, 2: 4] .+ 1;
    return nnodes, ncells, xy, icell;
end

struct Mesh
    nnodes::Int64
    ncells::Int64
    xy::Matrix{Float64}
    icell::Matrix{Int64}
end

function Mesh(file_name::String)::Mesh
    nnodes::Int64, ncells::Int64, xy::Matrix{Float64}, icell::Matrix{Int64} = read_mesh_ply(file_name);
    return Mesh(
        nnodes,
        ncells,
        xy,
        icell
    );
end