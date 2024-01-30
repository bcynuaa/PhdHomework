# author: bcynuaa
# date: 2022/09/27

module Grid # module begin

using ExportAll;
using SparseArrays;
using LinearAlgebra;

###################################################################################################

"""
# readfile_Vector_String(filename::String, skiprow::Int64=0, row::Int64=1)::Vector{String}

**To read a file by lines and shift it to a string Vector**

- `filename::String` -> the path of the file
- `skiprow::Int64` -> the rows to skip
- `row::Int64` -> how much row you need in the file
"""
function readfile_Vector_String(filename::String, skiprow::Int64=0, row::Int64=1)::Vector{String}
    file::IOStream = open(filename, "r";)
    data::Vector{String} = readlines(file)[skiprow+1: skiprow+row];
    close(file);
    return data;
end

"""
# readVector(filename::String, T::DataType=Float64, skiprow::Int64=0, row::Int64=1)::Vector{T}

**To read a file by lines and shift it to a string Vector**

- `filename::String` -> the path of the file
- `T::DataType` -> the Vector data type
- `skiprow::Int64` -> the rows to skip
- `row::Int64` -> how much row you need in the file
"""
function readVector(filename::String, T::DataType=Float64, skiprow::Int64=0, row::Int64=1)::Vector{T}
    data::Vector{String} = readfile_Vector_String(filename, skiprow, row);
    return parse.(T, data);
end

"""
# `readMatrix(filename::String, T::DataType=Float64, skiprow::Int64=0, row::Int64=1)::Matrix{T}`

**To read a file by lines and shift it to a string Vector**

- `filename::String` -> the path of the file
- `T::DataType` -> the Matrix data type
- `skiprow::Int64` -> the rows to skip
- `row::Int64` -> how much row you need in the file
"""
function readMatrix(filename::String, T::DataType=Float64, skiprow::Int64=0, row::Int64=1)::Matrix{T}
    data::Vector{String} = readfile_Vector_String(filename, skiprow, row);
    col::Int64 = length(parse.(T, split(data[1])));
    mat::Matrix{T} = zeros(T, row, col);
    for j=1:row
        mat[j, :] = parse.(T, split(data[j]));
    end
    return mat;
end

###################################################################################################

"""
# `Grid2D::DataType`

- `nnodes::Int64` -> number of node
- `nedges::Int64` -> number of edge
- `ncells::Int64` -> number of cell
- `xy::Matrix{Float64}` -> nnodes*2's matrix for x and y
- `iedge::Matrix{Int64}` -> nedges*4's matrix for each edge, idege[n1, n2, c1, c2] while 1 for wall bc and 2 for farfield bc
- `icell::Matrix{Int64}` -> ncells*3's matrix for each cell
- `vol::Vector{Float64}` -> volume
- `dx::Vector{Float64}` -> the dx of edge
- `dy::Vector{Float64}` -> the dy of edge
- `ds::Vector{Float64}` -> the ds of edge
- `vec_n::Vector{Float64}` -> the normal verticle direction of the edge
- `vec_t::Vector{Float64}` -> the direction from i to j
- `wall_cell_ID::Vector{Int64}` -> cell's ID by the wall
- `vec_n_wall_cell::SparseMatrixCSC{Float64, Int64}` -> cell's normal direction by the wall
- `vec_t_wall_cell::SparseMatrixCSC{Float64, Int64}` -> cell's edge direction by the wall
"""
struct Grid2D
    nnodes::Int64 # number of node
    nedges::Int64 # number of edge
    ncells::Int64 # number of cell
    xy::Matrix{Float64} # nnodes*2's matrix for x and y
    iedge::Matrix{Int64} # nedges*4's matrix for each edge, idege[n1, n2, c1, c2] while -1 for wall bc and -2 for farfield bc
    icell::Matrix{Int64} # ncells*3's matrix for each cell
    vol::Vector{Float64} # each cell's area
    dx::Vector{Float64}; # each edge's dx
    dy::Vector{Float64}; # each edge's dy
    ds::Vector{Float64}; # each edge's ds
    vec_n::Matrix{Float64}; # each edge's outer normal unit vector
    vec_t::Matrix{Float64}; # each edge's unit vector
    wall_cell_edge_ID::Vector{Int64}; # cell's ID by the wall
    wall_cell_ID::Vector{Int64}; # cell's ID by the wall
    vec_n_wall_cell::SparseMatrixCSC{Float64, Int64}; # cell's normal direction by the wall
    vec_t_wall_cell::SparseMatrixCSC{Float64, Int64}; # cell's edge direction by the wall
end

"""
# `Grid2D(filename::String)::Grid2D`

**To generate a 2D grid struct Grid2D**

- `filename::String` -> the path of the grid file
"""
function Grid2D(filename::String)::Grid2D
    file = open(filename, "r");
    info::Vector{Int64} = parse.(Int64, split(readline(file)));
    close(file);
    nnodes::Int64 = info[1];
    nedges::Int64 = info[2];
    ncells::Int64 = info[3];
    xy::Matrix{Float64} = readMatrix(filename, Float64, 1, nnodes);
    iedge::Matrix{Int64} = readMatrix(filename, Int64, 1+nnodes, nedges);
    icell::Matrix{Int64} = readMatrix(filename, Int64, 1+nnodes+nedges, ncells);
    # vol::Vector{Float64} = readVector(filename, Float64, 1+nnodes+nedges+ncells, ncells);

    vol::Vector{Float64} = zeros(Float64, ncells);
    for k = 1: ncells
        vol[k] = det([
            1. xy[icell[k, 1], 1] xy[icell[k, 1], 2];
            1. xy[icell[k, 2], 1] xy[icell[k, 2], 2];
            1. xy[icell[k, 3], 1] xy[icell[k, 3], 2]
        ]) / 2.;
    end

    dx::Vector{Float64} = zeros(Float64, nedges);
    dy::Vector{Float64} = zeros(Float64, nedges);
    i::Int64 = 0;
    j::Int64 = 0
    for k = 1: nedges
        i = iedge[k, 1];
        j = iedge[k, 2];
        dx[k] = xy[j, 1] - xy[i, 1];
        dy[k] = xy[j, 2] - xy[i, 2];
    end
    ds::Vector{Float64} = sqrt.(dx.^2 .+ dy.^2);
    vec_n::Matrix{Float64} = hcat(dy, -dx) ./ ds;
    vec_t::Matrix{Float64} = hcat(dx, dy) ./ ds;
    vec_n_wall_cell::SparseMatrixCSC{Float64, Int64} = spzeros(Float64, ncells, 2);
    vec_t_wall_cell::SparseMatrixCSC{Float64, Int64} = spzeros(Float64, ncells, 2);
    wall_cell_edge_ID::Vector{Int64} = findall(item->item==-1, iedge[:, 4]);
    wall_cell_ID::Vector{Int64} = iedge[wall_cell_edge_ID, 3];
    vec_n_wall_cell[wall_cell_ID, :] = vec_n[wall_cell_edge_ID, :];
    vec_t_wall_cell[wall_cell_ID, :] = vec_t[wall_cell_edge_ID, :];
    return Grid2D(
        nnodes, nedges, ncells, xy, iedge, icell, vol,
        dx, dy, ds, vec_n, vec_t,
        wall_cell_edge_ID,
        wall_cell_ID,
        vec_n_wall_cell, vec_t_wall_cell
    );
end

###################################################################################################

"""
# grid_for_plot(grid::Grid2D)::Tuple{Matrix{Float64},Matrix{Float64}}

**To generate a x and y for the grid plot**

- `grid::Grid.Grid2D` -> the grid struct for plot
"""
function grid_for_plot(grid::Grid2D)::Tuple{Matrix{Float64},Matrix{Float64}}
    x::Matrix{Float64} = zeros(Float64, 2, grid.nedges);
    y::Matrix{Float64} = zeros(Float64, 2, grid.nedges);
    for j=1:grid.nedges
        for i=1:2
            x[i, j] = grid.xy[grid.iedge[j, i], 1];
            y[i, j] = grid.xy[grid.iedge[j, i], 2];
        end
    end
    return x, y;
end

###################################################################################################

@exportAll;

end # module end