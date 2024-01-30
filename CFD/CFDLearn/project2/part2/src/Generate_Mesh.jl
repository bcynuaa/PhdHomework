# author: bcynuaa
# date: 2023/01/05

###################################################################################################

using Parsers;
using LinearAlgebra;
using DelimitedFiles;

###################################################################################################

const height::Float64 = 4.;
const width::Float64 = 10.;
const x_inclined::Float64 = 7.;
const y_inclined::Float64 = 1.;
const error::Float64 = 1e-5;
const file_name::String = "..//data//mesh.ply";
const output_file_name::String = "..//data//mesh.grd";
const deli::String = "      ";

###################################################################################################

function readfile_Vector_String(filename::String, skiprow::Int64=0, row::Int64=1)::Vector{String}
    file::IOStream = open(filename, "r");
    data::Vector{String} = readlines(file)[skiprow+1: skiprow+row];
    close(file);
    return data;
end

function readVector(filename::String, T::DataType=Float64, skiprow::Int64=0, row::Int64=1)::Vector{T}
    data::Vector{String} = readfile_Vector_String(filename, skiprow, row);
    return parse.(T, data);
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

###################################################################################################

function vol(xy::Matrix{Float64}, icell::Matrix{Int64}, ID::Int64)
    x1::Float64, y1::Float64 = xy[icell[ID, 1], :];
    x2::Float64, y2::Float64 = xy[icell[ID, 2], :];
    x3::Float64, y3::Float64 = xy[icell[ID, 3], :];
    return det([
        1 x1 y1;
        1 x2 y2;
        1 x3 y3
    ]);
end

###################################################################################################

function on_the_incline(x::Float64, y::Float64)::Bool
    y_std::Float64 = y_inclined / (width - x_inclined) * (x - x_inclined);
    return isapprox(y_std, y; rtol = error);
end

function configure_mesh_ply(file_name::String)::Int64
    nnodes::Int64, ncells::Int64, xy::Matrix{Float64}, icell::Matrix{Int64} = read_mesh_ply(file_name);
    nedges::Int64 = nnodes + ncells - 1;
    iedge::Matrix{Int64} = zeros(Int64, nedges, 4);
    count::Int64 = 1;
    tmp_ID_list::Vector{Int64} = [];
    for k = 1: ncells
        tmp_ID_list = icell[k, :];
        push!(tmp_ID_list, tmp_ID_list[1]);
        for ID = 1: 3
            i = tmp_ID_list[ID];
            j = tmp_ID_list[ID+1];
            xi, yi = xy[i, :];
            xj, yj = xy[j, :];
            p = -1;
            if yi == 0. && yj == 0.
                if xj > xi
                    nothing;
                else
                    i, j = j, i;
                end
                p = -2;
            elseif on_the_incline(xi, yi) == true && on_the_incline(xj, yj) == true
                if xj > xi
                    nothing;
                else
                    i, j = j, i;
                end
                p = -1;
            elseif xi == width && xj == width
                if yj > yi
                    nothing;
                else
                    i, j = j, i
                end
                p = -2;
            elseif yi == height && yj == height
                if xj < xi
                    nothing;
                else
                    i, j = j, i;
                end
                p = -2;
            elseif xi == 0. && xj == 0.
                if yj < yi
                    nothing;
                else
                    i, j = j, i;
                end
                p = -2;
            else
                set = collect( intersect( Set( findall(item->item==i, iedge[:, 1]) ), Set( findall(item->item==j, iedge[:, 2]) ) ) );
                if isempty(set) == false
                    iedge[set[1], 4] = k;
                    continue;
                else
                    set = collect( intersect( Set( findall(item->item==j, iedge[:, 1]) ), Set( findall(item->item==i, iedge[:, 2]) ) ) );
                    if isempty(set) == false
                        iedge[set[1], 4] = k;
                        continue;
                    end
                end
            end
            iedge[count, :] = [i, j, k, p];
            count += 1;
        end
    end
    area::Vector{Float64} = zeros(Float64, ncells);
    for j = 1: ncells
        area[j] = vol(xy, icell, j);
    end
    open(output_file_name, "w") do io
        write(io, "$nnodes  $nedges  $ncells\n");
        writedlm(io, xy, deli);
        writedlm(io, iedge, deli);
        writedlm(io, icell, deli);
        writedlm(io, area, deli);
    end;
    return 0;
end

###################################################################################################

configure_mesh_ply(file_name);