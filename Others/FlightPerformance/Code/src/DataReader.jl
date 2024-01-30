# author: bcynuaa
# date: 2022/10/14

module DataReader # begin module

using Parsers;
using XMLDict;
# using ExportAll;

mytag::String = "performance"
split_char::Vector{Char} = [' ', '\n'];

struct Table
    name::String;
    x_name::String;
    y_name::String;
    x::Vector{Float64};
    y::Vector{Float64};
    data::Matrix{Float64};
    row::Int64;
    col::Int64;
end

struct Table3D
    name::String;
    x_name::String;
    y_name::String;
    z_name::String;
    x::Vector{Float64};
    y::Vector{Float64};
    z::Vector{Float64};
    data::Array{Float64, 3};
    row::Int64;
    col::Int64;
    num::Int64;
end

function Table(
    name::String,
    x_name::String,
    y_name::String,
    x::Vector{Float64},
    y::Vector{Float64},
    data::Matrix{Float64}
)::Table
    row::Int64 = length(y);
    col::Int64 = length(x);
    return Table(
        name, x_name, y_name,
        x, y, data,
        row, col
    );
end

function Table3D(
    name::String,
    x_name::String,
    y_name::String,
    z_name::String,
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    data::Array{Float64, 3}
)::Table3D
    row::Int64 = length(y);
    col::Int64 = length(x);
    num::Int64 = length(z);
    return Table3D(
        name,
        x_name, y_name, z_name,
        x, y, z, data,
        row, col, num
    )
end

struct OriginData
    file_name::String;
    mass_kg::Float64;
    CL_base::Table;
    CD_base::Table;
    thrustAB::Table;
    thrust::Table;
    cft::Table3D;
    cftAB::Table;
end

function OriginData(file_name::String)
    dict = xml_dict(String(read(file_name));)[mytag];
    mass_kg::Float64 = parse(Float64, dict["property"][""]);
    tables::Vector{Table} = Vector{Table}(undef, 0);
    for i = 1:length(dict["table"])
        if i == 5
            continue;
        end
        x_name = dict["table"][i]["breakpoints"][1][:property];
        y_name = dict["table"][i]["breakpoints"][2][:property];
        x = parse.(Float64, split(dict["table"][i]["breakpoints"][1][""], split_char; keepempty=false));
        y = parse.(Float64, split(dict["table"][i]["breakpoints"][2][""], split_char; keepempty=false));
        data = Matrix(
            transpose(
                reshape(
                    parse.(
                        Float64,
                        split(dict["table"][i]["tableData"], split_char; keepempty=false)
                    ),
                    (length(x), length(y))
                )
            )
        );
        push!(
            tables,
            Table(
                dict["table"][i][:name],
                x_name,
                y_name,
                x,
                y,
                data
            );
        );
    end
    x_name = dict["table"][5]["breakpoints"][1][:property];
    y_name = dict["table"][5]["breakpoints"][2][:property];
    z_name = dict["table"][5]["breakpoints"][3][:property];
    x = parse.(Float64, split(dict["table"][5]["breakpoints"][1][""], split_char; keepempty=false));
    y = parse.(Float64, split(dict["table"][5]["breakpoints"][2][""], split_char; keepempty=false));
    z = parse.(Float64, split(dict["table"][5]["breakpoints"][3][""], split_char; keepempty=false));
    data = permutedims(
        reshape(
            parse.(
                Float64,
                split(dict["table"][5]["tableData"], split_char; keepempty=false)
            ),
            (length(x), length(y), length(z))
        ),
        [3, 2, 1]
    );
    table5::Table3D = Table3D(
        dict["table"][5][:name],
        x_name,
        y_name,
        z_name,
        x,
        y,
        z,
        data
    );
    
    return OriginData(
        file_name,
        mass_kg,
        tables[1],
        tables[2],
        tables[3],
        tables[4],
        table5,
        tables[5]
    );
end

export Table, OriginData;

end # end module