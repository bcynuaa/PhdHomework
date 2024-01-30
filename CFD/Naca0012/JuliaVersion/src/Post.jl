# author: bcynuaa
# date: 2022/10/04

module Post # module begin

include("RK4.jl");
include("Visualization.jl");
using .RK4;
using .Visualization;

default_image_relative_path::String = "plot//"; # the default relative path to save plot images

###################################################################################################

"""
# function save_Cp_CL_CD(case::Case, grid::Grid2D, P::Vector{Float64}; plot::Bool=false)::Int64

**save the Cp, CL, CD data**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `P::Vector{Float64}` -> the cell's pressure value
;
- `plot::Bool` -> whether to do the plot, default to false

return `0`;
"""
function save_Cp_CL_CD(case::Case, grid::Grid2D, P::Vector{Float64}; plot::Bool=false)::Int64
    vec_CL::Vector{Float64} = [-sin(case.theta), cos(case.theta)];
    vec_CD::Vector{Float64} = [cos(case.theta), sin(case.theta)];
    q_inf::Float64 = case.rho_inf * case.C_inf^2 * case.Ma_inf^2 / 2.;
    
    Cp_upper::Vector{Float64} = [];
    Cp_lower::Vector{Float64} = [];
    x_upper::Vector{Float64} = [];
    x_lower::Vector{Float64} = [];
    
    cell_ID::Int64 = 0;
    x_c::Float64 = 0.;
    y_c::Float64 = 0.;
    Cp_cell::Float64 = 0.;
    
    fx::Float64 = 0.;
    fy::Float64 = 0.;
    CL::Float64 = 0.;
    CD::Float64 = 0.;
    
    for k = 1: length(grid.wall_cell_ID)
        cell_ID = grid.wall_cell_ID[k];
        x_c = (grid.xy[grid.icell[cell_ID, 1], 1] + grid.xy[grid.icell[cell_ID, 2], 1] + grid.xy[grid.icell[cell_ID, 3], 1]) / 3.;
        y_c = (grid.xy[grid.icell[cell_ID, 1], 2] + grid.xy[grid.icell[cell_ID, 2], 2] + grid.xy[grid.icell[cell_ID, 3], 2]) / 3.;
        Cp_cell = (P[cell_ID] - case.P_inf) / q_inf;
        if y_c > 0
            push!(Cp_upper, Cp_cell);
            push!(x_upper, x_c);
        else
            push!(Cp_lower, Cp_cell);
            push!(x_lower, x_c);
        end
        fx, fy = grid.vec_n_wall_cell[cell_ID, :] .* P[cell_ID] * grid.ds[grid.wall_cell_edge_ID[k]];
        CL += vec_CL' * [fx, fy] / q_inf;
        CD += vec_CD' * [fx, fy] / q_inf;
    end

    file_name::String = case.data_name;
    open(case.data_path * "Cp_CL_CD_" * file_name, "w") do io
        write(io, "CL = $CL\nCD = $CD\n\n");
        
        write(io, "x_upper:\n");
        for xu in x_lower
            write(io, "$xu    ");
        end
        write(io, "\n");
        
        write(io, "Cp_upper:\n");
        for Cpu in Cp_upper
            write(io, "$Cpu    ");
        end
        write(io, "\n");

        write(io, "x_lower:\n");
        for xl in x_lower
            write(io, "$xl    ");
        end
        write(io, "\n");

        write(io, "Cp_lower:\n");
        for Cpl in Cp_lower
            write(io, "$Cpl    ");
        end
        write(io, "\n");
    end

    println("successfully save Cp CL CD data to $(case.data_path * "Cp_CL_CD_" * file_name)");

    if plot == true
        Cp_plot(Cp_upper, x_upper, Cp_lower, x_lower, CL, CD, case.data_path * default_image_relative_path);
    end
    return 0;
end

###################################################################################################

"""

# function save_Tecplot_Paraview(case::Case, grid::Grid2D, W_expand::Matrix{Float64}; plot::Bool=false)::Int64

**save the result as Tecplot format, as well as available for paraview (FEM format), 
moreover, the value is X, Y, density, pressure, Ma, U, V.**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `W_expand::Matrix{Float64}` -> the expand conservative value w of all cells
;
- `plot::Bool` -> whether to do the plot, default to false

return `0`;
"""
function save_Tecplot_Paraview(case::Case, grid::Grid2D, W_expand::Matrix{Float64}; plot::Bool=false)::Int64
    VALUE::Matrix{Float64} = zeros(Float64, grid.nnodes, 7);
    COUNT::Vector{Float64} = zeros(Float64, grid.nnodes);
    cell_ID::Int64 = 0;
    k::Int64 = 0;
    for cell_ID = 1: grid.ncells
        for k = 1: 3
            VALUE[grid.icell[cell_ID, k], :] .+= W_expand[cell_ID, :] .* grid.vol[cell_ID];
            COUNT[grid.icell[cell_ID, k]] += grid.vol[cell_ID];
        end
    end
    VALUE ./= COUNT;

    file_name::String = case.data_name;
    MA::Vector{Float64} = sqrt.(VALUE[:, 2].^2 .+ VALUE[:, 3].^2) ./ VALUE[:, 7];
    RES::Matrix{Float64} = zeros(Float64, grid.nnodes, 5);
    RES[:, 1] = VALUE[:, 1];
    RES[:, 2] = VALUE[:, 5];
    RES[:, 3] = MA;
    RES[:, 4] = VALUE[:, 2];
    RES[:, 5] = VALUE[:, 3];
    open(case.data_path * "Tec_Para_" * file_name, "w") do io
        write(io, "TITLE=\"$(case.grid_file)\"\nVARIABLES=\"X\", \"Y\", \"density\", \"pressure\", \"MA\", \"U\", \"V\"\nZONE N=$(grid.nnodes), E=$(grid.ncells), F=FEPOINT, ET=TRIANGLE\n");
        for j = 1: grid.nnodes
            write(io, "$(grid.xy[j, 1])    $(grid.xy[j, 2])    $(RES[j, 1])    $(RES[j, 2])    $(RES[j, 3])    $(RES[j, 4])   $(RES[j, 5])\n");
        end
        for j = 1: grid.ncells
            write(io, "$(grid.icell[j, 1])   $(grid.icell[j, 2])   $(grid.icell[j, 3])\n");
        end
    end

    println("successfully save contour data (FEM) to $(case.data_path * "Tec_Para_" * file_name)");
    
    if plot == true
        value_name::Vector{String} = [
            "density", "pressure", "MA", "U", "V"
        ];
        title_name::Vector{String} = [
            "\$\\rho\$", "\$P\$", "\$Ma\$", "\$U\$", "\$V\$"
        ];
        tri = plt.matplotlib.tri.Triangulation(grid.xy[:, 1], grid.xy[:, 2], grid.icell .- 1);
        grid_plot(tri, case.grid_file, case.data_path * default_image_relative_path);
        for k = 1: 5
            contour_plot(tri, RES[:, k], value_name[k], title_name[k], case.data_path * default_image_relative_path);
        end
    end
    return 0;
end

###################################################################################################

"""
# function save_residual(case::Case, residual::Matrix{Float64}; plot::Bool=false)::Int64

**save the residual data**

- `case::Case` -> the Case contains parameters
- `residual::Matrix{Float64}` -> the residual Matrix
;
- `plot::Bool` -> whether to do the plot, default to false

return `0`;
"""
function save_residual(case::Case, residual::Matrix{Float64}; plot::Bool=false)::Int64
    file_name::String = case.data_name;
    res::Vector{Float64} = residual[:, 1];
    open(case.data_path * "residual_" * file_name, "w") do io
        for k = 1: case.STEP
            write(io, "$(res[k])\n");
        end
    end

    println("successfully save residual data to $(case.data_path * "residual_" * file_name)");
    
    if plot == true
        residual_plot(res, case.data_path * default_image_relative_path);
    end
    return 0;
end

###################################################################################################

"""
# function save_data(case::Case, grid::Grid2D, W::Matrix{Float64})::Int64

**use for contour data tests, only parse W without residual**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `W::Matrix{Float64}` -> the result of conservative value W

return `0`;
"""
function save_data(case::Case, grid::Grid2D, W::Matrix{Float64})::Int64
    W_expand::Matrix{Float64} = expand_W(case, W);
    if isdir(case.data_path) == false
        mkpath(case.data_path);
    end
    save_Cp_CL_CD(case, grid, W_expand[:, 5]);
    save_Tecplot_Paraview(case, grid, W_expand);
    return 0;
end

"""
# save_data(case::Case, grid::Grid2D, W::Matrix{Float64}, residual::Matrix{Float64}; plot::Bool=false)::Int64

**use to save data and do the plot**

- Cp, CL, CD data;
- contour data of density, pressure, Ma, U, V
- residual change data

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `W::Matrix{Float64}` -> the result of conservative value W
- `residual::Matrix{Float64}` -> the residual Matrix
;
- `plot::Bool` -> whether to do the plot, default to false

return `0`;
"""
function save_data(case::Case, grid::Grid2D, W::Matrix{Float64}, residual::Matrix{Float64}; plot::Bool=false)::Int64
    W_expand::Matrix{Float64} = expand_W(case, W);
    if isdir(case.data_path) == false
        mkpath(case.data_path);
    end
    if plot == true
        mkpath(case.data_path * default_image_relative_path);
    end
    save_Cp_CL_CD(case, grid, W_expand[:, 5]; plot=plot);
    save_Tecplot_Paraview(case, grid, W_expand; plot=plot);
    save_residual(case, residual; plot=plot);
    return 0;
end

###################################################################################################

export save_data;
export runge_kutta4;
export Case, Grid2D, expand_w, expand_W;

end # module end