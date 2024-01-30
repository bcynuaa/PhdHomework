# author: bcynuaa
# date: 2023/01/05

module Visualization # module end

using PyCall;
import PyPlot as plt;

const precision::Int64 = 2; # default data precision
const dpi::Int64 = 200; # default dpi
const label_fontsize::Float64 = 20.; # default label and text fontsize
const title_fontsize::Float64 = 25.;  # default title fontsize

###################################################################################################

const figsize_contour::Tuple{Int64, Int64} = (10, 10); # default contour figsize
const xlim::Vector{Float64} = [6., 10.]; # default xlim
const ylim::Vector{Float64} = [0., 4.]; # default ylim
const delta::Float64 = 0.1; # default delta value of colorbar
const cmap::String = "coolwarm"; # default cmap
const color_layer::Int64 = 15; # default cmap layer number
const edge_color::String = "k"; # default contour line color
const lw::Float64 = 0.1; # default contour line width
const mesh_color::String = "white"; # default mesh(grid) line color
const mesh_lw::Float64 = 0.15; # default mesh(grid) line width
const shrink::Float64 = 0.8; # the shrink scale of colorbar

"""
# function contour_plot(triangle, value::Vector{Float64}, value_name::String, title_name::String, save_path::String)::Int64

**do the contour plot**

- `triangle::PyCall.PyObject` -> the Triangulation from matplotlib.tri.Triangulation is a python object
- `value::Vector{Float64}` -> the value to do the contour plot, which size should be the same as tiangle's nodes
- `value_name::String` -> the value name like density, pressure, Ma, U, V
- `title_name::String` -> the title name usually a LaTex string like \$\\rho\$
- `save_path::String` -> the path to save the image

return `0`;
"""
function contour_plot(
    triangle,
    value::Vector{Float64},
    value_name::String,
    title_name::String,
    save_path::String
)::Int64
    # vmin::Float64 = minimum(value);
    # vmax::Float64 = maximum(value);
    # vmin_cb::Float64 = floor(vmin, digits=precision);
    # # vmax_cb::Float64 = ceil(vmax, digits=precision);
    # vmax_cb::Float64 = floor(vmax, digits=precision) + 10. ^ (-precision);
    # contour_level::Vector{Float64} = Vector(vmin_cb: delta: vmax_cb);
    plt.figure(figsize=figsize_contour, facecolor="white", dpi=dpi);
    plt.gca().set_aspect(1);
    # plt.tricontour(triangle, value, colors=edge_color, levels=contour_level, linewidths=lw);
    # plt.tricontourf(triangle, value, cmap=cmap, levels=contour_level, vmin=vmin, vmax=vmax);
    plt.tricontour(triangle, value, color_layer, colors=edge_color, linewidths=lw);
    plt.tricontourf(triangle, value, color_layer, cmap=cmap);
    # plt.triplot(triangle, color=mesh_color, lw=mesh_lw);
    plt.colorbar(shrink=shrink);
    plt.xlim(xlim[1], xlim[2]);
    plt.ylim(ylim[1], ylim[2]);
    plt.xlabel("\$X\$", fontsize=label_fontsize);
    plt.ylabel("\$Y\$", fontsize=label_fontsize);
    plt.title("Contour of " * title_name, fontsize=title_fontsize);
    plt.savefig(save_path * "contour_" * value_name * ".png", bbox_inches="tight");
    plt.savefig(save_path * "contour_" * value_name * ".pdf", bbox_inches="tight");
    plt.show();
    println("successfully save figure $value_name !");
    return 0;
end

###################################################################################################

const color_grid::String = "r"; # default grid mesh line color
const lw_grid::Float64 = 0.3; # default grid mesh line width

"""
# function grid_plot(triangle::PyCall.PyObject, title_name::String, save_path::String)

**do the grid plot**

- `triangle::PyCall.PyObject` -> the Triangulation from matplotlib.tri.Triangulation is a python object
- `value::Vector{Float64}` -> the value to do the contour plot, which size should be the same as tiangle's nodes
- `value_name::String` -> the value name like density, pressure, Ma, U, V
- `title_name::String` -> the title name usually a LaTex string like \$\\rho\$
- `save_path::String` -> the path to save the image

return `0`;
"""
function grid_plot(triangle::PyCall.PyObject, title_name::String, save_path::String)
    plt.figure(figsize=figsize_contour, facecolor="white", dpi=dpi);
    plt.gca().set_aspect(1);
    plt.triplot(triangle, color=color_grid, lw=lw_grid);
    plt.xlim(xlim[1], xlim[2]);
    plt.ylim(ylim[1], ylim[2]);
    plt.xlabel("\$X\$", fontsize=label_fontsize);
    plt.ylabel("\$Y\$", fontsize=label_fontsize);
    plt.title(title_name, fontsize=title_fontsize);
    plt.savefig(save_path * "grid.png", bbox_inches="tight");
    plt.savefig(save_path * "grid.pdf", bbox_inches="tight");
    plt.show();
    println("successfully save figure grid !");
end

###################################################################################################

const figsize_cp::Tuple{Int64, Int64} = (10, 6); # default cp figsize
const upper_color::String = "blue"; # default upper scatter color
const upper_marker::String = "."; # default upper scatter marker
const lower_color::String = "orange"; # default lower scatter color
const lower_marker::String = "v"; # default lower scatter marker
const xticks_cp::Vector{Float64} = Vector(0.: 0.1: 1.); # default cp plot xticks
const label_pos::String = "lower right"; # where the label lie
const text_pos::Vector{Float64} = [0.5, 0.8]; # where the text (CL, CD) lie
const text_color::String = "k"; # default text color
const bbox::Dict{String, String} = Dict("facecolor"=>"yellow", "boxstyle"=>"round"); # default text style

"""
# function Cp_plot(Cp_upper::Vector{Float64}, x_upper::Vector{Float64}, Cp_lower::Vector{Float64}, x_lower::Vector{Float64}, CL::Float64, CD::Float64, save_path::String)::Int64

**do the plot for Cp distribution**

- `Cp_upper::Vector{Float64}` -> the upper Cp value
- `x_upper::Vector{Float64}` -> the upper Cp's x position
- `Cp_lower::Vector{Float64}` -> the lower Cp value
- `x_lower::Vector{Float64}` -> the lower Cp's x position
- `CL::Float64` -> CL value
- `CD::Float64` -> CD value
- `save_path::String` -> the path to save the image

return `0`;
"""
function Cp_plot(
    Cp_upper::Vector{Float64},
    x_upper::Vector{Float64},
    Cp_lower::Vector{Float64},
    x_lower::Vector{Float64},
    CL::Float64,
    CD::Float64,
    save_path::String
)::Int64
    plt.figure(figsize=figsize_cp, facecolor="white", dpi=dpi);
    plt.gca().invert_yaxis();
    plt.scatter(x_upper, Cp_upper, label="upper \$C_p\$", color=upper_color, marker=upper_marker);
    plt.scatter(x_lower, Cp_lower, label="lower \$C_p\$", color=lower_color, marker=lower_marker);
    plt.legend(loc=label_pos, fontsize=label_fontsize);
    plt.grid(true);
    plt.xticks(xticks_cp);
    plt.xlabel("\$X\$", fontsize=label_fontsize);
    plt.ylabel("\$C_p\$", fontsize=label_fontsize);
    plt.title("Distribution of \$C_P\$", fontsize=title_fontsize);
    plt.text(
        text_pos[1], text_pos[2],
        "\$C_L=$(round(CL, digits=precision))\$\n\$C_D=$(round(CD, digits=precision))\$",
        fontsize=label_fontsize,
        bbox=bbox,
        color=text_color
    );
    plt.savefig(save_path * "distribution_Cp.png", bbox_inches="tight");
    plt.savefig(save_path * "distribution_Cp.pdf", bbox_inches="tight");
    plt.show()
    println("successfully save figure distribution of Cp !");
    return 0;
end

###################################################################################################

const figsize_res::Tuple{Int64, Int64} = (10, 6); # default residual figsize
const lw_res::Float64 = 3.; # default residual plot line width
const color_res::String = "blue"; # default residual plot line color
const label_pos_res::String = "upper right"; # where the residual label lie

"""
# function residual_plot(residual::Vector{Float64}, save_path::String)::Int64

**do the plot for residual**

- `residual::Vector{Float64}` -> the residual on each step
- `save_path::String` -> the path to save the image

return `0`;
"""
function residual_plot(residual::Vector{Float64}, save_path::String)::Int64
    plt.figure(figsize=figsize_res, facecolor="white", dpi=dpi);
    plt.plot(residual, lw=lw_res, color=color_res, label="residuals");
    plt.legend(loc=label_pos_res, fontsize=label_fontsize);
    plt.yscale("log");
    plt.grid(true);
    plt.xlabel("Step", fontsize=label_fontsize);
    plt.ylabel("Residual", fontsize=label_fontsize);
    plt.title("Residuals-Step", fontsize=title_fontsize);
    plt.savefig(save_path * "residuals.png", bbox_inches="tight");
    plt.savefig(save_path * "residuals.pdf", bbox_inches="tight");
    plt.show()
    println("successfully save figure distribution of residuals !");
    return 0;
end

###################################################################################################

export contour_plot, grid_plot, Cp_plot, residual_plot, plt;

end # module end