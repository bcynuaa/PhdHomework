# author: bcynuaa
# date: 2023/01/13

include("Solver.jl");
import PyPlot as plt;

mesh_file::String = "..//data//mesh.ply"; # define the grid file name
mesh::Mesh = Mesh(mesh_file); # construct the mesh file
case::Case = Case(); # construct a default case
contour_figure::String = "..//image//solution_julia."; # contour figure's path
mesh_figure::String = "..//image//mesh_figure."; # mesh figure's path
contour_layer::Int64 = 20; # how much contour line be included

triangle = plt.matplotlib.tri.Triangulation(mesh.xy[:, 1], mesh.xy[:, 2], mesh.icell .- 1); # construct a triangle object for plot

@time u::Vector{Float64} = solve(case, mesh); # solution under this mesh and case

# plot the contour of u
plt.figure(figsize=(12, 8), facecolor="white", dpi=200);
plt.gca().set_aspect(1);
plt.tricontourf(triangle, u, contour_layer, cmap="jet");
plt.colorbar();
plt.xlabel("\$x\$", fontsize=20);
plt.ylabel("\$y\$", fontsize=20);
plt.title("Contour of \$u\$", fontsize=20);
plt.savefig(contour_figure * "png", bbox_inches="tight");
plt.savefig(contour_figure * "pdf", bbox_inches="tight");
plt.show();

plt.figure(figsize=(12, 8), facecolor="white", dpi=200);
plt.gca().set_aspect(1);
plt.triplot(triangle, color="k", lw=0.1);
plt.xlabel("\$x\$", fontsize=20);
plt.ylabel("\$y\$", fontsize=20);
plt.title("Triangle Mesh", fontsize=20);
plt.savefig(mesh_figure * "png", bbox_inches="tight");
plt.savefig(mesh_figure * "pdf", bbox_inches="tight");
plt.show();