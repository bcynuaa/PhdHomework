# author: bcynuaa
# date: 2023/01/05

include("..//head.jl");
include("..//..//src//EulerEq.jl");
using .EulerEq;

Ma = 4.;

data_name::String = "MA$Ma.dat";

case::Case = Case(Ma, theta; grid_file=grid_file, data_path=data_path, data_name=data_name, STEP=STEP);
grid::Grid2D = Grid2D(case.grid_file);

solve(case, grid; plot=plot);