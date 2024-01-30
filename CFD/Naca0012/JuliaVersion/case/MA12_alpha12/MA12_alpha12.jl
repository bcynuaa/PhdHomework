# author: bcynuaa
# date: 2022/11/12

include("..//head.jl");
include("..//..//src//EulerEq.jl");
using .EulerEq;

Ma = 1.2;
theta = 1.2;
data_name::String = generate_data_name(Ma, theta);
theta *= pi/180.;

case::Case = Case(Ma, theta; grid_file=grid_file, data_path=data_path, data_name=data_name, STEP=STEP);
grid::Grid2D = Grid2D(case.grid_file);

solve(case, grid; plot=plot);