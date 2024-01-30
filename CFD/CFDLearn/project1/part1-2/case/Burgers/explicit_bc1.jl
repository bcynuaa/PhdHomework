# auther: bcynuaa
# date: 2022/10/19

include("..//..//src//Burgers1D.jl");
using .Burgers1D;
using CSV;
using DataFrames;

c::Vector{Float64} = [0.1, 0.5, 1., 5., 10.];
func1(x::Float64) = -1<x<1 ? 1-abs(x) : 0.;
func2(x::Float64) = x>0 ? 1. : 0.;
func::Array{Function} = [func1, func2];
c_name::Vector{String} = "nu=" .* string.(c);
func_name::Vector{String} = ["Triangle", "Step"];
save_path::String = "..//..//data//Burgers//explicit_bc1//";

if isdir(save_path) != true
    mkpath(save_path)
end

dx::Float64 = 0.1;
dt::Float64 = 0.005;
t_end::Float64 = 5.;
x::Vector{Float64} = Vector(-10.: dx: 10.);
t::Vector{Float64} = Vector(0.: dt: t_end);

shared_grid::Grid = Grid(x, t);

func_bc(t::Float64) = 0.;

for j = 1:length(c)
    for i = 1:length(func)
        println("c=$(c[j]), $(func_name[i]) IC");
        file_name = func_name[i] * "_" * c_name[j];
        solver = Solver(shared_grid, c[j], func[i], func_bc; bc_type=1);
        @time res = solve(solver; scheme_t=1);
        @time CSV.write(save_path * file_name * ".csv", DataFrame(res, :auto));
    end
end

println(save_path * "  finish!");