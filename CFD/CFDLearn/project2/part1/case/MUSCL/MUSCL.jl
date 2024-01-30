# author: bcynuaa
# date: 2022/11/25

include("..//head.jl");

method = "Limiter";
Limiter = "MUSCL";

@time result::Vector{Matrix{Float64}} = solve(case, grid; method=method, Limiter=Limiter);
@time generate_gif(grid, result);