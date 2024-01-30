# author: bcynuaa
# date: 2022/11/25

include("..//head.jl");

method = "FVS";

@time result::Vector{Matrix{Float64}} = solve(case, grid; method = method);
@time generate_gif(grid, result);