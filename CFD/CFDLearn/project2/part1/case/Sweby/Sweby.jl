# author: bcynuaa
# date: 2022/11/25

include("..//head.jl");

method = "Limiter";
Limiter = "Sweby";
beta::Float64 = 1.5;

@time result::Vector{Matrix{Float64}} = solve(case, grid; method=method, Limiter=Limiter, Sweby_beta = beta);
@time generate_gif(grid, result);