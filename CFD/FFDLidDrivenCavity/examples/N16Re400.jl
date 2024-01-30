"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-13 14:48:26
 # @ license: MIT
 # @ description: N = 16, Re = 400
 """

include("SharedConfig.jl");

N::Int64 = 16;
Re::Float64 = 400.;

case::Case = makeCaseLight(N, Re);
data_path::String = makeDataPath(case);

@time simulation(case, total_time, data_path);