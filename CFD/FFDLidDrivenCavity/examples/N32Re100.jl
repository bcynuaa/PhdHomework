"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-13 14:46:15
 # @ license: MIT
 # @ description: N = 32, Re = 100
 """

include("SharedConfig.jl");

N::Int64 = 32;
Re::Float64 = 100.;

case::Case = makeCaseLight(N, Re);
data_path::String = makeDataPath(case);

@time simulation(case, total_time, data_path);