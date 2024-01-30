"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-13 14:52:16
 # @ license: MIT
 # @ description: N = 64, Re = 400
 """

 include("SharedConfig.jl");

 N::Int64 = 64;
 Re::Float64 = 400.;
 
 case::Case = makeCaseLight(N, Re);
 data_path::String = makeDataPath(case);
 
 @time simulation(case, total_time, data_path);