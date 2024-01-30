"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-13 14:23:53
 # @ license: MIT
 # @ description: shared configuration file
 """

include("..//src//Simulation.jl");

const cavity_len::Float64 = 1.;
const dt::Float64 = 0.005;
const u_upper::Float64 = 1.;
const CFL::Float64 = 0.9;
const total_time::Float64 = 10.;

function makeCaseLight(N::Int64, Re::Float64)::Case
    return makeCase(cavity_len, N, u_upper, Re, CFL, dt);
end

function makeDataPath(case::Case)::String
    data_path::String = "..//data//N$(case.N_)Re$(Int64(floor(case.Re_)))//";
    if !isdir(data_path)
        mkpath(data_path);
    end
    return data_path;
end