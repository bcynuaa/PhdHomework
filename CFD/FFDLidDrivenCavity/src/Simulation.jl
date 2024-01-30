"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-12 20:46:16
 # @ license: MIT
 # @ description: contains the main function
 """

using LinearAlgebra;
using CSV;
using DataFrames;

include("Case.jl");
include("SparseMatrixInIteration.jl");
include("Advection.jl");
include("Diffusion.jl");
include("Projection.jl");

function iteration!(
    case::Case, step::Int64, 
    u::Matrix{Float64}, v::Matrix{Float64},
    matrixA::SparseMatrixCSC{Float64, Int64},
    matrixA_projection::SparseMatrixCSC{Float64, Int64},
    data_path::String,
)::Int64
    advection!(case, u, v);
    diffusion!(case, u, v, matrixA);
    projection!(case, u, v, matrixA_projection);
    filename1::String = data_path * "u_" * string(step) * ".csv";
    filename2::String = data_path * "v_" * string(step) * ".csv";
    CSV.write(filename1, DataFrame(u, :auto));
    CSV.write(filename2, DataFrame(v, :auto));
    print(step, " ");
    return 0;
end

function simulation(case::Case, total_time::Float64, data_path::String)::Int64
    total_step::Int64 = floor(total_time / case.dt_);
    u::Matrix{Float64}, v::Matrix{Float64} = initializeVelcoity(case);
    CSV.write(data_path * "u_0.csv", DataFrame(u, :auto));
    CSV.write(data_path * "v_0.csv", DataFrame(v, :auto));

    println("start simulation");
    println("total step: $total_step");
    println("make matrix A");
    @time matrixA::SparseMatrixCSC{Float64, Int64} = makeMatrixA(case);
    @time matrixA_projection::SparseMatrixCSC{Float64, Int64} = makeMatrixAProjection(case);
    println("make matrix A and A projection done, start simulating");
    for step = 1: total_step
        iteration!(case, step, u, v, matrixA, matrixA_projection, data_path);
    end
    println("\nsimulation done");
    return 0;
end