"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-13 14:57:46
 # @ license: MIT
 # @ description: the reader
 """

using DelimitedFiles;
import PyPlot as plt;

rcParams = plt.PyDict(plt.matplotlib."rcParams");
rcParams["text.usetex"] = true;
const figsize::Tuple{Float64, Float64} = (8, 6);
const labelsize::Float64 = 25;
const titlesize::Float64 = 30;
const levels::Int64 = 50;
const cmap::String = "jet";

const N_list::Vector{Int64} = [16, 32, 64];
const flag_list::Vector{String} = ["u", "v"];
const last_step::Int64 = 2000;
const dt::Float64 = 0.005;
const skip_step::Int64 = 20;
const cavity_len::Float64 = 1.0;

function getAtStep(N::Int64, Re::Float64, flag::String, step::Int64)
    data_path::String = "..//data//N$(N)Re$(Int64(floor(Re)))//"
    data::Matrix{Float64} = readdlm(data_path * flag * "_$(step).csv", ',', Float64, '\n'; skipstart=1);
    return data;
end

function getXY(N::Int64)::Tuple{Vector{Float64}, Vector{Float64}}
    x::Vector{Float64} = LinRange(0, cavity_len, N);
    y::Vector{Float64} = LinRange(0, cavity_len, N);
    return x, y
end

function getXorY(N::Int64)::Vector{Float64}
    return LinRange(0, cavity_len, N);
end

function getUMiddle(N::Int64, Re::Float64)::Vector{Float64}
    u::Matrix{Float64} = getAtStep(N, Re, "u", last_step);
    return u[:, Int64(N/2)];
end

function getVMiddle(N::Int64, Re::Float64)::Vector{Float64}
    v::Matrix{Float64} = getAtStep(N, Re, "v", last_step);
    return v[Int64(N/2), :];
end

const ref_data_U::Matrix{Float64} = readdlm("..//data//reference//Umiddle.dat", Float64; skipstart=1);
const ref_y::Vector{Float64} = ref_data_U[:, 1];
const ref_u_Re100::Vector{Float64} = ref_data_U[:, 2];
const ref_u_Re400::Vector{Float64} = ref_data_U[:, 3];

const ref_data_V::Matrix{Float64} = readdlm("..//data//reference//Vmiddle.dat", Float64; skipstart=1);
const ref_x::Vector{Float64} = ref_data_V[:, 1];
const ref_v_Re100::Vector{Float64} = ref_data_V[:, 2];
const ref_v_Re400::Vector{Float64} = ref_data_V[:, 3];