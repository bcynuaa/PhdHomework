# author: bcynuaa
# author: 2022/11/12

function generate_data_name(Ma::Float64, theta::Float64)::String
    return "MA$(Int64(round(Ma*10)))_alpha$(Int64(round(theta*10))).dat";;
end

Ma::Float64 = 0.;
theta::Float64 = 0.;
grid_file::String = "..//..//naca0012.grd";
data_path::String = ".//data//";

STEP::Int64 = 20000;
plot::Bool = true;