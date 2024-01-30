# auther: bcynuaa
# date: 2022/10/12

using CSV;
using DataFrames;

c::Vector{Float64} = [0.1, 0.5, 1., 5., 10.];
c_name::Vector{String} = "c=" .* string.(c);
func_name::Vector{String} = ["Triangle", "Step"];
save_path::String = "..//..//data//SingleWave//implicit_backward//";

dx::Float64 = 0.1;
dt::Float64 = 0.005;
t_end::Float64 = 5.;
x::Vector{Float64} = Vector(-10.: dx: 10.);
t::Vector{Float64} = Vector(0.: dt: t_end);

###################################################################################################

using Plots;
gr();
image_path::String = "..//..//image//SingleWave//implicit_backward//"

if isdir(image_path) != true
    mkpath(image_path)
end

skip::Int64 = 10;

for j = 1:length(c)
    for i = 1:length(func_name)
        println("c=$(c_name[j]), $(func_name[i]) IC");
        file_name = func_name[i] * "_" * c_name[j];
        @time data = Matrix( CSV.read(save_path * file_name * ".csv", DataFrame) )
        @time anim = @animate for k = 1:skip:size(data)[1]
            tk = string(round(t[k]; digits=2));
            plot(
                x, data[k, :],
                color=:blue, lw=2,
                title="$(func_name[i]) IC \$ t= $tk s\$",
                label="\$ \\frac{\\partial u}{\\partial t} + $(c[j])\\frac{\\partial u}{\\partial x}=0 \$"
            );
            xlabel!("\$ x \$")
            ylabel!("\$ u \$")
        end
        @time gif(anim, image_path * file_name * ".gif", fps=50)
    end
end

println(image_path * " finish");