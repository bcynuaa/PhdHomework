# author: bcynuaa
# date: 2022/11/25

include("..//src//ShockTube.jl");
using .ShockTube;
using Plots;

x_begin::Float64 = -2.;
x_end::Float64 = 2.;
dx::Float64 = 0.01;
dt::Float64 = 0.0001;
t_end::Float64 = 0.1;

case::Case = Case();
grid::UniformGrid1D = UniformGrid1D(-2., 2., 0.01, 0.1, 0.0001);

method::String = "FVS";
Limiter::String = "Minmod";

const param_name::Vector{String} = ["\$ \\rho \$", "\$ u \$", "\$ p \$"];
const gif_name::Vector{String} = ["rho.gif", "u.gif", "p.gif"];
const png_name::Vector{String} = ["rho.png", "u.png", "p.png"];
const pdf_name::Vector{String} = ["rho.pdf", "u.pdf", "p.pdf"];
const fps::Int64 = 30;

function generate_gif(
    grid::UniformGrid1D,
    result::Vector{Matrix{Float64}}
)::Int64
    for k = 1: 3
        anim = @animate for j = 1: grid.Step+1
            plot(
                grid.x, 
                result[j][:, k],
                color=:blue,
                lw=2,
                label=param_name[k]
            );
            xlabel!("\$x\$");
            ylabel!(param_name[k]);
            title!("\$ t = $(round( grid.t[j]; digits=3 )) s\$");
        end
        gif(anim, gif_name[k], fps=fps);
        plot();
        for i = 1: 200: grid.Step+1
            plot!(
                grid.x,
                result[i][:, k],
                label="\$t=$(round(i*dt; digits=3))s\$",
                lw=3
            );
        end
        xlabel!("\$x\$");
        ylabel!(param_name[k]);
        title!(method * "  " * Limiter);
        savefig(png_name[k]);
        savefig(pdf_name[k]);
    end
    return 0;
end