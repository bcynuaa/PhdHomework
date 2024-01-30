# author: bcynuaa
# date: 2022/10/04

module EulerEq # module begin

include("Post.jl");
using .Post;

###################################################################################################

"""
# function initialize_field(case::Case, grid::Grid2D)::Matrix{Float64}

**generate the initial W of flow field**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation

return

`W::Matrix{Float64}` -> the initial W of flow field
"""
function initialize_field(case::Case, grid::Grid2D)::Matrix{Float64}
    return repeat(case.w_inf', grid.ncells);
end

###################################################################################################

"""
# function solve(case::Case, grid::Grid2D; plot::Bool=false)::Matrix{Float64}

**solve the Euler equation of the case and grid**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
;
- `plot::Bool` -> whether to do the plot, default to false

return 

`W::Matrix{Float64}` -> the final conservative value W
"""
function solve(case::Case, grid::Grid2D; plot::Bool=false)::Matrix{Float64}
    t::Vector{Float64} = zeros(Float64, case.STEP);
    residual::Matrix{Float64} = zeros(Float64, case.STEP, 4);
    initial_residual::Vector{Float64} = zeros(Float64, 4);
    W_cur::Matrix{Float64} = initialize_field(case, grid);
    t_cur::Vector{Float64} = zeros(Float64, grid.ncells);

    @time begin
    for k = 1: case.STEP
        println("STEP: $k    Ma=$(case.Ma_inf)    alpha=$(case.theta)");
        @time begin
            W_next, t_cur, dt = runge_kutta4(case, grid, W_cur);
            residual[k, :] = sqrt.(sum( ((W_next .- W_cur) ./ t_cur).^2; dims=1 ));
            if k == 1
                initial_residual = residual[k, :]
            end
            residual[k, :] ./= initial_residual;
            if maximum(residual[k, :]) < case.error
                println("satisfy error condition < $(case.error) at Step $k");
                break
            end
            t[k] = dt;
            W_cur = copy(W_next);
        println("dt = $(dt)\nresiduals:\n$(residual[k, :])");
        end
        println("----------------------------------------------------------------------------------------------------");
    end
    end
    @time save_data(case, grid, W_cur, residual; plot=plot);
    return W_cur;
end

###################################################################################################

export solve, Case, Grid2D;

end # module end