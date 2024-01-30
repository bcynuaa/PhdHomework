# author: bcynuaa
# date: 2022/10/04

module RK4 # begin module

include("Flux.jl");
using .Flux;

###################################################################################################

# the coefficient of runge kutta 4
const alpha_rk4::Vector{Float64} = [1/4., 1/3., 1/2., 1.];

"""
# function runge_kutta4(case::Case, grid::Grid2D, W_cur::Matrix{Float64})::Tuple{Matrix{Float64}, Vector{Float64}, Float64}

**the Runge Kutta 4 method**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `W_cur::Matrix{Float64}` -> the current step of W

return

- `W::Matrix{Float64}` -> next W
- `t::Matrix{Float64}` -> the time t
- `dt::Vector{Float64}` -> the dt of this time
"""
function runge_kutta4(case::Case, grid::Grid2D, W_cur::Matrix{Float64})::Tuple{Matrix{Float64}, Vector{Float64}, Float64}
    Q::Matrix{Float64}, D::Matrix{Float64}, t::Vector{Float64} = each_step(case, grid, W_cur);
    dt::Float64 = minimum(t);
    W::Matrix{Float64} = W_cur .- (Q .- D) ./ grid.vol .* alpha_rk4[1] * dt;
    for m = 2: 4
        Q = each_step_use_in_rk4(case, grid, W);
        W = W_cur .- (Q .- D) ./ grid.vol .* alpha_rk4[m] * dt;
    end
    return W, t, dt;
end

###################################################################################################

export runge_kutta4;
export Case, Grid2D, expand_w, expand_W;

end # end module