# author: bcynuaa
# date: 2022/11/10

module Flux # module begin

include("Grid.jl");
include("Param.jl");
using .Grid;
using .Param;

###################################################################################################

"""
# laplace_W(grid::Grid2D, W::Matrix{Float64})::Matrix{Float64}

**to attain the laplace W - ΔW**

- `grid::Grid2D` -> the grid information
- `W::Matrix{Float64}` -> the W conservative value matrix

return

`W_laplace::Matrix{Float64}` -> the ΔW (∇²W)
"""
function laplace_W(grid::Grid2D, W::Matrix{Float64})::Matrix{Float64}
    W_laplace::Matrix{Float64} = zeros(Float64, grid.ncells, 4);
    laplace::Vector{Float64} = zeros(Float64, 4);
    k::Int64 = 0;
    p::Int64 = 0;
    for edge = 1: grid.nedges
        k, p = grid.iedge[edge, 3: 4];
        if p in [-1, -2]
            continue
        else
            laplace = W[p, :] - W[k, :];
            W_laplace[k, :] .+= laplace;
            W_laplace[p, :] .-= laplace;
        end
    end
    return W_laplace;
end

###################################################################################################

"""
# flux_on_edge(grid::Grid2D, edge_ID::Int64, w_edge_expand::Vector{Float64})::Tuple{Float64, Vector{Float64}}

**calculate flux on the edge**

- `grid::Grid2D` -> the grid of the calculation
- `edge_ID::Int64` -> the ID of the edge
- `w_edge_expand::Vector{Float64}` -> expanded conservative value w

return

- `Z::Float64` -> the Z value on the edge
- `flux::Vector{Float64}` -> the flux on the edge
"""
function flux_on_edge(
    grid::Grid2D, edge_ID::Int64,
    w_edge_expand::Vector{Float64}
)::Tuple{Float64, Vector{Float64}}
    dx::Float64 = grid.dx[edge_ID];
    dy::Float64 = grid.dy[edge_ID];
    rho::Float64 = w_edge_expand[1];
    U::Float64 = w_edge_expand[2];
    V::Float64 = w_edge_expand[3];
    P::Float64 = w_edge_expand[5];
    H::Float64 = w_edge_expand[6];
    Z::Float64 = U*dy - V*dx;
    flux::Vector{Float64} = [
        Z*rho,
        Z*rho*U + P*dy,
        Z*rho*V - P*dx,
        Z*rho*H
    ];
    return Z, flux;
end

"""
# function value_Riemann(case::Case, Un1::Float64, C1::Float64, Un2::Float64, C2::Float64 )::Tuple{Float64, Float64}

**calculate the value of Riemann**

- `case::Case` -> the case of calculation
- `Un1::Float64` -> first value of Un
- `C1::Float64` -> first value of sound speed
- `Un2::Float64` -> second value of Un
- `C2::Float64` -> second value of sound speed

return

- `R_p::Float64` -> R+
- `R_m::Float64` -> R-
"""
function value_Riemann(
    case::Case,
    Un1::Float64, C1::Float64,
    Un2::Float64, C2::Float64 
)::Tuple{Float64, Float64}
    R_p::Float64 = Un1 + 2. * C1 / (case.gamma - 1.);
    R_m::Float64 = Un2 - 2. * C2 / (case.gamma - 1.);
    return R_p, R_m;
end

"""
# flux_by_Riemann(case::Case, grid::Grid2D, edge_ID::Int64, R_p::Float64, R_m::Float64, Ut::Float64, s::Float64)::Vector{Float64}

**calculate flux by Riemann method**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `edge_ID::Int64` -> the ID of the edge
- `C_edge::Float64` -> the sound speed at the edge
- `s_edge::Float64` -> the p/ρ^γ at the edge
- `Un::Float64` -> the edge's normal speed U
- `U::Float64` -> the chosen x direction speed U
- `V::Float64` -> the chosen y direction speed V

return

- `Z::Float64` -> the Z value on the edge by Riemann
- `flux::Vector{Float64}` -> the flux on the edge by Riemann
"""
function flux_by_Riemann(
    case::Case, grid::Grid2D,
    edge_ID::Int64,
    Un_edge::Float64,
    C_edge::Float64,
    s_edge::Float64,
    Un::Float64,
    U::Float64,
    V::Float64
)
    rho_edge::Float64 = (C_edge^2 / s_edge / case.gamma) ^ (1. / (case.gamma - 1.));
    P_edge::Float64 = s_edge * rho_edge^case.gamma;
    U_edge::Float64 = U + (Un_edge - Un) * grid.vec_t[edge_ID, 2];
    V_edge::Float64 = V + (Un_edge - Un) * grid.vec_t[edge_ID, 1];
    Z::Float64 = Un_edge * grid.ds[edge_ID];
    flux::Vector{Float64} = [
        Z * rho_edge,
        Z * rho_edge * U_edge + P_edge * grid.dy[edge_ID],
        Z * rho_edge * V_edge - P_edge * grid.dx[edge_ID],
        Z * (case.gamma / (case.gamma - 1.) * P_edge + rho_edge / 2. * (U_edge^2 + V_edge^2))
    ];
    return Z, flux;
end

###################################################################################################

"""
# flux_at_wall_boundary!(grid::Grid2D, edge_ID::Int64, W_cur_expand::Matrix{Float64}, Q::Matrix{Float64})::Float64

**calculate flux at wall boundary**

- `grid::Grid2D` -> the grid of the calculation
- `edge_ID::Int64` -> the ID of the edge
- `W_cur_expand::Matrix{Float64}` -> the expand conservative value w of all cells
- `Q::Matrix{Float64}` -> the matrix used to store flux

return:

- `C * ds::Float64` -> α for Δtₖ
"""
function flux_at_wall_boundary!(
    grid::Grid2D,
    edge_ID::Int64,
    W_cur_expand::Matrix{Float64},
    Q::Matrix{Float64}
)::Float64
    k::Int64 = grid.iedge[edge_ID, 3];
    P::Float64, C::Float64 = W_cur_expand[k, [5, 7]];
    dx::Float64 = grid.dx[edge_ID];
    dy::Float64 = grid.dy[edge_ID];
    ds::Float64 = grid.ds[edge_ID];
    flux::Vector{Float64} = [0., P*dy, -P*dx, 0.];
    Q[k, :] .+= flux;
    return C * ds;
end

"""
# function flux_at_farfield_boundary!(case::Case, grid::Grid2D, edge_ID::Int64, W_cur_expand::Matrix{Float64}, Q::Matrix{Float64})::Int64

**calculate flux at far field**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `edge_ID::Int64` -> the ID of the edge
- `W_cur_expand::Matrix{Float64}` -> the expand conservative value w of all cells
- `Q::Matrix{Float64}` -> the matrix used to store flux

return 

`abs(Z) + C*grid.ds[edge_ID]::Float64` -> α on the edge;
"""
function flux_at_farfield_boundary!(
    case::Case, grid::Grid2D,
    edge_ID::Int64,
    W_cur_expand::Matrix{Float64},
    Q::Matrix{Float64}
)::Float64
    k::Int64 = grid.iedge[edge_ID, 3];
    w_k_expand::Vector{Float64} = W_cur_expand[k, :];
    w_inf_expand::Vector{Float64} = case.w_inf_expand;

    C_k::Float64 = w_k_expand[7];
    rho_k::Float64 = w_k_expand[1];
    P_k::Float64 = w_inf_expand[5];

    C_inf::Float64 = w_inf_expand[7];
    rho_inf::Float64 = w_inf_expand[1];
    P_inf::Float64 = w_inf_expand[5];

    Un_k::Float64 = grid.vec_n[edge_ID, :]' * w_k_expand[2: 3];
    Un_inf::Float64 = grid.vec_n[edge_ID, :]' * w_inf_expand[2: 3];
    
    R_p::Float64, R_m::Float64 = value_Riemann(case, Un_k, C_k, Un_inf, C_inf);
    Un_edge::Float64 = (R_p + R_m) / 2.;
    C_edge::Float64 = (case.gamma - 1.) / 4. * (R_p - R_m);
    Man_inf::Float64 = abs(Un_inf / C_inf);
    
    Z::Float64 = 0.;
    C::Float64 = 0.;
    flux::Vector{Float64} = zeros(Float64, 4);
    if Un_edge <= 0.
        if Man_inf <= 1. # subsonic inflow
            Z, flux = flux_by_Riemann(
                case, grid,
                edge_ID,
                Un_edge,
                C_edge,
                P_inf / rho_inf^case.gamma,
                Un_inf,
                case.w_inf_expand[2],
                case.w_inf_expand[3]
            );
            C = C_edge;
        else
            Z, flux = flux_on_edge(grid, edge_ID, w_inf_expand);
            C = C_inf;
        end
    else
        if Man_inf <= 1. # subsonic outflow
            Z, flux = flux_by_Riemann(
                case, grid,
                edge_ID,
                Un_edge,
                C_edge,
                P_k / rho_k^case.gamma,
                Un_k,
                w_k_expand[2],
                w_k_expand[3]
            );
            C = C_edge;
        else # supersonic outflow
            Z, flux = flux_on_edge(grid, edge_ID, w_k_expand);
            C = C_k;
        end
    end

    Q[k, :] .+= flux;
    return abs(Z) + C*grid.ds[edge_ID];
end

"""
# flux_between_cells_use_in_rk4!(case::Case, grid::Grid2D, edge_ID::Int64, W_cur::Matrix{Float64}, Q::Matrix{Float64})::Tuple{Float64, Vector{Float64}}

**calculate flux between 2 cells only**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `edge_ID::Int64` -> the ID of the edge
- `W_cur::Matrix{Float64}` -> the conservative value w of all cells
- `Q::Matrix{Float64}` -> the matrix used to store flux

return:

`Z::Float64` -> Z value of the edge 
- `w_edge_expand::Vector{Float64}` -> the expand conservative value on the edge
"""
function flux_between_cells_use_in_rk4!(
    case::Case, grid::Grid2D,
    edge_ID::Int64,
    W_cur::Matrix{Float64},
    Q::Matrix{Float64}
)::Tuple{Float64, Vector{Float64}}
    k::Int64 = grid.iedge[edge_ID, 3];
    p::Int64 = grid.iedge[edge_ID, 4];
    w_k::Vector{Float64} = W_cur[k, :];
    w_p::Vector{Float64} = W_cur[p, :];
    w_edge::Vector{Float64} = (w_k .+ w_p) ./ 2.;
    w_edge_expand::Vector{Float64} = expand_w(case, w_edge);

    Z::Float64, flux::Vector{Float64} = flux_on_edge(grid, edge_ID, w_edge_expand);
    Q[k, :] .+= flux;
    Q[p, :] .-= flux;
    return Z, w_edge_expand;
end

"""
# function flux_between_cells!(case::Case, grid::Grid2D, edge_ID::Int64, W_cur::Matrix{Float64}, W_cur_expand::Matrix{Float64}, W_cur_laplace::Matrix{Float64}, Q::Matrix{Float64}, D2::Matrix{Float64}, D4::Matrix{Float64}, t::Vector{Float64})::Int64

**calculate the flux between 2 cells together with artifical viscosity and tk**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `edge_ID::Int64` -> the ID of the edge
- `W_cur::Matrix{Float64}` -> the conservative value w of all cells
- `W_cur_expand::Matrix{Float64}` -> the expand conservative value w of all cells
- `W_cur_laplace::Matrix{Float64}` -> the laplace conservative value w of all cells
- `Q::Matrix{Float64}` -> the matrix used to store flux
- `D2::Matrix{Float64}` -> 2 order artifical viscosity
- `D4::Matrix{Float64}` -> 4 order artifical viscosity
- `t::Vector{Float64}` -> dt vector

return `0`;
"""
function flux_between_cells!(
    case::Case, grid::Grid2D,
    edge_ID::Int64,
    W_cur::Matrix{Float64},
    W_cur_expand::Matrix{Float64},
    W_cur_laplace::Matrix{Float64},
    Q::Matrix{Float64},
    D2::Matrix{Float64},
    D4::Matrix{Float64},
    t::Vector{Float64}
)::Int64
    Z::Float64, w_edge_expand::Vector{Float64} = flux_between_cells_use_in_rk4!(
        case, grid,
        edge_ID, 
        W_cur,
        Q
    );
    k::Int64 = grid.iedge[edge_ID, 3];
    p::Int64 = grid.iedge[edge_ID, 4];
    P_k::Float64 = W_cur_expand[k, 5];
    P_p::Float64 = W_cur_expand[p, 5];
    C::Float64 = w_edge_expand[7];
    ds::Float64 = grid.ds[edge_ID];
    
    alpha_edge::Float64 = abs(Z) + C*ds;
    nu_edge::Float64 = abs( (P_p - P_k) / (P_p + P_k) );
    epsilon2_edge::Float64 = case.k2 * nu_edge;
    epsilon4_edge::Float64 = max(0., case.k4 - epsilon2_edge);
    
    d2::Vector{Float64} = alpha_edge * epsilon2_edge * (W_cur[p, :] .- W_cur[k, :]);
    d4::Vector{Float64} = -alpha_edge * epsilon4_edge * (W_cur_laplace[p, :] .- W_cur_laplace[k, :]);
    item::Int64 = 0;
    value_sin::Float64 = 1.;
    
    if k in grid.wall_cell_ID
        item = k;
        value_sin = abs(grid.vec_t[edge_ID, :]' * grid.vec_n_wall_cell[item, :]);
    elseif p in grid.wall_cell_ID
        item = p;
        value_sin = abs(grid.vec_t[edge_ID, :]' * grid.vec_n_wall_cell[item, :]);
    end
    d2 *= value_sin;
    d4 *= value_sin;
    D2[k, :] .+= d2;
    D2[p, :] .-= d2;
    D4[k, :] .+= d4;
    D4[p, :] .-= d4;

    t[k] += alpha_edge;
    t[p] += alpha_edge;
    return 0;
end

###################################################################################################

"""
# each_step_use_in_rk4(case::Case, grid::Grid2D, W_cur::Matrix{Float64})::Matrix{Float64}

**the each step use in Runge Kutta 4**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `W_cur::Matrix{Float64}` -> the current step of W

return:

`Q::Matrix{Float64}` -> the flux of all cells
"""
function each_step_use_in_rk4(
    case::Case, grid::Grid2D, W_cur::Matrix{Float64}
)::Matrix{Float64}
    W_cur_expand::Matrix{Float64} = expand_W(case, W_cur);
    Q::Matrix{Float64} = zeros(Float64, grid.ncells, 4);
    for edge_ID = 1: grid.nedges
        if grid.iedge[edge_ID, 4] == -1 # wall boundary
            flux_at_wall_boundary!(
                grid, edge_ID, W_cur_expand, Q
            );
        elseif grid.iedge[edge_ID, 4] == -2 # far field boundary
            flux_at_farfield_boundary!(
                case, grid,
                edge_ID,
                W_cur_expand,
                Q
            );
        else
            flux_between_cells_use_in_rk4!(
                case, grid,
                edge_ID,
                W_cur,
                Q
            );
        end
    end
    return Q;
end

"""
# each_step(case::Case, grid::Grid2D, W_cur::Matrix{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}

**each step include flux and artifical viscosity**

- `case::Case` -> the Case contains parameters
- `grid::Grid2D` -> the grid of the calculation
- `W_cur::Matrix{Float64}` -> the current step of W

return:

`Q::Matrix{Float64}` -> the flux of all cells
`D::Matrix{Float64}` -> the artifical viscosity of all cells
`t::Vector{Float64}` -> the local Δt of all cells
"""
function each_step(
    case::Case, grid::Grid2D, W_cur::Matrix{Float64}
)::Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}
    W_cur_expand::Matrix{Float64} = expand_W(case, W_cur);
    W_cur_laplace::Matrix{Float64} = laplace_W(grid, W_cur);
    Q::Matrix{Float64} = zeros(Float64, grid.ncells, 4);
    D2::Matrix{Float64} = zeros(Float64, grid.ncells, 4);
    D4::Matrix{Float64} = zeros(Float64, grid.ncells, 4);
    t::Vector{Float64} = zeros(Float64, grid.ncells);
    for edge_ID = 1: grid.nedges
        if grid.iedge[edge_ID, 4] == -1 # wall boundary
            t[grid.iedge[edge_ID, 3]] += flux_at_wall_boundary!(
                grid, edge_ID, W_cur_expand, Q
            );
        elseif grid.iedge[edge_ID, 4] == -2 # far field boundary
            t[grid.iedge[edge_ID, 3]] += flux_at_farfield_boundary!(
                case, grid,
                edge_ID,
                W_cur_expand,
                Q
            );
        else
            flux_between_cells!(
                case, grid,
                edge_ID,
                W_cur,
                W_cur_expand,
                W_cur_laplace,
                Q,
                D2, D4,
                t
            );
        end
    end
    t = case.CFL * grid.vol ./ t;
    D::Matrix{Float64} = D2 .+ D4;
    return Q, D, t;
end

###################################################################################################

export each_step, each_step_use_in_rk4;
export Case, Grid2D, expand_w, expand_W;

end # module end