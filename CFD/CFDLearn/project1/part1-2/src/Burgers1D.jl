# auther: bcynuaa
# date: 2022/10/19

module Burgers1D # module begin

include("Base1D.jl");
using .Base1D;

struct Solver
    grid::Grid
    nu::Float64
    u0::Vector{Float64}
    func_bc::Function
    bc_type::Int64
end

function Solver(
    grid::Grid, nu::Float64,
    func::Function,
    func_bc::Function; bc_type::Int64=1
)::Solver
    return Solver(grid, nu, func.(grid.x), func_bc, bc_type);
end

function Solver(
    x::Vector{Float64}, t::Vector{Float64},
    nu::Float64, u0::Vector{Float64},
    func_bc::Function; bc_type::Int64=1
)::Solever
    grid::Grid = Grid(x, t);
    return Solver(grid, nu, u0, func_bc, bc_type);
end

function Solver(
    x::Vector{Float64}, t::Vector{Float64},
    nu::Float64, func::Function,
    func_bc::Function; bc_type::Int64=1
)::Solever
    grid::Grid = Grid(x, t);
    return Solver(grid, nu, func.(grid.x), func_bc, bc_type);
end

###################################################################################################

function solve_implicit!(
    s::Solver, u::Matrix{Float64},
    A1::SparseMatrixCSC{Float64, Int64},
    A2::SparseMatrixCSC{Float64, Int64}
)::Int64
    dt::Float64 = 0.;
    if s.bc_type == 1
        for k = 1:s.grid.Nt-1
            dt = s.grid.t[k+1] - s.grid.t[k];
            u[k+1, 1] = s.func_bc(s.grid.t[k+1]);
            B = I + dt*(spdiagm(u[k, :])*A1 - s.nu*A2);
            u[k+1, 2:end] = B[2:end, 2:end] \ Vector(u[k, 2:end] - u[k+1, 1]*B[2:end, 1]);
        end
    end
    return 0;
end

function solve_explicit!(
    s::Solver, u::Matrix{Float64},
    A1::SparseMatrixCSC{Float64, Int64},
    A2::SparseMatrixCSC{Float64, Int64}
)::Int64
    dt::Float64 = 0.;
    if s.bc_type == 1
        for k = 1:s.grid.Nt-1
            dt = s.grid.t[k+1] - s.grid.t[k];
            u[k+1, :] = (I - dt * spdiagm(u[k, :]) * A1 + s.nu*dt*A2) * u[k, :];
            u[k+1, 1] = s.func_bc(s.grid.t[k+1]);
        end
    end
    return 0;
end

function solve_explicit_conservation!(
    s::Solver, u::Matrix{Float64},
    A1::SparseMatrixCSC{Float64, Int64},
    A2::SparseMatrixCSC{Float64, Int64}
)::Int64
    dt::Float64 = 0.;
    if s.bc_type == 1
        for k = 1:s.grid.Nt-1
            dt = s.grid.t[k+1] - s.grid.t[k];
            u[k+1, :] = (I + s.nu*dt*A2) * u[k, :] .- (dt/2.) * A1 * (u[k, :].^2);
            u[k+1, 1] = s.func_bc(s.grid.t[k+1]);
        end
    end
    return 0;
end

###################################################################################################

function solve(s::Solver; scheme_t::Int64=0, scheme_x::Int64=0, conservation::Int64=0)::Matrix{Float64}
    u::Matrix{Float64} = zeros(s.grid.Nt, s.grid.Nx);
    u[1, :] = s.u0;
    A1::SparseMatrixCSC{Float64, Int64} = spzeros(1);
    if scheme_x==0
        A1 = central_SparseMatrix(s.grid);
    elseif scheme_x==-1
        A1 = backward_SparseMatrix(s.grid);
    else
        A1 = forward_SparseMatrix(s.grid)l
    end
    A2::SparseMatrixCSC{Float64, Int64} = central_SparseMatrix_order2(s.grid);
    if scheme_t==0
        solve_implicit!(s, u, A1, A2);
    else
        if conservation==0
            solve_explicit!(s, u, A1, A2);
        else
            solve_explicit_conservation!(s, u, A1, A2);
        end
    end
    return u;
end

###################################################################################################

export Grid, Solver, solve;

end # module end