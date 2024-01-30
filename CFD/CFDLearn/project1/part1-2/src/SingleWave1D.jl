# auther: bcynuaa
# date: 2022/10/12

###################################################################################################

module SingleWave1D # module begin

include("Base1D.jl");
using .Base1D;

###################################################################################################

struct Solver
    grid::Grid
    c::Float64
    u0::Vector{Float64}
end

function Solver(grid::Grid, c::Float64, func::Function)::Solver
    return Solver(grid, c, func.(grid.x));
end

function Solver(x::Vector{Float64}, t::Vector{Float64}, c::Float64, u0::Vector{Float64})::Solver
    grid::Grid = Grid(x, t);
    return Solver(grid, c, u0);
end

function Solver(x::Vector{Float64}, t::Vector{Float64}, c::Float64, func::Function)::Solver
    grid::Grid = Grid(x, t);
    return Solver(grid, c, func);
end

###################################################################################################

function solve_implicit!(s::Solver, u::Matrix{Float64}, A::SparseMatrixCSC{Float64, Int64})::Int64
    dt::Float64 = 0.;
    for k = 1:s.grid.Nt-1
        dt = s.grid.t[k+1] - s.grid.t[k];
        u[k+1, :] = (I + s.c*dt*A) \ u[k, :];
    end
    return 0;
end

function solve_explicit!(s::Solver, u::Matrix{Float64}, A::SparseMatrixCSC{Float64, Int64})::Int64
    dt::Float64 = 0.;
    for k = 1:s.grid.Nt-1
        dt = s.grid.t[k+1] - s.grid.t[k];
        u[k+1, :] = (I - s.c*dt*A) * u[k, :];
    end
    return 0;
end

###################################################################################################

function solve(s::Solver; scheme_t::Int64=0, scheme_x=-1)::Matrix{Float64}
    u::Matrix{Float64} = zeros(s.grid.Nt, s.grid.Nx);
    u[1, :] = s.u0;
    A::SparseMatrixCSC{Float64, Int64} = spzeros(1);
    if scheme_x == -1
        A = backward_SparseMatrix(s.grid);
    elseif scheme_x == 0
        A = central_SparseMatrix(s.grid);
    elseif scheme_x == 1
        A = forward_SparseMatrix(s.grid);
    else
        A = backward_SparseMatrix(s.grid);
    end
    if scheme_t==0
        solve_implicit!(s, u, A);
    else
        solve_explicit!(s, u, A);
    end
    return u;
end

export Grid, Solver, solve;

end # module end