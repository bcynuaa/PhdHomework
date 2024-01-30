# auther: bcynuaa
# date: 2022/10/18

###################################################################################################

module Heat1D # module begin

include("Base1D.jl");
using .Base1D;

###################################################################################################

struct Solver
    grid::Grid
    c::Float64
    u0::Vector{Float64}
    func_bc::Function
    bc_type::Int64
end

function Solver(
    grid::Grid, c::Float64,
    func::Function,
    func_bc::Function; bc_type::Int64=1
)::Solver
    return Solver(grid, c, func.(grid.x), func_bc, bc_type);
end

function Solver(
    x::Vector{Float64}, t::Vector{Float64},
    c::Float64, u0::Vector{Float64},
    func_bc::Function; bc_type::Int64=1
)::Solever
    grid::Grid = Grid(x, t);
    return Solver(grid, c, u0, func_bc, bc_type);
end

function Solver(
    x::Vector{Float64}, t::Vector{Float64},
    c::Float64, func::Function,
    func_bc::Function; bc_type::Int64=1
)::Solever
    grid::Grid = Grid(x, t);
    return Solver(grid, c, func.(grid.x), func_bc, bc_type);
end

###################################################################################################

function modify_A_bc2!(s::Solver, A::SparseMatrixCSC{Float64, Int64})
    modify::Matrix{Float64} = zeros(Float64, 3, 4);
    modify[2, 1:3] = interpolation_coefficient3(s.grid.x, 2; flag=0);
    modify[3, 2:4] = interpolation_coefficient3(s.grid.x, 3; flag=0);
    A[1, 1:4] = transpose( transpose( interpolation_coefficient3(s.grid.x, 2; flag=-1) ) * modify );
end

###################################################################################################

function solve_implicit!(
    s::Solver, u::Matrix{Float64},
    A::SparseMatrixCSC{Float64, Int64}
)::Int64
    dt::Float64 = 0.;
    if s.bc_type == 1
        for k = 1:s.grid.Nt-1
            dt = s.grid.t[k+1] - s.grid.t[k];
            u[k+1, 1] = s.func_bc(s.grid.t[k+1]);
            u[k+1, 2:end] = (I - s.c*dt*A[2:end, 2:end]) \ Vector(u[k, 2:end] - (-s.c*dt*u[k+1, 1]) * A[2:end, 1]);
        end
    elseif s.bc_type == 2
        modify_A_bc2!(s, A);
        for k = 1:s.grid.Nt-1
            dt = s.grid.t[k+1] - s.grid.t[k];
            tmp = u[k, :];
            tmp[1] += s.c*dt*interpolation_coefficient3(s.grid.x, 2; flag=-1)[1] * s.func_bc(s.grid.t[k+1]);
            u[k+1, :] = (I - s.c*dt*A) \ tmp;
        end
    end
    return 0;
end

function solve_explicit!(
    s::Solver, u::Matrix{Float64},
    A::SparseMatrixCSC{Float64, Int64}
)::Int64
    dt::Float64 = 0.;
    if s.bc_type == 1
        for k = 1:s.grid.Nt-1
            dt = s.grid.t[k+1] - s.grid.t[k];
            u[k+1, 1] = s.func_bc(s.grid.t[k+1]);
            u[k+1, 2:end] = (dt*s.c*u[k, 1]) * A[2:end, 1] + (I + s.c*dt*A[2:end, 2:end]) * u[k, 2:end];
        end
    elseif s.bc_type == 2
        modify_A_bc2!(s, A);
        for k = 1:s.grid.Nt-1
            dt = s.grid.t[k+1] - s.grid.t[k];
            u[k+1, :] = (I + s.c*dt*A) * u[k, :];
            u[k+1, 1] += s.c*dt*interpolation_coefficient3(s.grid.x, 2; flag=-1)[1] * s.func_bc(s.grid.t[k]);
        end
    end
    return 0;
end

###################################################################################################

function solve(s::Solver; scheme_t::Int64=0)::Matrix{Float64}
    u::Matrix{Float64} = zeros(s.grid.Nt, s.grid.Nx);
    u[1, :] = s.u0;
    A::SparseMatrixCSC{Float64, Int64} = central_SparseMatrix_order2(s.grid);
    if scheme_t==0
        solve_implicit!(s, u, A);
    else
        solve_explicit!(s, u, A);
    end
    return u;
end

###################################################################################################

export Grid, Solver, solve;

end # module end