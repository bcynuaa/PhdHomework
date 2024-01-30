# auther: bcynuaa
# date: 2022/10/31

###################################################################################################

module BurgersInviscid1D # module begin

using ExportAll;
include("Base1D.jl");
using .Base1D;

struct Solver
    grid::Grid;
    u0::Vector{Float64};
    func_bc::Function;
    bc_type::Int64;
end

function Solver(
    grid::Grid,
    func::Function,
    func_bc::Function; bc_type::Int64=1
)::Solver
    return Solver(grid, func.(grid.x), func_bc, bc_type);
end

function Solver(
    x::Vector{Float64}, t::Vector{Float64},
    u0::Vector{Float64},
    func_bc::Function; bc_type::Int64=1
)::Solever
    grid::Grid = Grid(x, t);
    return Solver(grid, u0, func_bc, bc_type);
end

function Solver(
    x::Vector{Float64}, t::Vector{Float64},
    func::Function,
    func_bc::Function; bc_type::Int64=1
)::Solever
    grid::Grid = Grid(x, t);
    return Solver(grid, func.(grid.x), func_bc, bc_type);
end

###################################################################################################

function residual_Upwind_Nonconservative(x::Vector{Float64}, u::Vector{Float64})::Vector{Float64}
    r::Vector{Float64} = zeros(Float64, length(u));
    r[2: end] = [-(u[j] - u[j-1]) * u[j] / (x[j] - x[j-1]) for j = 2: length(u)];
    r[1] = -(u[2] - u[1]) * u[1] / (x[2] - x[1]);
    return r;
end

function residual_Upwind_Conservative(x::Vector{Float64}, u::Vector{Float64})::Vector{Float64}
    r::Vector{Float64} = zeros(Float64, length(u));
    r[2: end] = [-(u[j]^2 - u[j-1]^2) / (x[j] - x[j-1]) for j = 2: length(u)];
    r[1] = -(u[2]^2 - u[1]^2) / (x[2] - x[1]);
    return r;
end

function residual_Godunov(x::Vector{Float64}, u::Vector{Float64})::Vector{Float64}
    N::Int64 = length(x);
    r::Vector{Float64} = zeros(Float64, N);
    c::Vector{Float64} = (u[2: end] .+ u[1: end-1]) ./ 2.;
    f::Vector{Float64} = zeros(Float64, N-1);
    for j = 1: N-1
        if u[j] >= u[j+1] # shock waves
            if c[j] > 0
                f[j] = u[j]^2 / 2.;
            else
                f[j] = u[j+1]^2 / 2.;
            end
        else # expansion waves
            if u[j] < 0. < u[j+1]
                f[j] = 0;
            elseif c[j] > 0
                u[j]^2 / 2.;
            else
                u[j+1]^2 / 2.;
            end
        end
    end
    r[2: end-1] = -(f[2: end] .- f[1: end-1]) ./ (x[3: end] .- x[1: end-2]) .*2.;
    r[1] = -(u[2]^2 / 2. - u[1]^2 / 2.) / (x[2] - x[1]);
    r[end] = -(u[end]^2 / 2. - u[end-1]^2 / 2.) / (x[end] - x[end-1]);
    return r;
end

function residual_Roe(x::Vector{Float64}, u::Vector{Float64})::Vector{Float64}
    N::Int64 = length(x);
    r::Vector{Float64} = zeros(Float64, N);
    u_::Vector{Float64} = [u[j]==u[j+1] ? u[j] : (u[j] + u[j+1]) / 2. for j = 1: N-1];
    f::Vector{Float64} = [ (u[j+1]^2/2. + u[j]^2/2.) / 2. - abs(u_[j]) * (u[j+1] - u[j]) / 2. for j = 1: N-1];
    r[2: end-1] = -(f[2: end] .- f[1: end-1]) ./ (x[3: end] .- x[1: end-2]) .*2.;
    r[1] = -(u[2]^2 / 2. - u[1]^2 / 2.) / (x[2] - x[1]);
    r[end] = -(u[end]^2 / 2. - u[end-1]^2 / 2.) / (x[end] - x[end-1]);
    return r;
end

function residual_MUSCL(kappa::Float64)::Function
    function residual(x::Vector{Float64}, u::Vector{Float64})::Vector{Float64}
        N::Int64 = length(x);
        r::Vector{Float64} = zeros(Float64, N);
        f::Vector{Float64} = (u.^2) ./ 2.;
        fL::Vector{Float64} = f[2:end-1] .+ (f[2:end-1] .- f[1:end-2]) .* (1-kappa)/4. + (f[3:end] .- f[2:end-1]) .* (1+kappa)/4.;
        fR::Vector{Float64} = f[1:end-2] .- (f[2:end-1] .- f[1:end-2]) .* (1+kappa)/4. - (f[3:end] .- f[2:end-1]) .* (1-kappa)/4.;
        r[2: end-1] = -(fL - fR) ./ (x[3: end] .- x[1: end-2]) .*2.;
        r[1] = -(u[2]^2 / 2. - u[1]^2 / 2.) / (x[2] - x[1]);
        r[end] = -(u[end]^2 / 2. - u[end-1]^2 / 2.) / (x[end] - x[end-1]);
        return r;
    end
    return residual;
end

###################################################################################################

function bc1!(iter::Function, R::Function, s::Solver, u::Matrix{Float64})::Int64
    for k = 1: s.grid.Nt-1
        dt = s.grid.t[k+1] - s.grid.t[k];
        u[k+1, :] = iter(R, s.grid.x, u[k, :], dt);
        u[k+1, 1] = s.func_bc(s.grid.t[k+1]);
    end
    return 0;
end

function bc2!(iter::Function, R::Function, s::Solver, u::Matrix{Float64})::Int64
    for k = 1: s.grid.Nt-1
        dt = s.grid.t[k+1] - s.grid.t[k];
        u[k+1, :] = iter(R, s.grid.x, u[k, :], dt);
        u[k+1, 1] = u[k, 1] - u[k, 1] * dt * s.func_bc(s.grid.t[k+1]);
    end
    return 0;
end

###################################################################################################

MUSCL_KAPPA::Dict = Dict(
    "Upwind_2nd_Order" => -1.,
    "Fromm" => 0.,
    "QUICK"=> 1/2.,
    "Upwind_3nd_Order" => 1/3.,
    "Centered" => 1.
);

function solve(s::Solver; scheme::String="Upwind_Conservative", use_rk4::Bool=false)::Matrix{Float64}
    u::Matrix{Float64} = zeros(s.grid.Nt, s.grid.Nx);
    u[1, :] = s.u0;
    if scheme == "Upwind_Nonconservative"
        R = residual_Upwind_Nonconservative;
    elseif scheme == "Upwind_Conservative"
        R = residual_Upwind_Conservative;
    elseif scheme == "Godunov"
        R = residual_Godunov;
    elseif scheme == "Roe"
        R = residual_Roe;
    else
        R = residual_MUSCL(MUSCL_KAPPA[scheme]);
    end
    if use_rk4 == false
        iter = non_RK4;
    else
        iter = RK4;
    end
    if s.bc_type == 1
        bc1!(iter, R, s, u);
    elseif s.bc_type == 2
        bc2!(iter, R, s, u);
    end
    return u;
end

###################################################################################################

@exportAll();
export Grid, Solver, solve;

end # module end