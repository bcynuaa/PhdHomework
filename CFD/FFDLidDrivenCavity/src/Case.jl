"""
 # @ author: bcynuaa | bcynuaa@163.com
 # @ date: 2023-06-12 19:08:00
 # @ license: MIT
 # @ description: contains the case struct and some functions
 """

struct Case
    cavity_len_::Float64 # cavity length
    N_::Int64 # grid number
    dx_::Float64 # grid size

    u_upper_::Float64 # upper wall velocity

    Re_::Float64 # Reynolds number
    nu_::Float64 # viscosity

    CFL_::Float64 # CFL number
    dt_::Float64 # time step
    lambda_::Float64 # relaxation factor
end

function makeCase(
    cavity_len::Float64 = 1.,
    N::Int64 = 101,
    u_upper::Float64 = 1.,
    Re::Float64 = 100.,
    CFL::Float64 = 0.9,
    dt::Float64 = 0.005,
)::Case
    dx::Float64 = cavity_len / (N-1);
    nu::Float64 = u_upper * cavity_len / Re;
    lambda::Float64 = nu * dt / dx^2;
    return Case(cavity_len, N, dx, u_upper, Re, nu, CFL, dt, lambda);
end

function initializeVelcoity(case::Case)::Tuple{Matrix{Float64}, Matrix{Float64}}
    N::Int64 = case.N_;
    u0::Matrix{Float64} = zeros(Float64, N, N);
    v0::Matrix{Float64} = zeros(Float64, N, N);
    u0[N, :] .= case.u_upper_;
    return u0, v0;
end

