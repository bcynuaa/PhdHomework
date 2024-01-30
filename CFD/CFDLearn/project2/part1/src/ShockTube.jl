# author: bcynuaa
# date: 2022/11/23

module ShockTube # module begin

include("EulerEq1D.jl");
include("Grid1D.jl");
using .EulerEq1D;
using .Grid1D;

###################################################################################################

struct Case
    gamma::Float64

    rhoL::Float64
    uL::Float64
    pL::Float64

    rhoR::Float64
    uR::Float64
    pR::Float64

    x_discontinuous::Float64

    conL::Vector{Float64}
    conR::Vector{Float64}

    initial_func::Function
end

function Case(
    gamma::Float64 = 1.4,
    rhoL::Float64 = 1e3,
    uL::Float64 = 0.,
    pL::Float64 = 1e5,
    rhoR::Float64 = 0.125e3,
    uR::Float64 = 0.,
    pR::Float64 = 1e4,
    x_discontinuous::Float64 = 0.
)::Case
    f(x::Float64)::Vector{Float64} = x < x_discontinuous ? base2converse(gamma, rhoL, uL, pL) : base2converse(gamma, rhoR, uR, pR);
    function f(x::Vector{Float64})::Matrix{Float64}
        M::Int64 = length(x);
        init::Matrix{Float64} = zeros(Float64, M, num_converse);
        for j = 1: M
            init[j, :] = f(x[j]);
        end
        return init;
    end
    return Case(
        gamma,
        rhoL,
        uL,
        pL,
        rhoR,
        uR,
        pR,
        x_discontinuous,
        base2converse(gamma, rhoL, uL, pL),
        base2converse(gamma, rhoR, uR, pR),
        f
    );
end

###################################################################################################

const alpha_rk4::Vector{Float64} = [1. / 4., 1. / 3., 1. / 2., 1.];

function runge_kutta4(
    case::Case, grid::UniformGrid1D,
    cur::Matrix{Float64},
    Res::Function
)::Matrix{Float64}
    next::Matrix{Float64} = cur .+ Res(case, grid, cur) .* (alpha_rk4[1] * grid.dt );
    for j = 2: 4
        next = cur .+ Res(case, grid, next) .* (alpha_rk4[j] * grid.dt);
    end
    return next;
end

###################################################################################################

function FVS_residuals(
    case::Case,
    grid::UniformGrid1D,
    cur::Matrix{Float64}
)::Matrix{Float64}
    residuals::Matrix{Float64} = zeros(Float64, size(cur));
    F_plus::Matrix{Float64} = F_plus_converse(case.gamma, cur);
    F_minus::Matrix{Float64} = F_minus_converse(case.gamma, cur);
    residuals[2: end-1, :] = -((F_plus[2: end-1, :] .- F_plus[1: end-2, :]) .+ (F_minus[3: end, :] - F_minus[2: end-1, :]) ) / grid.dx;
    return residuals;
end

###################################################################################################

# limiter

const MUSCL_KAPPA::Dict{String, Float64}= Dict(
    "Upwind_2nd_Order" => -1.,
    "Fromm" => 0.,
    "QUICK"=> 1/2.,
    "Upwind_3nd_Order" => 1/3.,
    "Centered" => 1.
);

const small_amount::Float64 = 1e-7;

function psi_MUSCL_KAPPA(r::Float64, scheme::String)::Float64
    kappa::Float64 = MUSCL_KAPPA[scheme];
    return (1. + kappa) / 2. * r + (1. - kappa) / 2.;
end

function generate_psi_MUSCL_KAPPA(scheme::String="Centered")::Function
    function psi(r::Float64)::Float64
        return psi_MUSCL_KAPPA(r, scheme);
    end
    return psi;
end

function psi_Minmod(r::Float64)::Float64
    return max(0., min(1., r));
end

function psi_Superbee(r::Float64)::Float64
    return max(0., min(2. * r, 1.), min(r, 2.));
end

function psi_VanLeer(r::Float64)::Float64
    return (r + abs(r)) / (1. + r);
end

function psi_MUSCL(r::Float64)::Float64
    return max(0., min(2. * r, (r + 1.) / 2., 2.));
end

function psi_Sweby(r::Float64, beta::Float64 = 1.5)
    return max(0., min(beta * r, 1.), min(r, beta));
end

function generate_psi_Sweby(beta::Float64 = 1.5)
    function psi(r::Float64)::Float64
        return psi_Sweby(r, beta);
    end
    return psi;
end

function Limiter_residuals(
    case::Case,
    grid::UniformGrid1D,
    cur::Matrix{Float64},
    psi::Function
)::Matrix{Float64}
    residuals::Matrix{Float64} = zeros(Float64, size(cur));
    cur_diff::Matrix{Float64} = cur[2: end, :] .- cur[1: end-1, :];
    
    rL::Matrix{Float64} = (cur_diff[2: end, :] .* cur_diff[1: end-1, :]) ./ (cur_diff[1: end-1, :].^2 .+ small_amount);
    rR::Matrix{Float64} = (cur_diff[1:end-1, :] .* cur_diff[2:end, :]) ./ (cur_diff[2: end, :].^2 .+ small_amount);
    
    fL::Matrix{Float64} = cur[2: end-1, :] .+ 1. / 2. * psi.(rL) .* cur_diff[1: end-1, :];
    fR::Matrix{Float64} = cur[2: end-1, :] .- 1. / 2. * psi.(rR) .* cur_diff[2: end, :];
    
    FL::Matrix{Float64} = F_plus_converse(case.gamma, fL);
    FR::Matrix{Float64} = F_minus_converse(case.gamma, fR);
    
    residuals[3: end-2, :] = -( FL[2: end-1, :] .+ FR[3: end, :] .- FL[1: end-2, :] .- FR[2: end-1, :] ) ./ grid.dx;
    # residuals[2, :] = -(FL[1, :] .+ FR[2, :] .- F_converse(case.gamma, fR[1, :])) ./ grid.dx;
    # residuals[end-1, :] = -(F_converse(case.gamma, fL[end, :]) .- FL[end-1, :] .- FR[end, :]) ./ grid.dx;
    residuals[2, :] = -(F_plus_converse(case.gamma, cur[2, :]) .- F_plus_converse(case.gamma, cur[1, :]) .+ F_minus_converse(case.gamma, cur[3, :]) .- F_minus_converse(case.gamma, cur[2, :])) ./ grid.dx;
    residuals[end-1, :] = -(F_plus_converse(case.gamma, cur[end-1, :]) .- F_plus_converse(case.gamma, cur[end-2, :]) .+ F_minus_converse(case.gamma, cur[end, :]) .- F_minus_converse(case.gamma, cur[end-1, :])) ./ grid.dx;
    return residuals;
end

function generate_Limiter_residuals(
    Limiter::String="Minmod";
    MUSCL_scheme::String = "Centered",
    Sweby_beta::Float64 = 1.5
)::Function
    if Limiter == "MUSCL_KAPPA"
        psi = generate_psi_MUSCL_KAPPA(MUSCL_scheme);
    elseif Limiter == "Minmod"
        psi = psi_Minmod;
    elseif Limiter == "Superbee"
        psi = psi_Superbee
    elseif Limiter == "VanLeer"
        psi = psi_VanLeer;
    elseif Limiter == "MUSCL"
        psi = psi_MUSCL;
    elseif Limiter == "Sweby"
        psi = generate_psi_Sweby(Sweby_beta);
    end
    function limiter_res(
        case::Case,
        grid::UniformGrid1D,
        cur::Matrix{Float64}
    )::Matrix{Float64}
        return Limiter_residuals(case, grid, cur, psi);
    end
    return limiter_res;
end

###################################################################################################

```
method:

- FVS
- Limiter

Limiter:

- MUSCL_KAPPA: should specify the MUSCL_scheme, default to "Centered" √
- Minmod √
- Superbee √
- VanLeer √
- MUSCL √
- Sweby: should specify the Sweby_beta, default to 1.5 √
```
function solve(
    case::Case, 
    grid::UniformGrid1D;
    method::String = "FVS",
    Limiter::String = "Minmod",
    MUSCL_scheme::String = "Centered",
    Sweby_beta::Float64 = 1.5
)::Vector{Matrix{Float64}}
    result::Vector{Matrix{Float64}} = fill(zeros(Float64, grid.N, num_converse), grid.Step+1);
    result[1] = case.initial_func(grid.x);
    Res::Function = FVS_residuals;
    if method == "FVS"
    elseif method == "Limiter"
        Res = generate_Limiter_residuals(Limiter; MUSCL_scheme=MUSCL_scheme, Sweby_beta=Sweby_beta);
    end
    for j = 1: grid.Step
        result[j+1] = runge_kutta4(case, grid, result[j], Res);
    end
    return result;
end

###################################################################################################

export UniformGrid1D;
export Case;
export solve;

end # module end