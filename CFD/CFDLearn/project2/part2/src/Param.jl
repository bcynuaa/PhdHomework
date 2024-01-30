# author: bcynuaaa_inf
# date: 2022/10/04

module Param # module begin

using ExportAll;

###################################################################################################

"""
# Case::DataType

- `Ma_inf::Float64` -> # speed of incoming flow
- `theta::Float64` -> angle of coming flow

- `grid_file::String` -> the grid file
- `data_path::String` -> the data saving path
- `data_name::String` -> the data file name

- `gamma::Float64` -> ratio of specific heat
- `R::Float64` -> gas constant
- `Cv::Float64` -> gas heat capacity

- `rho_inf::Float64` -> gas density
- `P_inf::Float64` -> gas pressure
- `T_inf::Float64` -> gas temperature
- `C_inf::Float64` -> gas speed

- `k2::Float64` -> coefficient chosen for 2nd order artifical viscosity
- `k4::Float64` -> coefficient chosen for 4th order artifical viscosity
- `CFL::Float64` -> number of CFL
- `error::Float64` -> the tolerant error

- `STEP::Int64` -> total step
    
- `w_inf::Vector{Float64}` -> conservative value of w
- `w_inf_expand::Vector{Float64}` -> expand the conservative value of w
"""
struct Case
    Ma_inf::Float64; # speed of incoming flow
    theta::Float64; # angle of coming flow

    grid_file::String; # the grid file
    data_path::String; # the data saving path
    data_name::String; # the data file name

    gamma::Float64; # ratio of specific heat
    R::Float64; # gas constant
    Cv::Float64; # gas heat capacity

    rho_inf::Float64; # gas density
    P_inf::Float64; # gas pressure
    T_inf::Float64; # gas temperature
    C_inf::Float64; # gas speed

    k2::Float64; # coefficient chosen for 2nd order artifical viscosity
    k4::Float64; # coefficient chosen for 4th order artifical viscosity
    CFL::Float64; # number of CFL
    error::Float64; # the tolerant error

    STEP::Int64; # total step
    
    w_inf::Vector{Float64}; # conservative value of w
    w_inf_expand::Vector{Float64}; # expand the conservative value of w
end

###################################################################################################

"""
# expand_w(gamma::Float64, w::Vector{Float64})::Vector{Float64}

**to attain variables from conservative value**

- `gamma::Float64` -> γ of the gas
- `w::Vector{Float64}` -> the vector of [ρ, ρU, ρV, ρE]

return

- `[ρ, U, V, E, P, H, C]::Vector{Float64}`

   [1, 2, 3, 4, 5, 6, 7]
"""
function expand_w(gamma::Float64, w::Vector{Float64})::Vector{Float64}
    rho::Float64 = w[1];
    U::Float64 = w[2] / rho;
    V::Float64 = w[3] / rho;
    E::Float64 = w[4] / rho;
    P::Float64 = (E - (U^2 + V^2)/2.) * (gamma - 1.) * rho;
    H::Float64 = E + P/rho;
    C::Float64 = sqrt(gamma * P / rho);
    return [rho, U, V, E, P, H, C];
    #         1, 2, 3, 4, 5, 6, 7
end

"""
# expand_w(case::Case, w::Vector{Float64})::Vector{Float64}

**to attain variables from conservative value**

- `case::Case` -> the Case struct
- `w::Vector{Float64}` -> the vector of [ρ, ρU, ρV, ρE]

return

- `[ρ, U, V, E, P, H, C]::Vector{Float64}`

   [1, 2, 3, 4, 5, 6, 7]
"""
function expand_w(case::Case, w::Vector{Float64})::Vector{Float64}
    return expand_w(case.gamma, w);
end

"""
# expand_W(gamma::Float64, W::Matrix{Float64})::Matrix{Float64}

**to attain variables from conservative value on all cells**

- `gamma::Float64` -> γ of the gas
- `w::Matrix{Float64}` -> the Matrix of all [ρ, ρU, ρV, ρE]

return:

all

- `[ρ, U, V, E, P, H, C]::Vector{Float64}`

   [1, 2, 3, 4, 5, 6, 7]
"""
function expand_W(gamma::Float64, W::Matrix{Float64})::Matrix{Float64}
    N::Int64 = size(W)[1];
    W_expand::Matrix{Float64} = zeros(Float64, N, 7);
    for j = 1: N
        # println("$j  $(W[j, :])");
        W_expand[j, :] = expand_w(gamma, W[j, :]);
    end
    return W_expand;             
end

"""
# expand_W(case::Case, W::Matrix{Float64})::Matrix{Float64}

**to attain variables from conservative value on all cells**

- `case::Case` -> the Case struct
- `w::Matrix{Float64}` -> the Matrix of all [ρ, ρU, ρV, ρE]

return:

all

- `[ρ, U, V, E, P, H, C]::Vector{Float64}`

   [1, 2, 3, 4, 5, 6, 7]
"""
function expand_W(case::Case, W::Matrix{Float64})::Matrix{Float64}
    return expand_W(case.gamma, W);         
end

###################################################################################################

"""
# Case(Ma_inf::Float64, theta::Float64; ...)

**get a Case by specifying parameters**

return

`case::Case`
"""
function Case(
    Ma_inf::Float64, theta::Float64;
    grid_file::String = "naca0012.grd",
    data_path::String = ".//",
    data_name::String = "result.dat",
    gamma::Float64 = 1.4,
    R::Float64 = 1. / 1.4,
    Cv::Float64 = 2.5/1.4,
    rho_inf::Float64 = 1.4, P_inf::Float64 = 1.0,
    T_inf::Float64 = 1.0, C_inf::Float64 = 1.0,
    k2::Float64 = 0.9, k4::Float64 = 0.0255,
    CFL::Float64 = 2.0, error::Float64 = 1e-5,
    STEP::Int64 = 10000
)::Case
    U_inf::Float64 = C_inf * Ma_inf * cos(theta);
    V_inf::Float64 = C_inf * Ma_inf * sin(theta);
    T_inf::Float64 = T_inf;
    E_inf::Float64 = Cv * T_inf + (U_inf^2 + V_inf^2)/2.;
    w_inf::Vector{Float64} = [rho_inf, rho_inf * U_inf, rho_inf * V_inf, rho_inf * E_inf];
    w_inf_expand::Vector{Float64} = expand_w(gamma, w_inf);
    return Case(
        Ma_inf, theta,
        grid_file, data_path, data_name,
        gamma, R, Cv,
        rho_inf, P_inf, T_inf, C_inf,
        k2, k4,
        CFL, error, STEP,
        w_inf, w_inf_expand
    );
end

###################################################################################################

@exportAll;

end # module end