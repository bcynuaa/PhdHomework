# author: bcynuaa
# date: 2022/11/22

module EulerEq1D # module begin

###################################################################################################

const num_base::Int64 = 6;
const num_converse::Int64 = 3;

###################################################################################################

function base2converse(
    gamma::Float64,
    rho::Float64,
    u::Float64,
    p::Float64
)::Vector{Float64}
    rho_u::Float64 = rho * u;
    rho_E::Float64 = p / (gamma - 1.) + 1. / 2. * rho * u^2;
    return [rho, rho_u, rho_E];
end

function base2converse(
    gamma::Float64,
    vec_base::Vector{Float64}
)::Vector{Float64}
    return base2converse(gamma, vec_base[1], vec_base[2], vec_base[3]);
end

function base2converse(
    gamma::Float64,
    mat_base::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_base)[1];
    mat_converse::Matrix{Float64} = zeros(Float64, M, num_converse);
    for j = 1: M
        mat_converse[j, :] = base2converse(gamma, mat_base[j, :]);
    end
    return mat_converse;
end

###################################################################################################

function converse2base(
    gamma::Float64,
    rho::Float64,
    rho_u::Float64,
    rho_E::Float64
)::Vector{Float64}
    u::Float64 = rho_u / rho;
    E::Float64 = rho_E / rho;
    p::Float64 = (E - u^2 / 2.) * rho * (gamma - 1.);
    c::Float64 = sqrt(gamma * p / rho);
    H::Float64 = E + p / rho;
    return [rho, u, p, c, E, H];
    #      [  1, 2, 3, 4, 5, 6];
end

function converse2base(
    gamma::Float64,
    vec_converse::Vector{Float64}
)::Vector{Float64}
    return converse2base(gamma, vec_converse[1], vec_converse[2], vec_converse[3]);
end

function converse2base(
    gamma::Float64,
    mat_converse::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_converse)[1];
    mat_base::Matrix{Float64} = zeros(Float64, M, num_base);
    for j = 1: M
        mat_base[j, :] = converse2base(gamma, mat_converse[j, :]);
    end
    return mat_base;
end

###################################################################################################

function F_converse(
    gamma::Float64,
    rho::Float64,
    rho_u::Float64 = Float64,
    rho_E::Float64 = Float64
)::Vector{Float64}
    rho::Float64, u::Float64, p::Float64, c::Float64, E::Float64, H::Float64 = converse2base(gamma, rho, rho_u, rho_E);
    return [
        rho * u,
        rho * u^2 + p,
        rho * u * H
    ];
end

function F_converse(
    gamma::Float64,
    vec_con::Vector{Float64}
)::Vector{Float64}
    return F_converse(gamma, vec_con[1], vec_con[2], vec_con[3]);
end

function F_converse(
    gamma::Float64,
    mat_con::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_con)[1];
    F::Matrix{Float64} = zeros(Float64, M, num_converse);
    for j = 1: M
        F[j, :] = F_converse(gamma, mat_con[j, :]);
    end
    return F;
end

###################################################################################################

function F_base(
    gamma::Float64,
    rho::Float64,
    u::Float64,
    p::Float64
)::Vector{Float64}
    return F_converse(gamma, base2converse(gamma, rho, u, p));
end

function F_base(
    gamma::Float64,
    vec_base::Vector{Float64}
)::Vector{Float64}
    return F_base(gamma, vec_base[1], vec_base[2], vec_base[3]);
end

function F_base(
    gamma::Float64,
    mat_base::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_base)[1];
    F::Matrix{Float64} = zeros(Float64, M, num_converse);
    for j = 1: M
        F[j, :] = F_base(gamma, mat_base[j, :]);
    end
    return F;
end

###################################################################################################

function F_plus_converse(
    gamma::Float64,
    rho::Float64,
    rho_u::Float64,
    rho_E::Float64
)::Vector{Float64}
    rho::Float64, u::Float64, p::Float64, c::Float64, E::Float64, H::Float64 = converse2base(gamma, rho, rho_u, rho_E);
    Ma::Float64 = u / c;
    return rho * c / 4. * (Ma + 1.)^2 * [
        1.,
        2. * c / gamma * ( 1. + (gamma - 1.) / 2. * Ma ),
        2. * c^2 / (gamma^2 - 1.) * ( 1. + (gamma - 1.) / 2. * Ma )^2
    ];
end

function F_plus_converse(
    gamma::Float64,
    vec_con::Vector{Float64}
)::Vector{Float64}
    return F_plus_converse(gamma, vec_con[1], vec_con[2], vec_con[3]);
end

function F_plus_converse(
    gamma::Float64,
    mat_con::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_con)[1];
    F_plus::Matrix{Float64} = zeros(Float64, M, num_converse);
    for j = 1: M
        F_plus[j, :] = F_plus_converse(gamma, mat_con[j, :]);
    end
    return F_plus;
end

###################################################################################################

function F_plus_base(
    gamma::Float64,
    rho::Float64,
    u::Float64,
    p::Float64
)::Vector{Float64}
    return F_plus_converse(gamma, base2converse(gamma, rho, u, p));
end

function F_plus_base(
    gamma::Float64,
    vec_base::Vector{Float64}
)::Vector{Float64}
    return F_plus_base(gamma, vec_base[1], vec_base[2], vec_base[3]);
end

function F_plus_base(
    gamma::Float64,
    mat_base::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_base)[1];
    F_plus::Matrix{Float64} = zeros(Float64, M, num_converse);
    for j = 1: M
        F_plus[j, :] = F_plus_base(gamma, mat_base[j, :]);
    end
    return F_plus;
end

###################################################################################################

function F_minus_converse(
    gamma::Float64,
    rho::Float64,
    rho_u::Float64,
    rho_E::Float64
)::Vector{Float64}
    rho::Float64, u::Float64, p::Float64, c::Float64, E::Float64, H::Float64 = converse2base(gamma, rho, rho_u, rho_E);
    Ma::Float64 = u / c;
    return -rho * c / 4. * (Ma - 1.)^2 * [
        1.,
        2. * c / gamma * ( -1. + (gamma - 1.) / 2. * Ma ),
        2. * c^2 / (gamma^2 - 1.) * ( -1. + (gamma - 1.) / 2. * Ma )^2
    ];
end

function F_minus_converse(
    gamma::Float64,
    vec_con::Vector{Float64}
)::Vector{Float64}
    return F_minus_converse(gamma, vec_con[1], vec_con[2], vec_con[3]);
end

function F_minus_converse(
    gamma::Float64,
    mat_con::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_con)[1];
    F_minus::Matrix{Float64} = zeros(Float64, M, num_converse);
    for j = 1: M
        F_minus[j, :] = F_minus_converse(gamma, mat_con[j, :]);
    end
    return F_minus;
end

###################################################################################################

function F_minus_base(
    gamma::Float64,
    rho::Float64,
    u::Float64,
    p::Float64
)::Vector{Float64}
    return F_minus_converse(gamma, base2converse(gamma, rho, u, p));
end

function F_minus_base(
    gamma::Float64,
    vec_base::Vector{Float64}
)::Vector{Float64}
    return F_minus_base(gamma, vec_base[1], vec_base[2], vec_base[3]);
end

function F_minus_base(
    gamma::Float64,
    mat_base::Matrix{Float64}
)::Matrix{Float64}
    M::Int64 = size(mat_base)[1];
    F_minus::Matrix{Float64} = zeros(Float64, M, num_converse);
    for j = 1: M
        F_minus[j, :] = F_minus_base(gamma, mat_base[j, :]);
    end
    return F_minus;
end
###################################################################################################

export num_base, num_converse;
export base2converse, converse2base;
export F_base, F_plus_base, F_minus_base;
export F_converse, F_plus_converse, F_minus_converse;

end # module end