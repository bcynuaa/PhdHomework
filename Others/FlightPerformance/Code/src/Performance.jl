# author: bcynuaa
# date: 2022/10/22

module Performance # begin module

include("Atmosphere.jl");
include("DataReader.jl");
import .Atmosphere as atm;
using .DataReader;
using ExportAll;

using Dierckx;
using Roots;
using Optim;

default_len::Int64 = 100;
default_spline_order = 2;
default_spline_bc = "extrapolate"
alpha_fit_range::Vector{Float64} = [-4., 8.];
transonic_mach_range::Vector{Float64} = [0.2, 1.7];

colors = [
    :blue, :green, :red,
    :purple, :black, :orange,
    :pink, :darkblue, :wheat,
    :gray, :darkgreen, :navy,
    :rosybrown, :violet, :lightgreen,
    :yellow
];

function spline_minor(x, y)
    func = Spline1D(x, y; k=default_spline_order, bc=default_spline_bc);
    xx = LinRange(minimum(x), maximum(x), default_len);
    yy = func.(xx);
    return xx, yy;
end

###################################################################################################

function polar_curve_fit(CL::Vector{Float64}, CD::Vector{Float64})::Vector{Float64}
    N::Int64 = length(CL);
    sum_CL2::Float64 = sum(CL.^2);
    sum_CL4::Float64 = sum(CL.^4);
    lqA::Matrix{Float64} = [N sum_CL2; sum_CL2 sum_CL4];
    lqb::Vector{Float64} = [sum(CD), sum(CD .* CL.^2)];
    return lqA \ lqb;
end

function polar_curve_fit!(pfm::Dict)::Int64
    N::Int64 = length(pfm[:mach]);
    CD0::Vector{Float64} = zeros(N);
    A::Vector{Float64} = zeros(N);
    for i=1:N
        CD0[i], A[i] = polar_curve_fit(
            pfm[:CL_base][i, pfm[:alpha_fit_range_index]],
            pfm[:CD_base][i, pfm[:alpha_fit_range_index]]
        );
    end
    pfm[:CD0] = CD0;
    pfm[:A] = A;
    pfm[:A_mach_func] = Spline1D(pfm[:mach], pfm[:A]; k=default_spline_order, bc=default_spline_bc);
    pfm[:CD0_mach_func] = Spline1D(pfm[:mach], pfm[:CD0]; k=default_spline_order, bc=default_spline_bc);
    return 0
end

function polar_curve_function(CD0::Float64, A::Float64, CL::Float64)::Float64
    return CD0 + A * CL^2;
end

###################################################################################################

function Clalpha_fit(CL::Vector{Float64}, alpha::Vector{Float64}; is_rad::Int64=0)
    if is_rad == 0
        a = alpha .*pi / 180;
    else
        a = alpha;
    end
    N::Int64 = length(CL);
    sum_a::Float64 = sum(a);
    sum_a2::Float64 = sum(a.^2);
    lqA::Matrix{Float64} = [sum_a N; sum_a2 sum_a];
    lqb::Vector{Float64} = [sum(CL), sum(CL .* a)];
    Clalpha::Float64, b::Float64 = lqA \ lqb;
    alpha0 = - b / Clalpha;
    return Clalpha, alpha0;
end

function Clalpha_fit!(pfm::Dict)::Int64
    N::Int64 = length(pfm[:mach]);
    Clalpha::Vector{Float64} = zeros(N);
    alpha0::Vector{Float64} = zeros(N);
    for i=1:N
        Clalpha[i], alpha0[i] = Clalpha_fit(
            pfm[:CL_base][i, pfm[:alpha_fit_range_index]],
            pfm[:alpha_deg][pfm[:alpha_fit_range_index]];
            is_rad=0
        );
    end
    pfm[:Clalpha] = Clalpha;
    pfm[:alpha0] = alpha0;
    pfm[:Clalpha_func] = Spline1D(
        pfm[:mach], pfm[:Clalpha];
        k = default_spline_order,
        bc = default_spline_bc
    );
    return 0;
end

###################################################################################################

function Thrust_Require(h::Float64, ma::Float64, pfm::Dict)::Float64
    a::Float64 = atm.a(h);
    rho::Float64 = atm.rho(h);
    A::Float64 = pfm[:A_mach_func](ma);
    CD0::Float64 = pfm[:CD0_mach_func](ma);
    CL::Float64 = 2 * pfm[:W0] / rho / a^2 / ma^2 / pfm[:S_ref];
    CD::Float64 = polar_curve_function(CD0, A, CL);
    return CD * rho * ma^2 * a^2 * pfm[:S_ref] / 2;
end

function Thrust_Require!(pfm::Dict)::Int64
    function thrust_require_func(h::Float64, ma::Float64)
        return Thrust_Require(h, ma, pfm);
    end
    pfm[:thrust_require_func] = thrust_require_func;
    TRAB::Matrix{Float64} = pfm[:thrust_require_func].(pfm[:alt_m], transpose(pfm[:machthAB]));
    TR::Matrix{Float64} = pfm[:thrust_require_func].(pfm[:alt_m], transpose(pfm[:machth]));
    pfm[:thrustAB_require] = TRAB;
    pfm[:thrust_require] = TR;
    return 0;
end

function Thrust_Available!(pfm::Dict)::Int64
    Ta_func_list = []
    TaAB_func_list = []
    for j = 1:length(pfm[:alt_m])
        push!(Ta_func_list, Spline1D(pfm[:machth], pfm[:thrust][j, :]; k=default_spline_order, bc=default_spline_bc));
        push!(TaAB_func_list, Spline1D(pfm[:machthAB], pfm[:thrustAB][j, :]; k=default_spline_order, bc=default_spline_bc));
    end
    function Ta_func(h::Float64, ma::Float64)
        if h < minimum(pfm[:alt_m])
            return Ta_func_list[1](ma);
        elseif h > maximum(pfm[:alt_m])
            h1 = pfm[:alt_m][end];
            h2 = pfm[:alt_m][end-1];
            Ta1 = Ta_func_list[end](ma);
            Ta2 = Ta_func_list[end-1](ma);
            return (Ta1-Ta2) / (h1-h2) * (h-h1) + Ta1;
        elseif h in pfm[:alt_m]
            item = findfirst(x->x==h, pfm[:alt_m]);
            return Ta_func_list[item](ma);
        else
            item = findfirst(x->x>=h, pfm[:alt_m]);
            h1 = pfm[:alt_m][item];
            h2 = pfm[:alt_m][item-1];
            Ta1 = Ta_func_list[item](ma);
            Ta2 = Ta_func_list[item-1](ma);
            return (Ta2-Ta1) / (h2-h1) * (h-h1) + Ta1;
        end
    end
    pfm[:thrust_available_func] = Ta_func;
    function TaAB_func(h::Float64, ma::Float64)
        if h < minimum(pfm[:alt_m])
            return TaAB_func_list[1](ma);
        elseif h > maximum(pfm[:alt_m])
            h1 = pfm[:alt_m][end];
            h2 = pfm[:alt_m][end-1];
            Ta1 = TaAB_func_list[end](ma);
            Ta2 = TaAB_func_list[end-1](ma);
            return (Ta1-Ta2) / (h1-h2) * (h-h1) + Ta1;
        elseif h in pfm[:alt_m]
            item = findfirst(x->x==h, pfm[:alt_m]);
            return TaAB_func_list[item](ma);
        else
            item = findfirst(x->x>=h, pfm[:alt_m]);
            h1 = pfm[:alt_m][item];
            h2 = pfm[:alt_m][item-1];
            Ta1 = TaAB_func_list[item](ma);
            Ta2 = TaAB_func_list[item-1](ma);
            return (Ta2-Ta1) / (h2-h1) * (h-h1) + Ta1;
        end
    end
    pfm[:thrustAB_available_func] = TaAB_func;

    function Delta_T_func(h::Float64, ma::Float64)
        return pfm[:thrust_available_func](h, ma) - pfm[:thrust_require_func](h, ma);
    end
    pfm[:thrust_delta_func] = Delta_T_func;

    function Delta_TAB_func(h::Float64, ma::Float64)
        return pfm[:thrustAB_available_func](h, ma) - pfm[:thrust_require_func](h, ma);
    end
    pfm[:thrustAB_delta_func] = Delta_TAB_func;

    return 0;
end

###################################################################################################

function Speed_Range!(pfm::Dict)::Int64
    mach_range_alt::Matrix{Float64} = zeros(length(pfm[:alt_m]), 2);
    machAB_range_alt::Matrix{Float64} = zeros(length(pfm[:alt_m]), 2);
    mach_min::Float64 = 2e-3;
    mach_max::Float64 = 3.;
    zero_delta_T_mach_range = (mach_min, mach_max);

    for i = 1:length(pfm[:alt_m])
        function Delta_T_nAB(mach)
            return pfm[:thrust_delta_func](pfm[:alt_m][i], mach);
        end
        function Delta_T_AB(mach)
            return pfm[:thrustAB_delta_func](pfm[:alt_m][i], mach);
        end

        mach_range_alt[i, :] = find_zeros(Delta_T_nAB, zero_delta_T_mach_range);
        machAB_range_alt[i, :] = find_zeros(Delta_T_AB, zero_delta_T_mach_range);
    end

    pfm[:mach_range_alt] = mach_range_alt;
    pfm[:machAB_range_alt] = machAB_range_alt;

    pfm[:mach_min_alt_func] = Spline1D(pfm[:alt_m], mach_range_alt[:, 1]; k=default_spline_order, bc=default_spline_bc);
    pfm[:mach_max_alt_func] = Spline1D(pfm[:alt_m], mach_range_alt[:, 2]; k=default_spline_order, bc=default_spline_bc);

    pfm[:machAB_min_alt_func] = Spline1D(pfm[:alt_m], machAB_range_alt[:, 1]; k=default_spline_order, bc=default_spline_bc);
    pfm[:machAB_max_alt_func] = Spline1D(pfm[:alt_m], machAB_range_alt[:, 2]; k=default_spline_order, bc=default_spline_bc);
    return 0;
end

function Speed_OPT!(pfm::Dict)::Int64
    mach_opt::Vector{Float64} = zeros(length(pfm[:alt_m]));
    machAB_opt::Vector{Float64} = zeros(length(pfm[:alt_m]));
    # mach_opt_find0::Float64 = 2e-3;
    # mach_opt_find1::Float64 = 2.;
    for j = 1:length((pfm[:alt_m]))
        # f(mach) = -pfm[:thrust_delta_func](pfm[:alt_m][j], mach);
        # fAB(mach) = -pfm[:thrustAB_delta_func](pfm[:alt_m][j], mach);
        # mach_opt[j] = optimize(f, mach_opt_find0, mach_opt_find1).minimizer;
        # machAB_opt[j] = optimize(fAB, mach_opt_find0, mach_opt_find1).minimizer;
        f(mach) = pfm[:thrust_require_func](pfm[:alt_m][j], mach);
        fAB(mach) = pfm[:thrust_require_func](pfm[:alt_m][j], mach);
        mach_opt[j] = optimize(f, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer;
        machAB_opt[j] = optimize(fAB, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer;
    end

    pfm[:mach_alt_opt] = mach_opt;
    pfm[:machAB_alt_opt] = machAB_opt;

    pfm[:mach_opt_alt_func] = Spline1D(pfm[:alt_m], mach_opt; k=default_spline_order, bc=default_spline_bc);
    pfm[:machAB_opt_alt_func] = Spline1D(pfm[:alt_m], machAB_opt; k=default_spline_order, bc=default_spline_bc);
    return 0;
end

function Speed_gamma!(pfm::Dict)::Int64
    mach_gamma::Vector{Float64} = zeros(length(pfm[:alt_m]));
    machAB_gamma::Vector{Float64} = zeros(length(pfm[:alt_m]));
    # mach_opt_find0::Float64 = 2e-3;
    # mach_opt_find1::Float64 = 2.;
    for j = 1:length((pfm[:alt_m]))
        f(mach) = -pfm[:thrust_delta_func](pfm[:alt_m][j], mach);
        fAB(mach) = -pfm[:thrustAB_delta_func](pfm[:alt_m][j], mach);
        mach_gamma[j] = optimize(f, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer;
        machAB_gamma[j] = optimize(fAB, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer;
    end

    pfm[:mach_alt_gamma] = mach_gamma;
    pfm[:machAB_alt_gamma] = machAB_gamma;

    pfm[:mach_gamma_alt_func] = Spline1D(pfm[:alt_m], mach_gamma; k=default_spline_order, bc=default_spline_bc);
    pfm[:machAB_gamma_alt_func] = Spline1D(pfm[:alt_m], machAB_gamma; k=default_spline_order, bc=default_spline_bc);
    return 0;
end

function Speed_Vvmax!(pfm::Dict)::Int64
    function Vv(h::Float64, mach::Float64)::Float64
        return pfm[:thrust_delta_func](h, mach) * mach * atm.a(h) / pfm[:W0];
    end
    function VvAB(h::Float64, mach::Float64)::Float64
        return pfm[:thrustAB_delta_func](h, mach) * mach * atm.a(h) / pfm[:W0];
    end

    pfm[:Vv_alt_mach_func] = Vv;
    pfm[:VvAB_alt_mach_func] = VvAB;

    mach_opt_find0::Float64 = 0.2;
    mach_opt_find1::Float64 = 1.0;
    machAB_opt_find0::Float64 = 0.8;
    machAB_opt_find1::Float64 = 1.8;
    V_qc::Vector{Float64} = zeros(length(pfm[:alt_m]));
    VAB_qc::Vector{Float64} = zeros(length(pfm[:alt_m]));
    ma_qc::Vector{Float64} = zeros(length(pfm[:alt_m]));
    maAB_qc::Vector{Float64} = zeros(length(pfm[:alt_m]));
    for j = 1:length((pfm[:alt_m]))
        f(mach) = -pfm[:thrust_delta_func](pfm[:alt_m][j], mach) * mach * atm.a(pfm[:alt_m][j]);
        fAB(mach) = -pfm[:thrustAB_delta_func](pfm[:alt_m][j], mach) * mach * atm.a(pfm[:alt_m][j]);
        res1 = optimize(f, mach_opt_find0, mach_opt_find1);
        res2 = optimize(fAB, machAB_opt_find0, machAB_opt_find1);
        V_qc[j] = - res1.minimum / pfm[:W0];
        VAB_qc[j] = - res2.minimum / pfm[:W0];
        ma_qc[j] = res1.minimizer;
        maAB_qc[j] = res2.minimizer;
    end

    pfm[:Vvmax_alt] = V_qc;
    pfm[:VvmaxAB_alt] = VAB_qc;
    pfm[:mach_qc_alt] = ma_qc;
    pfm[:machAB_qc_alt] = maAB_qc;

    pfm[:Vvmax_alt_func] = Spline1D(pfm[:alt_m], V_qc; k=default_spline_order, bc=default_spline_bc);
    pfm[:VvmaxAB_alt_func] = Spline1D(pfm[:alt_m], VAB_qc; k=default_spline_order, bc=default_spline_bc);
    
    pfm[:mach_qc_alt_func] = Spline1D(pfm[:alt_m], ma_qc; k=default_spline_order, bc=default_spline_bc);
    pfm[:machAB_qc_alt_func] = Spline1D(pfm[:alt_m], maAB_qc; k=default_spline_order, bc=default_spline_bc);
    return 0;
end

###################################################################################################

function Hmax!(pfm::Dict)::Int64
    f0(height) = pfm[:Vvmax_alt_func](height);
    fAB0(height) = pfm[:VvmaxAB_alt_func](height);
    pfm[:Hmax0] = find_zero(f0, 15000.);
    pfm[:Hmax0AB] = find_zero(fAB0, 15000.);

    Vv_rest = 5.;
    f0_practice(height) = pfm[:Vvmax_alt_func](height) - Vv_rest;
    fAB0_practice(height) = pfm[:VvmaxAB_alt_func](height) - Vv_rest;
    pfm[:Hmax0_practice] = find_zero(f0_practice, 15000.);
    pfm[:Hmax0AB_practice] = find_zero(fAB0_practice, 15000.);
    return 0;
end

###################################################################################################

function cft!(pfm::Dict)::Int64;
    cft_f_list = [];
    for k = 1:length(pfm[:eta])
        for j = 1:length(pfm[:alt_m])
            push!(
                cft_f_list,
                Spline1D(
                    pfm[:machcft], pfm[:cft][k, j, :];
                    k = default_spline_order,
                    bc = default_spline_bc
                )
            );
        end
    end
    cft_f = reshape(
        cft_f_list, (length(pfm[:alt_m]), length(pfm[:eta]))
    );
    function cft_func_all(eta::Float64, h::Float64, mach::Float64)
        item1 = findfirst(x->x>=eta, pfm[:eta]);
        item2 = findfirst(x->x>=h, pfm[:alt_m]);
        if item1 == 1
            item1 += 1;
        elseif isnothing(item1) == true
            item1 = length(pfm[:eta]);
        end
        if item2 == 1
            item2 += 1;
        elseif isnothing(item2) == true
            item2 = length(pfm[:alt_m]);
        end
        mach3 = cft_f[item2, item1](mach);
        mach4 = cft_f[item2, item1-1](mach);
        mach2 = cft_f[item2-1, item1](mach);
        mach1 = cft_f[item2-1, item1-1](mach);
        mach12 = ( mach1 * (pfm[:eta][item1] - eta) + mach2 * (eta - pfm[:eta][item1-1]) ) / (pfm[:eta][item1] - pfm[:eta][item1-1]);
        mach34 = ( mach4 * (pfm[:eta][item1] - eta) + mach3 * (eta - pfm[:eta][item1-1]) ) / (pfm[:eta][item1] - pfm[:eta][item1-1]);
        mach1234 = ( mach12 * (pfm[:alt_m][item2] - h) + mach34 * (h - pfm[:alt_m][item2-1]) ) / (pfm[:alt_m][item2] - pfm[:alt_m][item2-1]);
        return mach1234;
    end
    function cft_func(h::Float64, mach::Float64)
        eta = pfm[:thrust_require_func](h, mach) / pfm[:thrust_available_func](h, mach);
        return cft_func_all(eta, h, mach);
    end

    function cfR_func_all(eta::Float64, h::Float64, mach::Float64)
        return cfR_func_all(eta, h, mach) / mach / atm.a(h);
    end
    function cfR_func(h::Float64, mach::Float64)
        return cft_func(h, mach) / mach / atm.a(h);
    end

    pfm[:cft_func] = cft_func;
    pfm[:cft_func_all] = cft_func_all;

    pfm[:cfR_func] = cfR_func;
    pfm[:cfR_func_all] = cfR_func_all;
    return 0;
end

function cftAB!(pfm::Dict)::Int64
    cft_f = [];
    for j = 1:length(pfm[:alt_m])
        push!(
            cft_f,
            Spline1D(
                pfm[:machcftAB], pfm[:cftAB][j, :];
                k = default_spline_order,
                bc = default_spline_bc
            );
        );
    end
    function cftAB_func(h::Float64, mach::Float64)
        if h <= minimum(pfm[:alt_m])
            return cft_f[1](mach);
        elseif h in pfm[:alt_m]
            item = findfirst(x->x==h, pfm[:alt_m]);
            return cft_f[item](mach);
        elseif h > maximum(pfm[:alt_m])
            cft1 = cft_f[end](mach);
            h1 = pfm[:alt_m][end];
            cft2 = cft_f[end-1](mach);
            h2 = pfm[:alt_m][end-1];
            return (cft1 - cft2) / (h1 - h2) * (h - h2) + cft2;
        else
            item = findfirst(x->x>=h, pfm[:alt_m]);
            cft1 = cft_f[item](mach);
            h1 = pfm[:alt_m][item];
            cft2 = cft_f[item-1](mach);
            h2 = pfm[:alt_m][item-1];
            return (cft1 - cft2) / (h1 - h2) * (h - h2) + cft2;
        end
    end

    function cfRAB_func(h::Float64, mach::Float64)
        return cftAB_func(h, mach) / mach / atm.a(h);
    end

    pfm[:cftAB_func] = cftAB_func;
    pfm[:cfRAB_func] = cfRAB_func;
    return 0;
end

function Speed_Long!(pfm::Dict)::Int64
    # mach_long::Vector{Float64} = zeros(length(pfm[:alt_m]));
    # mach_long_find0::Float64 = 0.2;
    # mach_long_find1::Float64 = 1.0;
    mach_long = [];
    for j = 1:length(pfm[:alt_m])
        if pfm[:alt_m][j] in [4000., 5000., 6000.]
            continue
        end
        f(mach) = pfm[:cft_func](pfm[:alt_m][j], mach);
        push!(
            mach_long,
            optimize(f, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer
        );
        # mach_long[j] = optimize(f, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer;
    end

    pfm[:mach_long_alt] = mach_long;
    pfm[:mach_long_alt_func] = Spline1D(
        vcat(pfm[:alt_m][1:4], pfm[:alt_m][8:end]),
        mach_long;
        k=default_spline_order, bc=default_spline_bc
    );
    return 0;
end

function Speed_Far!(pfm::Dict)::Int64
    # mach_far::Vector{Float64} = zeros(length(pfm[:alt_m]));
    # mach_far_find0::Float64 = 0.2;
    # mach_far_find1::Float64 = 1.0;
    mach_far = [];
    for j = 1:length(pfm[:alt_m])
        if pfm[:alt_m][j] in [4000., 5000., 6000., 7000.]
            continue
        end
        f(mach) = pfm[:cft_func](pfm[:alt_m][j], mach) / atm.a(pfm[:alt_m][j]) / mach;
        push!(
            mach_far,
            optimize(f, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer
        );
        # mach_far[j] = optimize(f, pfm[:mach_range_alt][j, 1], pfm[:mach_range_alt][j, 2]).minimizer;
    end

    pfm[:mach_far_alt] = mach_far;
    pfm[:mach_far_alt_func] = Spline1D(
        vcat(pfm[:alt_m][1:4], pfm[:alt_m][9:end]),
        mach_far;
        k=default_spline_order, bc=default_spline_bc
    );
    return 0;
end

###################################################################################################

function initial_performance(file_name)::Dict
    od::OriginData = OriginData(file_name)
    pfm::Dict = Dict(
        :mass_kg => od.mass_kg,
        :W0 => od.mass_kg * atm.g0,
        :S_ref => 1.,

        :alpha_deg => od.CL_base.x,
        :alpha_rad => od.CL_base.x .* (pi/180),
        :mach => od.CL_base.y,
        :CL_base => od.CL_base.data,
        :CD_base => od.CD_base.data,

        :machthAB => od.thrustAB.x,
        :alt_m => od.thrustAB.y,
        :thrustAB => od.thrustAB.data,
        :machth => od.thrust.x,
        :thrust => od.thrust.data,

        :cftAB => od.cftAB.data,
        :machcftAB => od.cft.x,
        
        :cft => od.cft.data,
        :eta => od.cft.z,
        :machcft => od.cft.x
    );

    pfm[:alpha_fit_range_index] = findfirst(x->x==alpha_fit_range[1], od.CL_base.x): findfirst(x->x==alpha_fit_range[2], od.CL_base.x);
    pfm[:alpha_fit_range] = LinRange(alpha_fit_range[1], alpha_fit_range[2], default_len);
    
    pfm[:transonic_mach_range_index] = findfirst(x->x==transonic_mach_range[1], od.CL_base.y): findfirst(x->x==transonic_mach_range[2], od.CL_base.y);
    pfm[:transonic_mach_range] = LinRange(transonic_mach_range[1], transonic_mach_range[2], default_len);

    pfm[:CL_range] = LinRange.( minimum(pfm[:CL_base][:, pfm[:alpha_fit_range_index]], dims=2), maximum(pfm[:CL_base][:, pfm[:alpha_fit_range_index]], dims=2), 100);
    # pfm[:thrustAB_available_func] = Spline2D(pfm[:alt_m], pfm[:machthAB], pfm[:thrustAB]; kx=default_spline_order, ky=default_spline_order);
    # pfm[:thrust_available_func] = Spline2D(pfm[:alt_m], pfm[:machth], pfm[:thrust]; kx=default_spline_order, ky=default_spline_order);

    polar_curve_fit!(pfm); # 极曲线计算
    Clalpha_fit!(pfm); # Clalpha 拟合计算
    Thrust_Require!(pfm); # 需用推力计算
    Thrust_Available!(pfm); # 可用推力插值
    Speed_Range!(pfm); # 速度平飞范围计算
    Speed_OPT!(pfm); # 有利速度计算
    Speed_gamma!(pfm); # 陡升速度计算
    Speed_Vvmax!(pfm); # 计算Vvmax
    Hmax!(pfm); # 计算升限
    cft!(pfm); # 计算cft
    cftAB!(pfm); # 计算加力cft
    Speed_Far!(pfm); # 远航速度
    Speed_Long!(pfm); # 久行速度
    return pfm;
end

###################################################################################################

@exportAll;
export Table, OriginData, atm, Dierckx;

end # end module