# author: bcynuaa
# date: 2022/10/22

module Atmosphere # module begin

const gamma::Float64 = 1.4; # ratio of specific heat
const R::Float64 = 287.06; # gas constant
const T_absolute::Float64 = -273.15 # absolute zero
const T_h::Float64 = -56.5 - T_absolute; # temperature in stratosphere

# standard parameters (0 altitude)
const T0::Float64 = 15. - T_absolute;
const p0::Float64 = 1013.25e2;
const rho0::Float64 = 1.225;
const c0::Float64 = 340.294;
const g0::Float64 = 9.80665;

# heat descending per meter
const lambda::Float64 = 6.5 / 1000;
const ratio::Float64 = 5.494e-4 / -3.567e-3;
const lambda2::Float64 = lambda * ratio;

const h1::Float64 = 11e3; # troposphere height
const h2::Float64 = 20e3; # tropopause height
const h3::Float64 = 27432.; # stratosphere height

const T11::Float64 = 216.65; # temperature at troposphere height
const T20::Float64 = 216.65; # temperature at tropopause height
const p11::Float64 = 226.32e2; # pressure at trotosphere height
const p20::Float64 = 5475.06; # pressure at tropopause height
# const rho11::Float64 = 0.36392; # rho at trotosphere height

function T(h::Real)::Float64
    if h < 0.
        return T0;
    elseif 0<= h <= h1
        return T0 - lambda * h;
    elseif h1 < h <= h2
        return T_h;
    elseif h2 < h <= h3
        return T20 - lambda2 * (h - h2);
    else
        return T(h3);
    end 
end

function p(h::Real)::Float64
    if h < 0.
        return p0;
    elseif 0<= h <= h1
        return p0 * (T(h) / T0) ^ (g0 / R / lambda);
    elseif h1 < h <= h2
        return p11 * exp( -g0 * (h - h1) / R / T11 );
    elseif h2 < h <= h3
        return p20 * (T(h) / T20) ^ (g0 / R / lambda2);
    else
        return p(h3);
    end
end

function rho(h::Real)::Float64
    return p(h) / T(h) / R;
end

function a(h::Real)::Float64
    return air_speed(T(h));
end

function air_speed(T::Real)::Float64
    return sqrt(gamma * R * T);
end

function air_speed(p::Real, rho::Real)::Float64
    return sqrt(gamma * p / rho);
end

end # module end