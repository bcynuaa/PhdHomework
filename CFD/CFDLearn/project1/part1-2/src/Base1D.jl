# auther: bcynuaa
# date: 2022/10/12

###################################################################################################

using LinearAlgebra;
using SparseArrays;

###################################################################################################

module Base1D # module begin

using LinearAlgebra;
using SparseArrays;
using ExportAll;

###################################################################################################

struct Grid
    x::Vector{Float64}
    t::Vector{Float64}
    Nx::Int64
    Nt::Int64
end

function Grid(x::Vector{Float64}, t::Vector{Float64})::Grid
    return Grid(x, t, length(x), length(t));
end

function Grid(x_begin::Float64, x_end::Float64, dx::Float64, dt::Float64, Step::Int64)::Grid
    return Grid(collect(x_begin: dx: x_end), [i*dt for i=0:Step-1]);
end

###################################################################################################

function interpolation_coefficient3(x::Vector{Float64}, i::Int64; flag::Int64=0)
    d1::Float64 = x[i] - x[i-1];
    d2::Float64 = x[i+1] - x[i];
    if flag==0
        return [
            -d2 / (d1^2 + d1*d2),
            (d2 - d1) / d1 / d2,
            d1 / (d1*d2 + d2^2)
        ];
    elseif flag==-1
        return [
            (-2*d1-d2) / (d1^2 + d1*d2),
            (d1+d2) / d1 / d2,
            -d1 / (d1*d2 + d2^2)
        ];
    else
        return [
            d2 / (d1^2 + d1*d2),
            (-d1-d2) / d1 / d2,
            (d1+2*d2) / (d1*d2 + d2^2)
        ];
    end
end

function interpolation_coefficient3_order2(x::Vector{Float64}, i::Int64)
    d1::Float64 = x[i] - x[i-1];
    d2::Float64 = x[i+1] - x[i];
    return [
        2/d1/(d1+d2),
        -2/d1/d2,
        2/d2/(d1+d2)
    ];
end

###################################################################################################

function backward_SparseMatrix(grid::Grid)::SparseMatrixCSC{Float64, Int64}
    d1::Vector{Float64} = zeros(grid.Nx - 1);
    d0::Vector{Float64} = zeros(grid.Nx);
    for i = 2:grid.Nx
        d1[i-1] = - 1/(grid.x[i] - grid.x[i-1]);
        d0[i] = -d1[i-1];
    end
    d0[1] = -1 / (grid.x[2] - grid.x[1]);
    A::SparseMatrixCSC{Float64, Int64} = spdiagm(-1=>d1) + spdiagm(d0);
    A[1, 2] = -A[1, 1];
    #A[1, 1:3] = interpolation_coefficient3(grid.x, 2; flag=-1);
    return A;
end

function central_SparseMatrix(grid::Grid)::SparseMatrixCSC{Float64, Int64}
    d1::Vector{Float64} = zeros(grid.Nx - 1);
    d0::Vector{Float64} = zeros(grid.Nx);
    d_1::Vector{Float64} = zeros(grid.Nx - 1);
    for i = 2:grid.Nx-1
        d_1[i-1], d0[i], d1[i] = interpolation_coefficient3(grid.x, i);
    end
    A::SparseMatrixCSC{Float64, Int64} = spdiagm(-1=>d_1, 1=>d1) + spdiagm(d0);
    A[1, 1:3] = interpolation_coefficient3(grid.x, 2; flag=-1);
    A[end, end-2:end] = interpolation_coefficient3(grid.x, grid.Nx-1; flag=1);
    return A;
end

function forward_SparseMatrix(grid::Grid)::SparseMatrixCSC{Float64, Int64}
    d1::Vector{Float64} = zeros(grid.Nx - 1);
    d0::Vector{Float64} = zeros(grid.Nx);
    for i = 1:grid.Nx-1
        d1[i] = 1/(grid.x[i+1] - grid.x[i]);
        d0[i] = -d1[i];
    end
    d0[end] = 1 / (grid.x[end] - grid.x[end-1]);
    A::SparseMatrixCSC{Float64, Int64} = spdiagm(1=>d1) + spdiagm(d0);
    A[end, end-1] = -A[end, end];
    return A;
end

function central_SparseMatrix_order2(grid::Grid)::SparseMatrixCSC{Float64, Int64}
    d1::Vector{Float64} = zeros(grid.Nx - 1);
    d0::Vector{Float64} = zeros(grid.Nx);
    d_1::Vector{Float64} = zeros(grid.Nx - 1);
    for i = 2:grid.Nx-1
        d_1[i-1], d0[i], d1[i] = interpolation_coefficient3_order2(grid.x, i);
    end
    A::SparseMatrixCSC{Float64, Int64} = spdiagm(-1=>d_1, 1=>d1) + spdiagm(d0);
    A[1, 1:3] = interpolation_coefficient3_order2(grid.x, 2);
    A[end, end-2:end] = interpolation_coefficient3_order2(grid.x, grid.Nx-1);
    return A;
end

###################################################################################################

@exportAll();

end # mdule end