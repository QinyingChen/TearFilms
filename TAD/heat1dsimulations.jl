include("spectral.jl")
using DifferentialEquations, JLD2

# parameterize the diffusion coefficient using point values on a grid

n = 64  # space resolution
m = 80  # time resolution
N = 5000   # number of samples

# heat equation
function timederiv(u, p, t)
    return fderiv(p .* fderiv(u))
end
u₀ = x -> exp(sinpi(2x))

x = collect(2*(0:n-1)) / n
U = zeros(Float32, n, m + 1, N)     # store solutions
diff_coeff = zeros(Float32, n, N)   # store coefficients
for idx in 1:N
    coeff = [-1]
    # random coefficient
    while any(coeff .< 0)
        coeff = 1 .+ sum(sinpi.(k.*x) .* 0.4*randn()/k^2 for k in 1:15)
    end
    diff_coeff[:, idx] .= coeff
    ivp = ODEProblem(timederiv, u₀.(x), (0.0, 1.0), coeff)
    sol = solve(ivp, QNDF(autodiff=false), saveat=range(0, 1, length=m+1))
    # store solution
    for i in 1:m+1
        U[:, i, idx] .= sol.u[i]
    end
end

@save "heat1d_data.jld2" diff_coeff U
##

@load "heat1d_data.jld2" diff_coeff U
include("multi_pca.jl")
pca_sol = MultiPCA(U, rank=30, eps=0.01)
V = transform(pca_sol, U);
pca_coeff = MultiPCA(diff_coeff, rank=16, eps=0.01)
C = transform(pca_coeff, diff_coeff);
@save "heat1d_pca.jld2" pca_sol pca_coeff V C

# using HDF5
# h5open("heat1d_pca.h5", "w") do file
#     write(file, "V", V)
#     write(file, "C", C)
# end
