include("spectral.jl")
using DifferentialEquations, JLD2

# parameterize the diffusion coefficient using point values on a grid
# zero Dirichlet boundary conditions

n = 64  # space resolution
m = 80  # time resolution
N = 5000   # number of samples

D, x = cheb(n)
t = range(0, 1, length=m+1)
D1 = D[:, 2:end-1]
D2 = D[2:end-1, :]

# advection-diffusion equation
function timederiv(u, p, t)
    k, c = p
    uₓ = D1 * u
    return -c .* uₓ[2:end-1] + D2 * (k .* uₓ)
end
u₀ = x -> (1-x^4) * exp(cospi(3*x))

U = zeros(Float32, n - 1, m + 1, N)     # store solutions
Coeff_k = zeros(Float32, n + 1, N)   # store coefficients
Coeff_c = zeros(Float32, n - 1, N)   # store coefficients
for idx in 1:N
    θ = acos.(x)
    k = [-1]
    # random diffusion
    while any(k .< 0.05)
        k = 0.5 .+ sum(cos.(m.*θ) .* 0.1*randn()/m for m in 1:10)
    end
    θ = acos.(x̂)
    c = [-1]
    # random advection
    while any(c .< 0.05)
        c = 2 .+ sum(cos.(m.*θ) .* 0.5*randn()/m for m in 1:10)
    end
    Coeff_k[:, idx] .= k
    Coeff_c[:, idx] .= c
    ivp = ODEProblem(timederiv, u₀.(x̂), (0.0, 1.0), (k, c))
    sol = solve(ivp, QNDF(autodiff=false), saveat=t)
    # store solution
    for i in 1:m+1
        U[:, i, idx] .= sol.u[i]
    end
end

@save "advdiff1d_bc_data2.jld2" Coeff_k Coeff_c U
##

# # @load "heat1d_data2.jld2" coeff_k U
# # include("multi_pca.jl")
# # pca_sol = MultiPCA(U, rank=30, eps=0.01)
# # V = transform(pca_sol, U);
# # pca_coeff = MultiPCA(diff_coeff, rank=16, eps=0.01)
# # C = transform(pca_coeff, diff_coeff);
# # @save "heat1d_pca2.jld2" pca_sol pca_coeff V C

# using HDF5
# h5open("heat1d_pca2.h5", "w") do file
#     write(file, "V", V)
#     write(file, "C", C)
# end
