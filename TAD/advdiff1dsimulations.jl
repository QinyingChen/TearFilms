include("spectral.jl")
using DifferentialEquations, JLD2

# parameterize the diffusion coefficient using point values on a grid

n = 64  # space resolution
m = 80  # time resolution
N = 5000   # number of samples

F, F⁻¹, diff_mult1 = plan_fderiv(n)

# heat equation
function timederiv(u, p, t)
    k, c = p
    uₓ = fderiv(u, F, F⁻¹, diff_mult1)
    return -c .* uₓ + fderiv(k .* uₓ, F, F⁻¹, diff_mult1)
end
u₀ = x -> exp(sinpi(2x))

x = collect(2*(0:n-1)) / n
U = zeros(Float32, n, m + 1, N)     # store solutions
Coeff_k = zeros(Float32, n, N)   # store coefficients
Coeff_c = zeros(Float32, n, N)   # store coefficients
for idx in 1:N
    k = [-1]
    # random diffusion
    while any(k .< 0.05)
        k = 0.5 .+ sum(sinpi.(m.*x) .* 0.2*randn()/m^2 for m in 1:10)
    end
    c = [-1]
    # random advection
    while any(c .< 0.05)
        c = 2 .+ sum(sinpi.(m.*x) .* 1*randn()/m^2 for m in 1:10)
    end
    Coeff_k[:, idx] .= k
    Coeff_c[:, idx] .= c
    ivp = ODEProblem(timederiv, u₀.(x), (0.0, 1.0), (k, c))
    sol = solve(ivp, QNDF(autodiff=false), saveat=range(0, 1, length=m+1))
    # store solution
    for i in 1:m+1
        U[:, i, idx] .= sol.u[i]
    end
end

@save "advdiff1d_data2.jld2" Coeff_k Coeff_c U
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
