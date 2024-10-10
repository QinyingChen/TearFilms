include("spectral.jl")
using DifferentialEquations, JLD2

# parameterize the advection and diffusion coefficients using point values on a grid

n = 32  # space resolution
m = 60  # time resolution
N = 500   # number of samples

plan = plan_fderiv((n, n))

# heat equation
function timederiv(u, k, t)
    # k, c = p
    U = reshape(u, n, n)
    ux = fderiv(U, 1, 1, plan)
    uy = fderiv(U, 2, 1, plan)
    return vec(fderiv(k.x .* ux, 1, 1, plan) + fderiv(k.y .* uy, 2, 1, plan))
end
u₀ = (x, y) -> exp(cospi(2x) + sinpi(2y))

t = range(0, 1, length=m+1)
x = collect(2*(0:n-1)) / n
U = zeros(Float32, n, n, m + 1, N)     # store solutions
Coeff_kx = zeros(Float32, n, n, N)   # store coefficients
Coeff_ky = zeros(Float32, n, n, N)   # store coefficients
for idx in 1:N
    kx = [-1]
    # random diffusion
    while any(kx .< 0.05)
        kx = fill(0.5, n, n)
        for mx in 1:6, my in 1:6
            kx .+= 0.1*randn()/hypot(mx,my) * sinpi.(mx*x) .* sinpi.(my*x')
            kx .+= 0.1*randn()/hypot(mx,my) * cospi.(mx*x) .* sinpi.(my*x')
            kx .+= 0.1*randn()/hypot(mx,my) * sinpi.(mx*x) .* cospi.(my*x')
            kx .+= 0.1*randn()/hypot(mx,my) * cospi.(mx*x) .* cospi.(my*x')
        end
    end
    ky = [-1]
    while any(ky .< 0.05)
        ky = fill(0.5, n, n)
        for mx in 1:6, my in 1:6
            ky .+= 0.2*randn()/hypot(mx,my) * sinpi.(mx*x) .* sinpi.(my*x')
        end
    end
    Coeff_kx[:, :, idx] .= kx
    Coeff_ky[:, :, idx] .= ky
    init = u₀.(x, x')
    ivp = ODEProblem(timederiv, vec(init), (0.0, 1.0), (;x=kx, y=ky))
    sol = solve(ivp, QNDF(autodiff=false), saveat=t)
    # store solution
    for i in 1:m+1
        U[:, :, i, idx] .= reshape(sol.u[i], n, n)
    end
    @save "heat2d_data3.jld2" Coeff_kx Coeff_ky U x t
end

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
