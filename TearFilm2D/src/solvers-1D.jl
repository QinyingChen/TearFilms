# Fourier discretization over (-π, π]
struct Grid1D
    t::AbstractVector{Float64}
    D::AbstractMatrix{Float64}
    DD::AbstractMatrix{Float64}
end

function Grid1D(m, L=1)
    @assert iseven(m)
    n = m ÷ 2
    x = (1-n:n) * π / n
    entry(k) = k == 0 ? 0.0 : 0.5 * (-1)^k * cot(k * π / m)
    Dx = [entry(mod(i - j, m)) for i in 1:m, j in 1:m]
    entry2(k) = k == 0 ? -π^2 / 3(2π / m)^2 - 1 / 6 : -(-1)^k / (2 * (sin(k * π / m))^2)
    Dxx = [entry2(mod(i - j, m)) for i in 1:m, j in 1:m]
    return Grid1D(L * x, Dx / L, Dxx / L^2)
end

# Fourier discretization over [0, π] with even/odd parity
struct Grid1DSym
    t::AbstractVector{Float64}
    D_even::AbstractMatrix{Float64}
    D_odd::AbstractMatrix{Float64}
    DD::AbstractMatrix{Float64}
end

function Grid1DSym(m, L=1)
    @assert iseven(m)
    n = m ÷ 2
    g = Grid1D(m, L)
    x = g.t[n:end]

    # When differentiating an even function, leave out the end values,
    # because they are zero.
    Dx_even = g.D[n+1:m-1, n:m]
    Dx_even[:, 2:n] .+= g.D[n+1:m-1, n-1:-1:1]

    # When differentiating an odd function, restore the end values.
    Dx_odd = g.D[n:m, n+1:m-1]
    Dx_odd .-= g.D[n:m, n-1:-1:1]

    Dxx = g.DD[n:m, n:m]
    Dxx[:, 2:n] .+= g.DD[n:m, n-1:-1:1]
    return Grid1DSym(x, Dx_even, Dx_odd, Dxx)
end

# Compute the radial Laplacian
function laplacian(u, grid)
    lapu = (grid.D_odd * (grid.t[2:end-1] .* u)) ./ grid.t
    lapu[1] = 2 * dot(grid.DD[1,2:end-1], u)
    return lapu
end

# time derivative function for the radial PDE
function radial_ode!(dv, v, params, t)
    Je, index, grid, constants = params
    @unpack ell, Pc, invPec = constants
    # use views to avoid allocation and copying
    h = view(v, index.h)
    p = view(v, index.p)
    c = view(v, index.c)

    c_x = grid.D_even * c
    ubar = ((-h[2:end-1] .^ 2 / 12) .* (grid.D_even * p))

    h_rr = laplacian(grid.D_even * h, grid)
    tmp = laplacian(h[2:end-1] .* ubar, grid)
    tmp2 = laplacian(h[2:end-1] .* c_x, grid)

    @. dv[index.h] = Pc * (c - 1) - tmp - Je
    @. dv[index.p] = -h_rr - p
    @. dv[index.c] = (invPec * tmp2 - Pc * (c - 1) * c + Je * c) / h
    @. dv[index.c[2:end-1]] -= ubar * c_x
    return dv
end

function radial_solve(n, tspan, rw=0.5, vb=0.1,d=4.5e-6)
    # discretization parameters
    grid = Grid1DSym(2n, sqrt(2))
    M = Diagonal([ones(n+1); zeros(n+1); ones(n+1)])
    u0 = [ones(n+1); zeros(n+1); ones(n+1)]
    index = (h = 1:n+1, p = n+2:2(n+1), c = 2(n+1)+1:3(n+1))

    # physical parameters
    consts = PhysicalConstants()
    pde_consts = PDEConstants(consts, 1e-6 / (60vb),d)

    # evaporative flux
    hump = @. exp(-(grid.t / rw)^2 / 2)
    Je =  vb .+ (1 - vb) * hump

    # solve
    params = (Je, index, grid, pde_consts)
    IVP = ODEProblem(ODEFunction(radial_ode!, mass_matrix=M), u0, tspan, params)
    return solve(IVP, reltol=1e-6, abstol=1e-6)
end
