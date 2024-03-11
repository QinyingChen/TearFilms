
function pack!(u, h, p, c)
    m, n = size(h)
    mn = m*n
    u[1:mn] .= vec(h)
    u[mn+1:2mn] .= vec(p)
    u[2mn+1:3mn] .= vec(c)
    return u
end

pack(h, p, c) = pack!(similar(h,(3*length(h),)), h, p, c)

function unpack!(h, p, c, u, index)
    h .= reshape(u[index.h], index.size)
    p .= reshape(u[index.p], index.size)
    c .= reshape(u[index.c], index.size)
    return h, p, c
end

function unpack(u, index)
    h = reshape(view(u, index.h), index.size )
    p = reshape(view(u, index.p), index.size )
    c = reshape(view(u, index.c), index.size )
    return h, p, c
end

function twodim_ode!(du, u, params, t)
    Je, index, gridx, gridy, constants = params
    @unpack Pc, invPec = constants
    h, p, c = unpack(u, index)

    tmp = -h .^ 2 / 12
    ubar = tmp .* (gridx.D * p)
    vbar = tmp .* (p * gridy.D')
    c_x = gridx.D * c
    c_y = c * gridy.D'

    osmo = @. Pc * (c - 1)
    h_lap = gridx.DD * h + h * gridy.DD'

    tmp = gridx.D * (h .* ubar) + (h .* vbar) * gridy.D'
    dh = osmo - tmp - Je

    dp = @. -h_lap - p
    tmp = gridx.D * (h .* c_x) + (h .* c_y) * gridy.D'
    dc = @. (invPec * tmp - osmo * c + Je * c) / h - (ubar * c_x + vbar * c_y)
    pack!(du, dh, dp, dc)
end

function twodim_solve(m, n, tspan, xw=0.5, yw=0.5, vb=0.1,d = 4.5e-6,
    solver=QNDF(linsolve=KrylovJL_GMRES()),
    tol=1e-6
    )

    # discretization
    gridx, gridy = Grid1D(m), Grid1D(n)
    M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
    u0 = pack(ones(m, n), zeros(m, n), ones(m, n))
    index = (h=1:m*n, p=m*n+1:2*m*n, c=2*m*n+1:3*m*n, size=(m,n))

    # physical parameters
    consts = PhysicalConstants()
    pde_consts = PDEConstants(consts, 1e-6 / (60vb),d)

    # evaporative flux
    hump = [ exp(-(x/xw)^2 / 2 - (y/yw)^2 / 2) for x in gridx.t, y in gridy.t ]
    center = argmax(hump)
    Je = vb .+ (1-vb) * hump

    # callback to halt integration when the solution is too thin
    condition(u,t,integrator)=(reshape(u[1:m*n],(m,n)))[center] < 1/4.5
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)

    params = (Je, index, gridx, gridy, pde_consts)
    dudt = ODEFunction(twodim_ode!, mass_matrix=M)
    prob_hpc = ODEProblem(dudt, u0, tspan, params)

    if isnothing(solver)
        sol_hpc = solve(prob_hpc, reltol=tol, abstol=tol, progress=true, progress_steps=20)
    else
        sol_hpc = solve(prob_hpc, solver, callback=cb, reltol=tol, abstol=tol, progress=true, progress_steps=20)
    end

    h = t -> reshape(sol_hpc(t, idxs=index.h), index.size)
    p = t -> reshape(sol_hpc(t, idxs=index.p), index.size)
    c = t -> reshape(sol_hpc(t, idxs=index.c), index.size)

    return gridx, gridy, sol_hpc, h, p, c
end
