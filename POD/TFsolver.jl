using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using Krylov
using LinearSolve
using Plots
using LaTeXStrings,Printf
using JLD2
function pack!(u, h, p, c) 
    m, n = size(h)
    mn = m*n
    u[1:mn] .= vec(h)
    u[mn+1:2mn] .= vec(p)
    u[2mn+1:3mn] .= vec(c)
    return u 
end

pack(h,p,c) = pack!(similar(h,(3*length(h),)), h, p, c)

function unpack!(h, p, c, u)
    sz = size(h)
    mn = length(h)
    h .= reshape( u[1:mn], sz )
    p .= reshape( u[mn+1:2mn], sz )
    c .= reshape( u[2mn+1:3mn], sz )
    return h, p, c 
end

unpack(u, sz) = unpack!( similar(u, sz), similar(u, sz), similar(u, sz), u )


function fourier(m, n,
    xw,
    yw,
    vb,
    Pc,
    invPec;
    tspan=(0.0,5.0),
    solver=QNDF(linsolve=KrylovJL_GMRES()),
    tol=1e-8
    )
    hx = 2π / m
    x = (-π .+ hx*(1:m))
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    y = (-π .+ hy*(1:n))
    entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

    

    function TF2d(du,u,params,t) 
        h, p, c = unpack(u, (m,n)) 
       
        ubar = (-h.^2/12).*(Dx*p)
        vbar = (-h.^2/12).*(p*Dy')
        c_x = Dx*c 
        c_y = c*Dy'
       
        osmo = Pc*(c .- 1)
        h_lap = Dxx*h + h*Dyy'
        
        tmp = Dx*(h.*ubar)+(h.*vbar)*Dy'
        dh = @. osmo - tmp - Jval
        
        dp = @. -h_lap  - p
        tmp = Dx*(h.*c_x) + (h.*c_y)*Dy'
        dc = @. (invPec*tmp - osmo*c + Jval*c)/h - (ubar*c_x + vbar*c_y)
        pack!(du, dh, dp, dc)
        update!(prog, tspan[2]-t)
    end
    
    M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
    dudt = ODEFunction(TF2d, mass_matrix=M)
    
    hump = [ exp(-(x/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ]
    Q = findmax(hump)[2]
    Jval = vb.+(1-vb)*hump
    
    u0 = pack(ones(m,n), Dxx*ones(m,n)+ones(m,n)*Dyy', ones(m,n))
    prob_hpc = ODEProblem(dudt, u0, tspan)
    
    
    condition(u,t,integrator)=(reshape(u[1:m*n],(m,n)))[Q[1],Q[2]] < 1/4.5
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
    prog = ProgressThresh(0.0, 0.5)
    if isnothing(solver)
        sol_hpc = solve(prob_hpc, reltol=tol, abstol=tol)
    else
        sol_hpc = solve(prob_hpc, solver, callback=cb, reltol=tol, abstol=tol)
    end
    update!(prog, 0.0)
       
    return x, y, sol_hpc
end



#m=40;n=40;xw=0.5;yw=0.5;vb=0.1;Pc=0.392;invPec=1/6.76;

 #x, y, sol_hpc = fourier(m,n,xw,yw,vb,Pc,invPec);

 