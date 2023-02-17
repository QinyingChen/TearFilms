function pack!(h, p, c) 
    
    return [vec(h);vec(p);vec(c)]
end

function unpack!(u)
    sz = (m,n)
    mn = m*n
   
    return reshape( u[1:mn], sz ),reshape( u[mn+1:2mn], sz ),reshape( u[2mn+1:3mn], sz )
     
end

m=42;
n=42;
function fourier(m, n;
    tspan=(0.0,3.0),
    solver=QNDF(linsolve=KrylovJL_GMRES()),
    tol=1e-8
    )
    
    hx = 2π / m
    x = -π .+ hx*(1:m)
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    y = -π .+ hy*(1:n)
    entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

    function TF2d(du,u,params,t) 
        h, p, c = unpack!(u)           
        ubar = (-h.^2/12).*(Dx*p)
        vbar = (-h.^2/12).*(p*Dy')
        c_x = Dx*c 
        c_y = c*Dy'
       
        osmo = params.Pc*(c .- 1)
        h_lap = Dxx*h + h*Dyy'
        
        tmp = Dx*(h.*ubar)+(h.*vbar)*Dy'
        dh = @. osmo - tmp - J
        dp = @. -h_lap - params.A/h^3 - p
        tmp = Dx*(h.*c_x) + (h.*c_y)*Dy'
        dc = @. (params.invPec*tmp - osmo*c + J*c)/h - (ubar*c_x + vbar*c_y)
        du .= pack!(dh, dp, dc)
        update!(prog, tspan[2]-t)
    end
    
    M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
    dudt = ODEFunction(TF2d, mass_matrix=M)
    
    hump = [ exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y ]
    vb = 0.1
    J = @. vb + (1-vb)*hump 

    constants = (J=J, ⍺=4.06e-2, A=0, Pc=0.392, invPec=1/6.76, invPecf=1/27.7)
    u0 = pack!(ones(m,n), constants.A*ones(m,n), ones(m,n)) 
    prob_hpc = ODEProblem(dudt, u0, tspan, constants)
    
    prog = ProgressThresh(0.0, 0.5)
    if isnothing(solver)
        sol_hpc = solve(prob_hpc, reltol=tol, abstol=tol)
    else
        sol_hpc = solve(prob_hpc, solver, reltol=tol, abstol=tol)
    end
    update!(prog, 0.0)
   
    function TF2df(df, f, params, t)
        f_x = Dx*f 
        f_y = f*Dy'
        h, p, c = unpack!(sol_hpc(t))  
        ubar = (-h.^2/12) .* (Dx*p)
        vbar = (-h.^2/12) .* (p*Dy')
        osmo = params.Pc*(c .- 1)
        tmp = Dx*(h.*f_x) + (h.*f_y)*Dy'
      
        @. df = (params.invPecf*tmp - osmo*f + J*f)/h - (ubar*f_x + vbar*f_y)
    end
    
    constants = (J=J, ⍺=4.06e-2, A=0, Pc=0.392, invPec=1/6.76, invPecf=1/27.7)

    dfdt = ODEFunction(TF2df)
    f0 = ones(m,n)
    prob_f = ODEProblem(dfdt, f0, tspan, constants)
    sol_f = solve(prob_f, reltol=tol, abstol=tol)

    return x, y, sol_hpc, sol_f
end
@elapsed x, y, sol_hpc, sol_f=fourier(42,42)
