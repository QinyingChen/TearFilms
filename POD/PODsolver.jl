include("TFsolver.jl")
m=40;n=40;xw=0.5;yw=0.5;vb=0.1;Pc=0.392;invPec=1/6.76;
x, y, sol_hpc = fourier(m,n,xw,yw,vb,Pc,invPec);

function PODsol(m,n,xw,yw,vb,Pc,invPec,sh,sp,sc; t=range(0,0.5,41),tspan=(0,3))
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
    
    
    
    W = zeros(m*n, length(t))
    for j in eachindex(t)
    W[:,j] = sol_hpc(t[j])[1:m*n]
    end
    U,σ,V = svd(W)
    σ
    Bh=U[:,1:sh]
    
    
    W = zeros(m*n, length(t))
    for j in eachindex(t)
    W[:,j] = sol_hpc(t[j])[m*n+1:m*n*2]
    end
    U,σ,V = svd(W)
    σ
    Bp=U[:,1:sp]
    
    W = zeros(m*n, length(t))
    for j in eachindex(t)
    W[:,j] = sol_hpc(t[j])[2*m*n+1:3*m*n]
    end
    U,σ,V = svd(W)
    σ
    Bc=U[:,1:sc]
    
    function PODfun(v, params, t)
        h,p,c = unpack([Bh*v[1:sh];Bp*v[sh+1:sh+sp];Bc*v[sh+sp+1:sh+sp+sc]], (m,n)) 
        hump = [ exp(-(x/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ]

        Jval = vb.+(1-vb)*hump
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
       
        return [Bh'*vec(dh);Bp'*vec(dp);Bc'*vec(dc)]
        
    end
    hump = [ exp(-(x/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ]
    Q = findmax(hump)[2]
    condition(v,t,integrator)=(reshape(Bh*v[1:sh],(m,n)))[Q[1],Q[2]] < 1/4.5
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
    u0 = pack(ones(m,n), Dxx*ones(m,n)+ones(m,n)*Dyy', ones(m,n))
    U1 = Bh'*u0[1:m*n]
    U2 = Bp'*u0[m*n+1:2*m*n]
    U3 = Bc'*u0[2*m*n+1:end]
    newu = [U1;U2;U3] 
    M2 = cat(Bh',Bp',Bc'; dims=(1,2))
    M3 = cat(Bh,Bp,Bc; dims=(1,2))
    
  
    M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
    POD = ODEFunction(PODfun, mass_matrix=M2*M*M3)
    
    reducedprob = ODEProblem(POD, newu, tspan)
    pod_sol = solve(reducedprob, callback=cb, reltol=1e-8, abstol=1e-8)
    
    return x, y, pod_sol, Bh, Bp, Bc
    end
  
    
    m=40;n=40;xw=0.5;yw=0.5;vb=0.1;Pc=0.392;invPec=1/6.76;sh=20;sp=30;sc=20;
    x, y, pod_sol, Bh, Bp, Bc = PODsol(m,n,xw,yw,vb,Pc,invPec,sh,sp,sc);