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


function fourier(m, n;
    tspan=(0.0,6.0),
    solver=QNDF(linsolve=KrylovJL_GMRES()),
    tol=1e-8
    )
    hx = 2π / m
    nx = Int(m/2)
    x = (-π .+ hx*(1:m))[nx:end]
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    ny = Int(n/2)
    y = (-π .+ hy*(1:n))[ny:end]
    entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

    Kx=zeros(m,nx+1)
Tx = zeros(nx+1,m)
for i = 1:nx
    Kx[m+1-i,nx+2-i]=1
    Kx[i,nx+1-i]=1    
end
K1 = zeros(m,nx+1)
for i = 1:nx
    K1[m+1-i,nx+2-i]=-1
    K1[m,nx+1]=1
    K1[i,nx+1-i]=1   
end
for i = 1:nx+1
    Tx[nx+2-i,m+1-i]=1
end
Dx_2 = -(Tx*Dx*K1)
Dx = Tx*Dx*Kx
Dxx = Tx*Dxx*Kx

Ky=zeros(n,ny+1)
Ty = zeros(ny+1,n)
for i = 1:ny
    Ky[n+1-i,ny+2-i]=1
    Ky[i,ny+1-i]=1    
end
K2 = zeros(n,ny+1)
for i = 1:ny
    K2[n+1-i,ny+2-i]=-1
    K2[n,ny+1]=1
    K2[i,ny+1-i]=1   
end
for i = 1:ny+1
    Ty[ny+2-i,n+1-i]=1
end
Dy_2 = -(Ty*Dy*K2)
Dy = Ty*Dy*Ky
Dyy = Ty*Dyy*Ky

    function TF2d(du,u,params,t) 
        h, p, c = unpack(u, (nx+1,ny+1)) 
        #hump = [ exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y ]
        J=(p)->constants.vb.+(1-constants.vb)*hump.+ constants.alpha*p 
        Jval = J(p)     
        #Jval=constants.vb.+(1-constants.vb)*hump     
        ubar = (-h.^2/12).*(Dx*p)
        vbar = (-h.^2/12).*(p*Dy')
        c_x = Dx*c 
        c_y = c*Dy'
       
        osmo = params.Pc*(c .- 1)
        h_lap = Dxx*h + h*Dyy'
        
        tmp = Dx_2*(h.*ubar)+(h.*vbar)*Dy_2'
        dh = @. osmo - tmp - Jval
        dp = @. -h_lap - p
        tmp = Dx_2*(h.*c_x) + (h.*c_y)*Dy_2'
        dc = @. (params.invPec*tmp - osmo*c + Jval*c)/h - (ubar*c_x + vbar*c_y)
        pack!(du, dh, dp, dc)
        update!(prog, tspan[2]-t)
    end
    
    M = Diagonal([ones((nx+1)*(ny+1)); zeros((nx+1)*(ny+1)); ones((nx+1)*(ny+1))])
    dudt = ODEFunction(TF2d, mass_matrix=M)
    #hump = [ exp(-x^2-y^2).*exp(-4*(y+x^2)^2) for x in x, y in y ]
    hump = [ exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y ]
   
    
    constants = (vb=0.1,alpha=4.06e-2, Pc=0.392, invPec=1/6.76, invPecf=1/27.7)
    u0 = pack(ones(nx+1,ny+1), zeros(nx+1,ny+1), ones(nx+1,ny+1)) 
    prob_hpc = ODEProblem(dudt, u0, tspan, constants)

    
    condition(u,t,integrator)=norm(abs.(reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1))) - reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1)))>0
   #condition(u,t,integrator)=findmin(reshape(u[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))) != ((reshape(u[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1)))[1,1],CartesianIndex(1,1))
    #condition(u,t,integrator)=norm(abs.(reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1))) - reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1)))>0 || findmin(reshape(u[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))) != ((reshape(u[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1)))[1,1],CartesianIndex(1,1))
    #condition(u,t,integrator)=norm(abs.(reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1))) - reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1)))>0 || (reshape(u[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1)))[1,1] > (reshape(u[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1)))[2,1]
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)

    prog = ProgressThresh(0.0, 0.5)
    if isnothing(solver)
        sol_hpc = solve(prob_hpc, reltol=tol, abstol=tol)
    else
        sol_hpc = solve(prob_hpc, solver, callback=cb,reltol=tol, abstol=tol)
        #sol_hpc = solve(prob_hpc, solver,reltol=tol, abstol=tol)
    end
    update!(prog, 0.0)
    

    return x, y, sol_hpc
end
x4, y4, sol_hpc4 = fourier(40,40);

m = 40;
nx = Int(m/2);
n = 40;
ny = Int(n/2);
H=reshape(sol_hpc4(1.7)[1:(nx+1)*(ny+1)],(nx+1,ny+1))
P=reshape(sol_hpc4(1.7)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
C=reshape(sol_hpc4(1.7)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
surface(x4,y4,P)
