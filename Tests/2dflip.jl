using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using Krylov
using LinearSolve
using Plots
using LaTeXStrings,Printf
function newcheb(N)
    x = [ 1+cos(pi*j/N) for j  in 0:N ]
    c(n) = (n==0) || (n==N) ? 2 : 1
    entry(i,j) = i==j ? 0 : c(i)/c(j) * (-1)^(i+j) / (x[i+1] - x[j+1])
    D = [ entry(i,j) for i in 0:N, j in 0:N ]
    D  = D - diagm(vec(sum(D,dims=2)));    # diagonal entries
    return D, x
end

function extendy(V)
    Uy = V*endvals_y'
    return [Uy[:,1] V Uy[:,2] ]
end
function extendx(V)
    Ux = endvals_x*V
    return [Ux[1,:]' ;V;Ux[2,:]' ]
end
m=16;
n=16;
Dx,x=newcheb(m)
ends_x = [1,m+1]
endvals_x = -Dx[ends_x, ends_x] \ Dx[ends_x,2:m]
Dy,y=newcheb(n)
ends_y = [1,n+1]
endvals_y = -Dy[ends_y, ends_y] \ Dy[ends_y,2:n]
Dxx=Dx^2
Dyy=Dy^2

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

function chop(x)
return x[2:end-1,2:end-1]
end

unpack(u, sz) = unpack!( similar(u, sz), similar(u, sz), similar(u, sz), u )

function cheby(m, n;
    tspan=(0.0,3.0),
    solver=QNDF(linsolve=KrylovJL_GMRES()),
    tol=1e-8
    )
  
    function TF2d(du,u,params,t) 
        h, p, c = unpack(u, (m,n)) 
        H=extendy(extendx(h))  
        P=extendy(extendx(p))
        C=extendy(extendx(c)) 
        Jval = J + params.alpha*P       
        ubar = (-H.^2/12).*(Dx*P)
        vbar = (-H.^2/12).*(P*Dy')
        c_x = Dx*C 
      
        c_y = C*Dy'
       
       
        osmo = params.Pc*(C .- 1)
        h_lap = Dxx*H + H*Dyy'
        
        tmp = Dx*(H.*ubar)+(H.*vbar)*Dy'
         dh =  chop(@. osmo - tmp - Jval)
        dp =  chop(@. -h_lap - params.A*H^(-3) - P)
        tmp = Dx*(H.*c_x) + (H.*c_y)*Dy'
        dc =  chop(@. (params.invPec*tmp - osmo*C + Jval*C)/H - (ubar*c_x + vbar*c_y))
        pack!(du, dh, dp, dc)
        update!(prog, tspan[2]-t)
    end
    
    M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
    dudt = ODEFunction(TF2d, mass_matrix=M)
    #hump = [ exp(-x^2-y^2).*exp(-4*(y+x^2)^2) for x in x, y in y ]
    
    hump = [ exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y ]
    vb = 0.1
    J = @. vb + (1-vb)*hump 

    constants = (J=J, A=5.5e-3, alpha=4.06e-2,Pc=0.392, invPec=1/6.76, invPecf=1/27.7)
    u0 = pack(ones(m,n), constants.A*ones(m,n), ones(m,n)) 
    prob_hpc = ODEProblem(dudt, u0, tspan, constants)
    
    prog = ProgressThresh(0.0, 0.5)
    if isnothing(solver)
        sol_hpc = solve(prob_hpc, reltol=tol, abstol=tol)
    else
        sol_hpc = solve(prob_hpc, solver, reltol=tol, abstol=tol)
    end
    update!(prog, 0.0) 

    return x, y, sol_hpc
end
x, y, sol_hpc = cheby(m-1,n-1);

X=[-x;flip(x)]
Y=[-y;flip(y)]

R=reshape(sol_hpc(2)[1:m*n],(m,n))
h1=extendy(extendx(reshape(sol_hpc(2)[1:(m-1)*(n-1)],(m-1,n-1))))
p1=extendy(extendx(reshape(sol_hpc(2)[(m-1)*(n-1)+1:2(m-1)*(n-1)],(m-1,n-1))))

H=vcat(h1,flipc(h1))
TR=vcat(flipr(h1),flipc(flipr(h1)))
surface(X,Y,[P TR])
surface(x,y,p1)
P=vcat(p1,flipc(p1))
TR=vcat(flipr(p1),flipc(flipr(p1)))



# H2=hcat(H,flipr(H))
# PLO=[H flipr(H)]
# HH=hcat(flipm(H),H)
# P=vcat(p1,flipc(p1))
# PLOO=[P flipr(P)]
# C=vcat(c,flipm(c))
# F=vcat(f,flipm(f))

















