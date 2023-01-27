using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings
using FFTW
using ToeplitzMatrices
L=0.540;
function fourier(m=16,n=10)
    unvec=z->reshape(z,m,n);
    hx = 2π / m
x = hx * (1:m)
entry(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hx / 2)
Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]

entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]

hy = 2π / n
y = hy * (1:n)
entry3(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hy / 2)
Dy = [ entry(mod(i-j,n)) for i in 1:n, j in 1:n ]

entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]


function TF2d(du,u,params,t)
    
    #vb=v_min/v_max  v_min=1μm/min, v_max=10μm/min 
    h=u[1:m*n]
    p=u[m*n+1:2m*n]
    c=u[2m*n+1:3m*n]
    H=unvec(h)
    P=unvec(p)
    C=unvec(c)
    ubar=(-H.^2/12).*(Dx*P)
    vbar = (-H.^2/12).*(P*Dy')
    c_x = Dx*C 
    c_y = C*Dy'
    
    osmo = params.Pc*(C .- 1)
    u_xx = Dxx*H + H*Dyy'
    
    tmp = Dx*(H.*ubar)+(H.*vbar)*Dy'
    J=(x,y,P)->constants.vb.+(1-constants.vb)*hump(x,y).+ constants.alpha*P 
    Jval = J(x,y,P)
    du[1:m*n] = vec(osmo - tmp .-Jval)
    du[m*n+1:2m*n] = vec(-u_xx - params.A./H.^3 - P)
    tmp2 = Dx*(H.*c_x)+(H.*c_y)*Dy'
    du[2m*n+1:3m*n] = vec((params.invPec*tmp2 - osmo.*C + Jval.*C)./H - (ubar.*c_x + vbar.*c_y))
    
end
    
M = Diagonal([ones(m*n);zeros(m*n);ones(m*n)])
f = ODEFunction(TF2d,mass_matrix=M)


hump(x,y)=[exp(-(x/constants.xw)^2/2)*exp(-(y/constants.xw)^2/2) for x in x, y in y]


 #J=(x,y,p)->constants.vb.+(1-constants.vb)*hump(x,y).+ constants.alpha*p 

u0 = ones(3m*n)

constants = (vb=0.1,xw=0.5/L,alpha=4.06e-2,A=5.5e-3,Pc=0.392,invPec=1/6.76)

tspan=(0.0,3.0)

prob_mm = ODEProblem(f,u0,tspan,constants)


return x,y, solve(prob_mm,reltol=1e-8,abstol=1e-8)
end

x,y,sol=fourier();

K=ones(10,10)
K=Diagonal(ones(10))
Dyy*K


m=16
n=10

hx = 2π / m
x = hx * (1:m)

x=-π.+hx*(0:m-1)
entry(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hx / 2)
Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]

entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]

hy = 2π / n
y = hy * (1:n)
entry3(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hy / 2)
Dy = [ entry(mod(i-j,n)) for i in 1:n, j in 1:n ]

entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]