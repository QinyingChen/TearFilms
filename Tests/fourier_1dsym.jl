using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf
m=64;
n=Int(m/2);
hx = 2π / m
    x = -π .+ hx*(1:m)
    x = x[n:end]  # [0,pi]
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]

# Modify Dx to T*Dx*K
K=zeros(m,n+1)
T = zeros(n+1,m)
for i = 1:Int(m/2)
    K[m+1-i,n+2-i]=1
    K[i,n+1-i]=1    
end
for i = 1:n+1
    T[n+2-i,m+1-i]=1
end
Dx = T*Dx*K
Dxx = T*Dxx*K

# check derivatives
Dx_1 = T*Dx*K
Dxx_1 = T*Dxx*K

U = rand(m,1)
for i = 1:n-1
    U[i]=U[m-i]
end

V = U[n:end]

norm((Dx*U)[n:end]-Dx_1*V,Inf)         

norm((Dxx*U)[n:end]-Dxx_1*V,Inf)

# checked: Dx,Dxx are ok.

    
function solve1d(n)

function TF1d(dv,v,params,t)  
    h=v[1:n]
    p=v[n+1:2n]
    c=v[2n+1:3n]

    c_x=Dx*c

    vb = 0.1
    J=(x,p)->vb.+(1-vb)*exp.(-(x/params.xw).^2/2).+ params.alpha*p

    Jval=J.(x,p)
    
    ubar=(-h.^2/12).*(Dx*p)
    
    dv[1:n].= params.Pc*(c.-1)-Dx*(h.*ubar).-Jval
   
  
    dv[n+1:2n].=-Dxx*h - p
 # @show -Dxx*h - p   close to zero. It is ok.
    dv[2n+1:3n].=(params.invPec*Dx*(h.*c_x)-params.Pc*(c.-1).*c+Jval.*c)./h-(ubar.*c_x)
       
end

M=Diagonal([ones(n);zeros(n);ones(n)])
f = ODEFunction(TF1d,mass_matrix=M)

u0=[ones(n);Dxx*ones(n);ones(n)]
tspan=(0.0,3.0)
constants = (alpha=4.06e-2, xw=0.5,A=5.5e-3, Pc=0.392, invPec=1/6.76)
prob_mm = ODEProblem(f,u0,tspan,constants)

return x,solve(prob_mm,reltol=1e-8,abstol=1e-8)
end
x,sol = solve1d(n+1);

H=sol(2)[1:n+1]
-Dxx*H
P=sol(2)[n+2:2(n+1)]

plot(x,sol(2)[1:n+1])  # h
plot(x,sol(2)[n+2:2(n+1)])  # p
plot(x,sol(2)[2(n+1)+1:3(n+1)])  # c