using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf
include("TF1D/utility.jl")

  
function solve1d(m,xw,vb,Pc,invPec,invPecf)  
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
K1 = zeros(m,n+1)
for i = 1:n
    K1[m+1-i,n+2-i]=-1
    K1[m,n+1]=1
    K1[i,n+1-i]=1   
end
Dx_2 = -(T*Dx*K1)
Dx = T*Dx*K
Dxx = T*Dxx*K
n=n+1
function TF1d(dv,v,p,t)  
    
    h=v[1:n]
    p=v[n+1:2n]
    c=v[2n+1:3n]
    f=v[3n+1:4n]

    c_x=Dx*c
    f_x=Dx*f

  
    Jval=vb.+(1-vb)*exp.(-(x/xw).^2/2)
    
    ubar=(-h.^2/12).*(Dx*p)
    
    dv[1:n].= Pc*(c.-1)-Dx_2*(h.*ubar)-Jval
   
  
    dv[n+1:2n].=-Dxx*h - p
   

    dv[2n+1:3n].=(invPec*Dx_2*(h.*c_x)-Pc*(c.-1).*c+Jval.*c)./h-(ubar.*c_x)

    dv[3n+1:4n].=(invPecf*Dx_2*(h.*f_x)-Pc*(c.-1).*f+Jval.*f)./h-(ubar.*f_x)

       
end

M=Diagonal([ones(n);zeros(n);ones(n);ones(n)])
f = ODEFunction(TF1d,mass_matrix=M)



u0=[ones(n);Dxx*ones(n);ones(n);ones(n)]
tspan=(0.0,5.0)
condition(v,t,integrator)=v[1] < 0.111
   
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
prob_mm = ODEProblem(f,u0,tspan)

return x,solve(prob_mm,callback=cb,reltol=1e-8,abstol=1e-8)
end
x10,sol10 = solve1d(80,0.5,0.1,0.392,1/6.76,1/27.7);