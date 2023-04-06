using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf
include("utility.jl")
# 1d spot with symmetry
m=64;  # change m, and run the whole thing
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
  
function solve1d(n)

function TF1d(dv,v,params,t)  
    h=v[1:n]
    p=v[n+1:2n]
    c=v[2n+1:3n]

    c_x=Dx*c

    vb = 0.1
    J=(x,p)->vb.+(1-vb)*exp.(-(x/params.xw).^2/2).+ params.alpha*p

    Jval=J.(x,p)
    
    ubar=((-h.^2/12).*(Dx*p))
    # Deal with singularity at 0
    h_xx = [(Dxx*h)[1];(Dx_2*(x.*(Dx*h)))[2:end]./x[2:end]]
    tmp = [(Dx_2*(h.*ubar))[1];(Dx_2*(x.*(h.*ubar)))[2:end]./x[2:end]]
    tmp2 = [(Dx_2*(h.*c_x))[1];(Dx_2*(x.*(h.*c_x)))[2:end]./x[2:end]]

    dv[1:n].= params.Pc*(c.-1)-tmp.-Jval
   
    dv[n+1:2n].=-h_xx - params.A ./h.^3 - p
 
    dv[2n+1:3n].=(params.invPec*tmp2-params.Pc*(c.-1).*c+Jval.*c)./h-(ubar.*c_x)
       
end

M=Diagonal([ones(n);zeros(n);ones(n);ones(n)])
f = ODEFunction(TF1d,mass_matrix=M)
constants = (alpha=4.06e-2, A=5.5e-3, xw=0.5,Pc=0.392, invPec=1/6.76)
u0=[ones(n);[(Dxx*ones(n))[1];(Dx_2*(x.*(Dx*ones(n))))[2:end]./x[2:end]]+constants.A ./ones(n);ones(n);ones(n)]
tspan=(0.0,10.0)

prob_mm = ODEProblem(f,u0,tspan,constants)

#condition(v,t,integrator)=norm(abs.(v[1:n]) - v[1:n])>0  
#condition(v,t,integrator)=v[n+1]>v[n+2]

  # affect!(integrator) = terminate!(integrator)
   # cb = DiscreteCallback(condition,affect!)

#return x,solve(prob_mm,callback=cb,reltol=1e-8,abstol=1e-8)
return x,solve(prob_mm,reltol=1e-8,abstol=1e-8)
end
x,sol = solve1d(n+1);

H=sol(10)[1:n+1]
P=sol(10)[n+2:2(n+1)]
C=sol(10)[2(n+1)+1:end]
-Dxx*H

plot(x,P)
plot(x,sol(2)[1:n+1])  # h
plot(x,sol(2)[n+2:2(n+1)])  # p
plot(x,sol(2)[2(n+1)+1:3(n+1)]) 

X = [-(flip(x))[2:end-1];x]
Hfull = [(flip(H))[2:end-1];H]
plot(X,Hfull)

Pfull = [(flip(P))[2:end-1];P]
plot(X,Pfull)


