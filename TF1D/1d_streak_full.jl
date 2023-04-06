using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf

# 1d streak full solution
function solve1d(n)
  

    hx = 2π / n
    x = -π .+ hx*(1:n)
    
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,n)) for i in 1:n, j in 1:n ]

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
   
    #dv[n+1:2n].=-Dxx*h - p
 
    dv[n+1:2n].=-Dxx*h - params.A ./h.^3 - p
 
    dv[2n+1:3n].=(params.invPec*Dx*(h.*c_x)-params.Pc*(c.-1).*c+Jval.*c)./h-(ubar.*c_x)
       
end

M=Diagonal([ones(n);zeros(n);ones(n)])
f = ODEFunction(TF1d,mass_matrix=M)
constants = (alpha=4.06e-2, A=5.5e-3,xw=0.5,Pc=0.392, invPec=1/6.76)
u0=[ones(n);Dxx*ones(n)+constants.A ./ones(n);ones(n)]
tspan=(0.0,10.0)

prob_mm = ODEProblem(f,u0,tspan,constants)



return x,solve(prob_mm,reltol=1e-8,abstol=1e-8)
end
x,sol= solve1d(64);

n=64;
H=sol(10)[1:n]
P=sol(10)[n+1:2n]
C=sol(10)[2n+1:end]


plot(x,P)

