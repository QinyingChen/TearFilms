using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf

using ComponentArrays 

function newcheb(N)
    x = [ 1+cos(pi*j/N) for j  in 0:N ]
    c(n) = (n==0) || (n==N) ? 2 : 1
    entry(i,j) = i==j ? 0 : c(i)/c(j) * (-1)^(i+j) / (x[i+1] - x[j+1])
    D = [ entry(i,j) for i in 0:N, j in 0:N ]
    D  = D - diagm(vec(sum(D,dims=2)));    # diagonal entries
    return D, x
end

function extend(V)
    Uy = endvals*V
    return [Uy[1]; V; Uy[2]]
end
N=40;
D,x=newcheb(N)
ends = [1, N+1]
endvals = -D[ends, ends] \ D[ends,2:N]
D2=D^2
Dxx=D2
Dx=D


function solve1d(n)
function TF1d(du,u,params,t)
   
    H=extend(u.h)
    P=extend(u.p)
    C=extend(u.c)
    F=extend(u.f)       
    Jval = J + params.alpha*P
    ubar = (-H.^2/12).*(Dx*P)
    c_x = Dx*C
    f_x = Dx*F
    osmo = params.Pc*(C .- 1)
    u_xx = Dxx*H
   
    tmp = Dx*(H.*ubar)
    @. du.h = (osmo - tmp - Jval)[2:end-1]
   
    @. du.p = (-u_xx*H - params.A*H^(-3)- P)[2:end-1]
       
    tmp = Dx*(H.*c_x)
    @. du.c = ((params.invPec*tmp - osmo*C + Jval*C)/H - (ubar*c_x))[2:end-1]
 
    tmp = Dx*(H.*f_x)
    @. du.f = ((params.invPec*tmp - osmo*F + Jval*F)/H - (ubar*f_x))[2:end-1]
    
end
    
M = Diagonal([ones(n);zeros(n);ones(n);ones(n)])
f = ODEFunction(TF1d,mass_matrix=M)
hump=[exp(-(x/0.5)^2/2) for x in x]
vb=0.1
J = @. vb + (1-vb)*hump
constants = (J=J,vb=0.1,xw=0.5,alpha=4.06e-2,A=5.5e-3,Pc=0.392,invPec=1/6.76)


#J = (x,p) -> constants.vb + (1-constants.vb)*hump(x) + constants.alpha*p


 
u0 = ComponentArray(h=ones(n),p=(Dxx*extend(ones(n))+constants.A*extend(ones(n)))[2:end-1],c=ones(n),f=ones(n))

tspan=(0.0,5.0)
prob_mm = ODEProblem(f,u0,tspan,constants)

return x,solve(prob_mm,reltol=1e-8,abstol=1e-8)
end
x,sol=solve1d(39);
Q=sol(4).p
plot(x,extend(Q))

t=range(0,3,21);
u = sol.(t);
h = hcat([extend(u.h) for u in u]...);
p = hcat([extend(u.p) for u in u]...);
c = hcat([extend(u.c) for u in u]...);
f = hcat([extend(u.f) for u in u]...);

function flip(x)
    z = zeros(length(x))   
       for i=1:length(x)
           z[i]=x[length(x)+1-i]
       end
   
   return z
   end
   function flipm(x)
   z = zeros(size(x))
   for j=1:size(x)[2]
   
       z[:,j]=flip(x[:,j])
   end
   return z
   end
   H=vcat(h,flipm(h))
   P=vcat(p,flipm(p))
   C=vcat(c,flipm(c))
   F=vcat(f,flipm(f))
   X=[-x;flip(x)]

label = hcat([@sprintf("t=%.2f",t) for t in t]...)
opts = (line_z=t',label=label,color=:viridis,legend=false,colorbar=true)
plot(layout=(4,1),link=:x)
plot!(X,H,subplot=1,xlabel="x",ylabel="thickness";opts...)
plot!(X,P,subplot=2,xlabel="x",ylabel="pressure";opts...)
plot!(X,C,subplot=3,xlabel="x",ylabel="osmolarity";opts...)
plot!(X,F,subplot=4,xlabel="x",ylabel="fluorescein";opts...)






        


