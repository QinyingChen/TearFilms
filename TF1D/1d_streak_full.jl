using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf
include("utility.jl")
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
constants = (alpha=4.06e-2, A=5.5e-3,xw=0.5,Pc=0.392, invPec=1/6.76, invPecf=1/27.7)
u0=[ones(n);Dxx*ones(n)+constants.A ./ones(n);ones(n)]
tspan=(0.0,10.0)

prob_mm = ODEProblem(f,u0,tspan,constants)



return x,solve(prob_mm,reltol=1e-8,abstol=1e-8)
end
x,sol= solve1d(128);

function solve1df(n)
  

    hx = 2π / n
    x = -π .+ hx*(1:n)
    
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,n)) for i in 1:n, j in 1:n ]

function TF1df(dv,v,params,t)  
    h=sol(t)[1:n]
    p=sol(t)[n+1:2n]
    c=sol(t)[2n+1:end]
    f=v[1:n]
    f_x=Dx*f

    vb = 0.1
    J=(x,p)->vb.+(1-vb)*exp.(-(x/params.xw).^2/2).+ params.alpha*p

    Jval=J.(x,p)
    
    ubar=(-h.^2/12).*(Dx*p)

    dv[1:n].=(params.invPecf*Dx*(h.*f_x)-params.Pc*(c.-1).*f+Jval.*f)./h-(ubar.*f_x)
       
end


f = ODEFunction(TF1df)
constants = (alpha=4.06e-2, A=5.5e-3,xw=0.5,Pc=0.392, invPec=1/6.76, invPecf=1/27.7)
u0=ones(n)
tspan=(0.0,10.0)

prob_mm = ODEProblem(f,u0,tspan,constants)



return x,solve(prob_mm,reltol=1e-8,abstol=1e-8)
end
x,solf= solve1df(128);

n=128;
H=sol(3)[1:n]
P=sol(10)[n+1:2n]
C=sol(5)[2n+1:3n]
F=solf(10)[1:n]

plot!(x,C)


d = 4.5e-6
ϵf = 1.75e+7  #Napierian extinction coefficient
fcr = 0.0053   #critical fluorescein concentration
Φ = ϵf * fcr * d
F0 = ones(n);
H0= ones(n);
FI = ((-exp.(-Φ * F0 .* H0)) .+ 1) ./ ((F0 .^ 2) .+ 1)
FI[1]
I0 = 1 / FI[1]
plot(x,I0*(((-exp.(-Φ * solf(10)[1:n] .* sol(10)[1:n])) .+ 1) ./ (((solf(10)[1:n]) .^ 2) .+ 1)))

anim = @animate for t in range(0,10,length=81)
    plot(x,sol(t)[1:n],
    xaxis=("x",(-3,3)),yaxis=("y",(0,1.5)),
    title=@sprintf("thickness, t=%.3f",t),
    dpi=100)
    end
    mp4(anim,"1dh_streak.mp4")

    anim = @animate for t in range(0,10,length=81)
        plot(x,sol(t)[n+1:2n],
        xaxis=("x",(-3,3)),yaxis=("y",(-11,1)),
        title=@sprintf("pressure, t=%.3f",t),
        dpi=100)
        end
        mp4(anim,"1dp_streak.mp4")


anim = @animate for t in range(0,10,length=81)
    plot(x,sol(t)[2n+1:end],
    xaxis=("x",(-3,3)),yaxis=("y",(0,4)),
    title=@sprintf("osmolarity, t=%.3f",t),
    dpi=100)
    end
    mp4(anim,"1dc_streak.mp4")


    anim = @animate for t in range(0,10,length=81)
        plot(x,solf(t)[1:n],
        xaxis=("x",(-3,3)),yaxis=("y",(0,8)),
        title=@sprintf("f, t=%.3f",t),
        dpi=100)
        end
        mp4(anim,"1df_streak.mp4")

        anim = @animate for t in range(0,10,length=81)
            plot(x,I0*(((-exp.(-Φ * solf(t)[1:n] .* sol(t)[1:n])) .+ 1) ./ (((solf(t)[1:n]) .^ 2) .+ 1)),
            xaxis=("x",(-3,3)),yaxis=("y",(0,1.5)),
            title=@sprintf("Intensity, t=%.3f",t),
            dpi=100)
            end
            mp4(anim,"1dI_streak.mp4")

   



