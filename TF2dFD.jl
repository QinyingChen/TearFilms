using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf
using LinearSolve
using ComponentArrays  # use field names instead of index ranges
using Krylov

L=0.54   

function solve_1d(a=-3,b=3,c=-3,d=3,m=15,n=10)
   

    hx = (b-a)/m
    hy=(d-c)/n
    x = @. a + hx*(0:m-1) 
    y=@. c + hy*(0:n-1)

    dp = fill(0.5/hx,m-1) 
    Dx = diagm(-1=>-dp,1=>dp)
    Dx[1,m] = -1/(2*hx)
    Dx[m,1] = 1/(2*hx)

    d0 = fill(-2/hx^2,m) 
    dp = ones(m-1)/hx^2 
    Dxx = diagm(-1=>dp,0=>d0,1=>dp)
    Dxx[1,m] = 1/(hx^2)
    Dxx[m,1] = 1/(hx^2)

    dp = fill(0.5/hy,n-1) 
    Dy = diagm(-1=>-dp,1=>dp)
    Dy[1,n] = -1/(2*hy)
    Dy[n,1] = 1/(2*hy)

    d0 = fill(-2/hy^2,n) 
    dp = ones(n-1)/hy^2 
    Dyy = diagm(-1=>dp,0=>d0,1=>dp)
    Dyy[1,n] = 1/(hy^2)
    Dyy[n,1] = 1/(hy^2)
    
    
    


    function TF1d(du,u,params,t)
    
        #vb=v_min/v_max  v_min=1μm/min, v_max=10μm/min 
        Jval = J(x,y,u.p)
            
        ubar=(-u.h.^2/12).*(Dx*u.p)
        vbar = (-u.h.^2/12).*(u.p*Dy')
        c_x = Dx*u.c 
        c_y = u.c*Dy'
        
        osmo = params.Pc*(u.c .- 1)
        u_xx = Dxx*u.h + u.h*Dyy'
        
        tmp = Dx*(u.h.*ubar)+(u.h.*vbar)*Dy'
        @. du.h = osmo - tmp - Jval
        @. du.p = -u_xx - params.A/u.h^3 - u.p
        tmp = Dx*(u.h.*c_x)+(u.h.*c_y)*Dy'
        @. du.c = (params.invPec*tmp - osmo*u.c + Jval*u.c)/u.h - (ubar*c_x + vbar*c_y)
    end
        
    M = Diagonal([ones(m*n);zeros(m*n);ones(m*n)])
    f = ODEFunction(TF1d,mass_matrix=M)

    
    hump(x,y)=[exp(-(x/constants.xw)^2/2)*exp(-(y/constants.xw)^2/2) for x in x, y in y]
   
   
     J=(x,y,p)->constants.vb.+(1-constants.vb)*hump(x,y).+ constants.alpha*p 
    
    u0 = ComponentArray(h=ones(m,n),p=ones(m,n),c=ones(m,n)) 
    constants = (vb=0.1,xw=0.5/L,alpha=4.06e-2,A=5.5e-3,Pc=0.392,invPec=1/6.76)
    
    tspan=(0.0,3.0)
    
    prob_mm = ODEProblem(f,u0,tspan,constants)
    
    return x,y, solve(prob_mm,reltol=1e-8,abstol=1e-8)
end


##

x,y,sol=solve_1d();

## solve for f

function solve_2df(a=-3,b=3,c=-3,d=3,m=15,n=10)
    hx = (b-a)/m
    hy=(d-c)/n
    x = @. a + hx*(0:m-1) 
    y=@. c + hy*(0:n-1)

    dp = fill(0.5/hx,m-1) 
    Dx = diagm(-1=>-dp,1=>dp)
    Dx[1,m] = -1/(2*hx)
    Dx[m,1] = 1/(2*hx)

    d0 = fill(-2/hx^2,m) 
    dp = ones(m-1)/hx^2 
    Dxx = diagm(-1=>dp,0=>d0,1=>dp)
    Dxx[1,m] = 1/(hx^2)
    Dxx[m,1] = 1/(hx^2)

    dp = fill(0.5/hy,n-1) 
    Dy = diagm(-1=>-dp,1=>dp)
    Dy[1,n] = -1/(2*hy)
    Dy[n,1] = 1/(2*hy)

    d0 = fill(-2/hy^2,n) 
    dp = ones(n-1)/hy^2 
    Dyy = diagm(-1=>dp,0=>d0,1=>dp)
    Dyy[1,n] = 1/(hy^2)
    Dyy[n,1] = 1/(hy^2)

function TF2df(du,u,params,t)
    f_x = Dx*u.f 
    f_y = u.f*Dy'
    w=sol(t)    
    Jval = J(x,y)    
    ubar=(-w.h.^2/12).*(Dx*w.p)
    vbar = (-w.h.^2/12).*(w.p*Dy')
    osmo = params.Pc*(w.c .- 1)
    tmp = Dx*(w.h.*f_x)+(w.h.*f_y)*Dy'
  
    @. du.f = (params.invPecf*tmp - osmo*u.f + Jval*u.f)/w.h - (ubar*f_x + vbar*f_y)
    end
M = Diagonal(ones(m*n))
f = ODEFunction(TF2df,mass_matrix=M)


hump(x,y)=[exp(-(x/constants.xw)^2/2)*exp(-(y/constants.xw)^2/2) for x in x, y in y]


 J=(x,y)->constants.vb.+(1-constants.vb)*hump(x,y)

u0 = ComponentArray(f=ones(m,n)) 
constants = (vb=0.1,xw=0.5/L,alpha=4.06e-2,A=5.5e-3,Pc=0.392,invPecf=1/27.7)

tspan=(0.0,3.0)

prob_mm = ODEProblem(f,u0,tspan,constants)

return x,y, solve(prob_mm,Rodas5(),reltol=1e-8,abstol=1e-8)
end
x,y,solf=solve_2df()







##
anim = @animate for t in range(0,3,length=81)
    surface(x,y,solf(t).f',color=:redsblues,
    xaxis=(L"x",(-3,3)),yaxis=(L"y",(-3,3)),zlims=(0,8),clims=(0,8),
    title=@sprintf("TF2Df, t=%.3f",t),
    colorbar=:none)
    end
    mp4(anim,"TF2df.mp4")
    

##
surface(x,y,sol(2).h',
xlabel=L"x",ylabel=L"y",zaxis=((0,1),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2Dh",dpi=100 )

surface(x,y,solf(2).f',
xlabel=L"x",ylabel=L"y",zaxis=((0,5),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2Df",dpi=100 )





t=range(0,3,21);
u = sol.(t)
h = hcat([u.h for u in u]...);
p = hcat([u.p for u in u]...);
c = hcat([u.c for u in u]...);
f = hcat([u.f for u in u]...);

##
##  # use Alt-enter to evaluate between the ## lines

label = hcat([@sprintf("t=%.2f",t) for t in t]...)
opts = (line_z=t',label=label,color=:viridis,legend=false,colorbar=true)
plot(layout=(4,1),link=:x)
plot!(x,h,subplot=1,title="FD vmax=10 vmin=1 xw=0.5mm",xlabel="x",ylabel="thickness";opts...)
plot!(x,p,subplot=2,xlabel="x",ylabel="pressure";opts...)
plot!(x,c,subplot=3,xlabel="x",ylabel="osmolarity";opts...)
plot!(x,f,subplot=4,xlabel="x",ylabel="fluorescein";opts...)

##

## Animations
name = Dict(:h=>"thicknessfd",:p=>"pressurefd",:c=>"osmolarityfd",:f=>"fluoresceinfd")
ylims = Dict(:h=>extrema(h),:p=>extrema(p),:c=>extrema(c),:f=>extrema(f))
which = :p
anim = @animate for t in range(0,3,length=81)
    plot(x,sol(t)[which],
        xaxis=("x",(-3.5,3.5)),yaxis=("y",ylims[which]),
        title=@sprintf("%s, t=%.2f",name[which],t)
    )
end
mp4(anim,"$(name[which]).mp4")


sol(2).h


d=4.5e-6
ϵf=1.75e+7  #Napierian extinction coefficient
fcr=0.0053   #critical fluorescein concentration
Φ=ϵf*fcr*d
FI=((-exp.(-Φ*c.*h)).+1)./((solf(2).f.^2).+1)
I0=1/FI[1]
FI=I0*((-exp.(-Φ*c.*h)).+1)./((f.^2).+1)
plot(x,FI,xlabel="x",ylabel="FI",label = false)

FI=((-exp.(-Φ*sol(0).c'.*sol(0).h')).+1)./((solf(0).f'.^2).+1)
FI[1]
I0=1/FI[1]



surface(x,y,I0*((-exp.(-Φ*sol(2).c'.*sol(2).h')).+1)./((solf(2).f'.^2).+1),
xlabel=L"x",ylabel=L"y",zaxis=((0,1),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2DI",dpi=100 )


##
anim = @animate for t in range(0,3,length=81)
    surface(x,y,I0*((-exp.(-Φ*sol(t).c'.*sol(t).h')).+1)./((solf(t).f'.^2).+1),color=:redsblues,
    xaxis=(L"x",(-3,3)),yaxis=(L"y",(-3,3)),zlims=(0,1),clims=(0,1),
    title=@sprintf("TF2DI, t=%.3f",t),
    colorbar=:none)
    end
    mp4(anim,"TF2dI.mp4")
##