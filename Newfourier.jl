using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf
using ComponentArrays

using FFTW
using ToeplitzMatrices 

L=0.540;

function fourier(m,n)
    hx = 2π / m
    x=-π.+hx*(0:m-1)
    entry(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    y=-π.+hy*(0:n-1)
    entry3(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]
    
    
    function TF2d(du,u,params,t)
    
        #vb=v_min/v_max  v_min=1μm/min, v_max=10μm/min 
        J=(x,y,p)->constants.vb.+(1-constants.vb)*hump(x,y).+ constants.alpha*p 
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
    f = ODEFunction(TF2d,mass_matrix=M)

    
    hump(x,y)=[exp(-(x/constants.xw)^2/2)*exp(-(y/constants.xw)^2/2) for x in x, y in y]
   
   
    #J=(x,y,p)->constants.vb.+(1-constants.vb)*hump(x,y).+ constants.alpha*p 
    
    u0 = ComponentArray(h=ones(m,n),p=ones(m,n),c=ones(m,n)) 
    constants = (vb=0.1,xw=0.5/L,alpha=4.06e-2,A=5.5e-3,Pc=0.392,invPec=1/6.76)
    
    tspan=(0.0,3.0)
    
    prob_mm = ODEProblem(f,u0,tspan,constants)
    
    return x,y, solve(prob_mm,QNDF(linsolve=LinSolveGMRES()),reltol=1e-8,abstol=1e-8)
end

x,y,sol=fourier(42,40);

## solve for f

function solve_2df(m,n)
    hx=2*π/m
    x=-π.+hx*(0:m-1)
    column=vcat(0,[0.5*(-1)^j.*cot(j*hx/2) for j in 1:m-1])
    Dx=Circulant(column)
    column2=vcat(-π^2/3hx^2-1/6,[-(-1)^j/(2*(sin(j*hx/2))^2) for j in 1:m-1])
    Dxx=Circulant(column2)

    hy=2*π/n
    y=-π.+hy*(0:n-1)
    column3=vcat(0,[0.5*(-1)^j.*cot(j*hy/2) for j in 1:n-1])
    Dy=Circulant(column3)
    column4=vcat(-π^2/3hy^2-1/6,[-(-1)^j/(2*(sin(j*hy/2))^2) for j in 1:n-1])
    Dyy=Circulant(column4)

function TF2df(du,u,params,t)
    f_x = Dx*u.f 
    f_y = u.f*Dy'
    w=sol(t)   
    J=(x,y)->constants.vb.+(1-constants.vb)*hump(x,y).+ constants.alpha.*w.p
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

#J=(x,y)->constants.vb.+(1-constants.vb)*hump(x,y).+ constants.alpha.*w.p 
 

u0 = ComponentArray(f=ones(m,n)) 
constants = (vb=0.1,xw=0.5/L,alpha=4.06e-2,A=5.5e-3,Pc=0.392,invPecf=1/27.7)

tspan=(0.0,3.0)

prob_mm = ODEProblem(f,u0,tspan,constants)

return x,y, solve(prob_mm,reltol=1e-8,abstol=1e-8)
end
x,y,solf=solve_2df(16,10);



t=range(0,3,21);


## try some surface plots for h,p,c


surface(x,y,sol(2).h',
xlabel=L"x",ylabel=L"y",zaxis=((0,1),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2Dh",dpi=100 )

surface(x,y,sol(2).p',
xlabel=L"x",ylabel=L"y",zaxis=((-2,1),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2Dp",dpi=100 )

surface(x,y,sol(2).c',
xlabel=L"x",ylabel=L"y",zaxis=((0,3),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2Dc",dpi=100 )

surface(x,y,solf(2).f',
xlabel=L"x",ylabel=L"y",zaxis=((0,5),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2Df",dpi=100 )


## plot FI


d=4.5e-6
ϵf=1.75e+7  #Napierian extinction coefficient
fcr=0.0053   #critical fluorescein concentration
Φ=ϵf*fcr*d
FI=((-exp.(-Φ*sol(0).f'.*sol(0).h')).+1)./((sol(0).f'.^2).+1)
FI[1]
I0=1/FI[1]

surface(x,y,I0*((-exp.(-Φ*sol(2).f'.*sol(2).h')).+1)./((sol(2).f'.^2).+1),
xlabel=L"x",ylabel=L"y",zaxis=((0,1),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-2,2),colorbar=:none,
title="TF2DI",dpi=100 )


## Animations


anim = @animate for t in range(0,3,length=81)
    surface(x,y,sol(t).f',color=:redsblues,
    xaxis=(L"x",(-3,3)),yaxis=(L"y",(-3,3)),zlims=(0,8),clims=(0,8),
    title=@sprintf("TF2Df, t=%.3f",t),
    colorbar=:none)
    end
    mp4(anim,"2df.mp4")

anim = @animate for t in range(0,3,length=81)
        surface(x,y,sol(t).h',color=:redsblues,
        xaxis=(L"x",(-3,3)),yaxis=(L"y",(-3,3)),zlims=(0,1),clims=(0,1),
        title=@sprintf("TF2Dh, t=%.3f",t),
        colorbar=:none)
        end
        mp4(anim,"2dh.mp4")

anim = @animate for t in range(0,3,length=81)
            surface(x,y,sol(t).c',color=:redsblues,
            xaxis=(L"x",(-3,3)),yaxis=(L"y",(-3,3)),zlims=(1,3),clims=(1,3),
            title=@sprintf("TF2Dc, t=%.3f",t),
            colorbar=:none)
            end
            mp4(anim,"2dc.mp4")

anim = @animate for t in range(0,3,length=81)
            surface(x,y,sol(t).p',color=:redsblues,
            xaxis=(L"x",(-3,3)),yaxis=(L"y",(-3,3)),zlims=(-2,1),clims=(-2,1),
            title=@sprintf("TF2Dp, t=%.3f",t),
            colorbar=:none)
            end
            mp4(anim,"2dp.mp4")

anim = @animate for t in range(0,3,length=81)
    surface(x,y,I0*((-exp.(-Φ*sol(t).f'.*sol(t).h')).+1)./((sol(t).f'.^2).+1),color=:redsblues,
    xaxis=(L"x",(-3,3)),yaxis=(L"y",(-3,3)),zlims=(0,1),clims=(0,1),
    title=@sprintf("TF2DI, t=%.3f",t),
    colorbar=:none)
    end
    mp4(anim,"2dFI.mp4")