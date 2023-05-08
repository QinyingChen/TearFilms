using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using Krylov
using LinearSolve
using Plots
using LaTeXStrings,Printf
using JLD2
include("utility.jl")
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

unpack(u, sz) = unpack!( similar(u, sz), similar(u, sz), similar(u, sz), u )


function fourier(m, n;
    tspan=(0.0,15.0),
    solver=QNDF(linsolve=KrylovJL_GMRES()),
    tol=1e-8
    )
    hx = 2π / m
    nx = Int(m/2)
    x = (-π .+ hx*(1:m))[nx:end]
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    ny = Int(n/2)
    y = (-π .+ hy*(1:n))[ny:end]
    entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

    Kx=zeros(m,nx+1)
Tx = zeros(nx+1,m)
for i = 1:nx
    Kx[m+1-i,nx+2-i]=1
    Kx[i,nx+1-i]=1    
end
K1 = zeros(m,nx+1)
for i = 1:nx
    K1[m+1-i,nx+2-i]=-1
    K1[m,nx+1]=1
    K1[i,nx+1-i]=1   
end
for i = 1:nx+1
    Tx[nx+2-i,m+1-i]=1
end
Dx_2 = -(Tx*Dx*K1)
Dx = Tx*Dx*Kx
Dxx = Tx*Dxx*Kx

Ky=zeros(n,ny+1)
Ty = zeros(ny+1,n)
for i = 1:ny
    Ky[n+1-i,ny+2-i]=1
    Ky[i,ny+1-i]=1    
end
K2 = zeros(n,ny+1)
for i = 1:ny
    K2[n+1-i,ny+2-i]=-1
    K2[n,ny+1]=1
    K2[i,ny+1-i]=1   
end
for i = 1:ny+1
    Ty[ny+2-i,n+1-i]=1
end
Dy_2 = -(Ty*Dy*K2)
Dy = Ty*Dy*Ky
Dyy = Ty*Dyy*Ky

    function TF2d(du,u,params,t) 
        h, p, c = unpack(u, (nx+1,ny+1)) 
        hump = [ exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y ]
        J=(p)->constants.vb.+(1-constants.vb)*hump.+ constants.alpha*p 
        Jval = J(p)     
        #Jval=constants.vb.+(1-constants.vb)*hump     
        ubar = (-h.^2/12).*(Dx*p)
        vbar = (-h.^2/12).*(p*Dy')
        c_x = Dx*c 
        c_y = c*Dy'
       
        osmo = params.Pc*(c .- 1)
        h_lap = Dxx*h + h*Dyy'
        
        tmp = Dx_2*(h.*ubar)+(h.*vbar)*Dy_2'
        dh = @. osmo - tmp - Jval
        #dp = @. -h_lap - params.A ./h.^3 - p
        dp = @. -h_lap  - p
        tmp = Dx_2*(h.*c_x) + (h.*c_y)*Dy_2'
        dc = @. (params.invPec*tmp - osmo*c + Jval*c)/h - (ubar*c_x + vbar*c_y)
        pack!(du, dh, dp, dc)
        update!(prog, tspan[2]-t)
    end
    
    M = Diagonal([ones((nx+1)*(ny+1)); zeros((nx+1)*(ny+1)); ones((nx+1)*(ny+1))])
    dudt = ODEFunction(TF2d, mass_matrix=M)
    #hump = [ exp(-x^2-y^2).*exp(-4*(y+x^2)^2) for x in x, y in y ]
    #hump = [ exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y ]
   
    
    #constants = (vb=0.05,alpha=2.87e-2, Pc=0.196, invPec=1/9.56, invPecf=1/39.17)
    constants = (vb=0.2,alpha=5.74e-2, Pc=0.784, invPec=1/4.78, invPecf=1/19.6)
    #u0 = pack(ones(nx+1,ny+1), Dxx*ones(nx+1,ny+1)+ones(nx+1,ny+1)*Dyy'+ constants.A ./ones(nx+1,ny+1), ones(nx+1,ny+1)) 
    u0 = pack(ones(nx+1,ny+1), Dxx*ones(nx+1,ny+1)+ones(nx+1,ny+1)*Dyy', ones(nx+1,ny+1))
    prob_hpc = ODEProblem(dudt, u0, tspan, constants)
    condition(u,t,integrator)=norm(abs.(reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1))) - reshape(u[1:(nx+1)*(nx+1)],(nx+1,nx+1)))>0
   
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
    prog = ProgressThresh(0.0, 0.5)
    if isnothing(solver)
        sol_hpc = solve(prob_hpc, reltol=tol, abstol=tol)
    else
        sol_hpc = solve(prob_hpc, solver, callback=cb,reltol=tol, abstol=tol)
        #sol_hpc = solve(prob_hpc, solver,reltol=tol, abstol=tol)
    end
    update!(prog, 0.0)
       
    return x, y, sol_hpc
end
x3, y3, sol_hpc3= fourier(80,80);
l=(0.045/1.3e-3/(1e-5/60))^(1/4)*4.5e-3;
coeff_l_20 = 0.5^(1/4)
coeff_l_5 = 2^(1/4)
(coeff_l_20)^2*4.06
(coeff_l_5)^2*4.06
6.76*(coeff_l_20)^2*2
(6.76*(coeff_l_5)^2)/2
27.7*(coeff_l_20)^2*2
(27.7*(coeff_l_5)^2)/2

@load "streakvm10.jld2" x4 y4 sol_hpc4
@load "e2vm10.jld2" x3 y3 sol_hpc3
@load "e1vm10.jld2" x2 y2 sol_hpc2
@load "spotvm10.jld2" x y sol_hpc


@load "spotvm20.jld2" x y sol_hpc
@load "e1vm20.jld2" x2 y2 sol_hpc2
@load "e2vm20.jld2" x3 y3 sol_hpc3



@load "spotvm5.jld2" x y sol_hpc
@load "e1vm5.jld2" x2 y2 sol_hpc2
@load "e2vm5.jld2" x3 y3 sol_hpc3




@load "flspot10.jld2" sol_f
@load "fle1v10.jld2" sol_f2
@load "fle2v10.jld2" sol_f3
@load "flstreakv10.jld2" sol_f4

@load "flspot20.jld2" sol_f
@load "fle1v20.jld2" sol_f2
@load "fle2v20.jld2" sol_f3

function fourierf(m, n;
    tspan=(0.0,12.0),
    tol=1e-8
    )
    hx = 2π / m
    nx = Int(m/2)
    x = (-π .+ hx*(1:m))[nx:end]
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    ny = Int(n/2)
    y = (-π .+ hy*(1:n))[ny:end]
    entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

    Kx=zeros(m,nx+1)
Tx = zeros(nx+1,m)
for i = 1:nx
    Kx[m+1-i,nx+2-i]=1
    Kx[i,nx+1-i]=1    
end
K1 = zeros(m,nx+1)
for i = 1:nx
    K1[m+1-i,nx+2-i]=-1
    K1[m,nx+1]=1
    K1[i,nx+1-i]=1   
end
for i = 1:nx+1
    Tx[nx+2-i,m+1-i]=1
end
Dx_2 = -(Tx*Dx*K1)
Dx = Tx*Dx*Kx
Dxx = Tx*Dxx*Kx

Ky=zeros(n,ny+1)
Ty = zeros(ny+1,n)
for i = 1:ny
    Ky[n+1-i,ny+2-i]=1
    Ky[i,ny+1-i]=1    
end
K2 = zeros(n,ny+1)
for i = 1:ny
    K2[n+1-i,ny+2-i]=-1
    K2[n,ny+1]=1
    K2[i,ny+1-i]=1   
end
for i = 1:ny+1
    Ty[ny+2-i,n+1-i]=1
end
Dy_2 = -(Ty*Dy*K2)
Dy = Ty*Dy*Ky
Dyy = Ty*Dyy*Ky
function TF2df(df, f, params, t)
    f_x = Dx*f 
    f_y = f*Dy'
    h, p, c = unpack(sol_hpc3(t), (nx+1,ny+1))  
    hump = [ exp(-(x/0.5)^2/2)*exp(-(y/4)^2/2) for x in x, y in y ]
    Jval=constants.vb.+(1-constants.vb)*hump.+ constants.alpha.*p
    ubar = (-h.^2/12) .* (Dx*p)
    vbar = (-h.^2/12) .* (p*Dy')
    osmo = params.Pc*(c .- 1)
    tmp = Dx_2*(h.*f_x) + (h.*f_y)*Dy_2'
  
    @. df = (params.invPecf*tmp - osmo*f + Jval*f)/h - (ubar*f_x + vbar*f_y)
end
constants = (vb=0.2,alpha=5.74e-2, Pc=0.784, invPec=1/4.78, invPecf=1/19.6)
#constants = (vb=0.05,alpha=2.87e-2, Pc=0.196, invPec=1/9.56, invPecf=1/39.17)

dfdt = ODEFunction(TF2df)
f0 = ones(nx+1,ny+1)
prob_f = ODEProblem(dfdt, f0, tspan, constants)
sol_f = solve(prob_f, reltol=tol, abstol=tol)


return sol_f
end
sol_f3 = fourierf(80,80);

m =80;
nx = Int(m/2);
n = 80;
ny = Int(n/2);
H=reshape(sol_hpc3(10.4)[1:(nx+1)*(ny+1)],(nx+1,ny+1))
heatmap(X, Y, Hfull, color=:viridis, clims=(0, 1))
contour(x,y,H')

P=reshape(sol_hpc3(10.4)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
C=reshape(sol_hpc(2)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
F=reshape(sol_f(2)[1:(nx+1)*(ny+1)],(nx+1,ny+1))

surface(X,Y,Pfull,zlims=(-5,1))
contour(x,y,P')
X = [(-flip(x))[2:end];x[2:end]]
Y = [(-flip(y))[2:end];y[2:end]]

surface(X,Y,Hfull)

Right = [(flipc(H))[2:end-1,:];H]
Hfull = [(flipr(Right))[:,2:end-1] Right] 

RightP = [(flipc(P))[2:end-1,:];P]
Pfull = [(flipr(RightP))[:,2:end-1] RightP] 

RightC = [(flipc(C))[2:end-1,:];C]
Cfull = [(flipr(RightC))[:,2:end-1] RightC]


RightF = [(flipc(F))[2:end-1, :]; F]
Ffull= [(flipr(RightF))[:, 2:end-1] RightF]

anim = @animate for t in range(0, 5.6, length=81)
    H = reshape(sol_hpc(t)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightH = [(flipc(H))[2:end-1, :]; H]
    Hfull = [(flipr(RightH))[:, 2:end-1] RightH]
    P=reshape(sol_hpc(t)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    RightP = [(flipc(P))[2:end-1,:];P]
Pfull = [(flipr(RightP))[:,2:end-1] RightP] 
C=reshape(sol_hpc(t)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
RightC = [(flipc(C))[2:end-1,:];C]
Cfull = [(flipr(RightC))[:,2:end-1] RightC]
    F = reshape(sol_f(t)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightF = [(flipc(F))[2:end-1, :]; F]
    Ffull= [(flipr(RightF))[:, 2:end-1] RightF]
    
    A = surface(X, Y, Hfull, color=:redsblues, dpi=200,
        xaxis=(L"x", (-pi, pi)), yaxis=(L"y", (-pi, pi)), zlims=(0, 1), clims=(0, 1),
        title=@sprintf("H, t=%.3f", t),
        colorbar=:none)
    B = surface(X, Y, Pfull, color=:redsblues, dpi=200,
        xaxis=(L"x", (-pi, pi)), yaxis=(L"y", (-pi, pi)), zlims=(-50, 1), clims=(0, 1),
        title=@sprintf("P, t=%.3f", t),
        colorbar=:none)
    O = surface(X, Y, Cfull, color=:redsblues, dpi=200,
        xaxis=(L"x", (-pi, pi)), yaxis=(L"y", (-pi, pi)), zlims=(0, 4), clims=(0, 1),
        title=@sprintf("C, t=%.3f", t),
        colorbar=:none)
    D = surface(X, Y, Ffull, color=:redsblues, dpi=200,
        xaxis=(L"x", (-pi, pi)), yaxis=(L"y", (-pi, pi)), zlims=(0, 26), clims=(0, 1),
        title=@sprintf("F, t=%.3f", t),
        colorbar=:none)
    layout = @layout [a b; c d]
    plot(A, B, O, D; layout)

end
mp4(anim, "2d_spot_noA.mp4")

anim = @animate for t in range(0, 6, length=81)
    H = reshape(sol_hpc(t)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightH = [(flipc(H))[2:end-1, :]; H]
    Hfull = [(flipr(RightH))[:, 2:end-1] RightH]
    P=reshape(sol_hpc(t)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    RightP = [(flipc(P))[2:end-1,:];P]
Pfull = [(flipr(RightP))[:,2:end-1] RightP] 
C=reshape(sol_hpc(t)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
RightC = [(flipc(C))[2:end-1,:];C]
Cfull = [(flipr(RightC))[:,2:end-1] RightC]
    F = reshape(sol_f(t)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightF = [(flipc(F))[2:end-1, :]; F]
    Ffull= [(flipr(RightF))[:, 2:end-1] RightF]
    
    A = heatmap(X, Y, Hfull, color=:viridis, aspect_ratio=1,clims=(0,1), dpi=200,
        title=@sprintf("H, t=%.3f", t),
    )
    B= heatmap(X, Y, Pfull, color=:redsblues, aspect_ratio=1,clims=(-50,1), dpi=200,
        title=@sprintf("P, t=%.3f", t),
    )
    O = heatmap(X, Y, Cfull, color=:viridis, aspect_ratio=1,clims=(0,4), dpi=200,
        title=@sprintf("C, t=%.3f", t),
    )
   D = heatmap(X, Y, Ffull, color=:viridis, aspect_ratio=1, clims=(0,27), dpi=200,
        title=@sprintf("F, t=%.3f", t),
    )
    layout = @layout [a b; c d]
    plot(A, B, O, D; layout)

end
mp4(anim, "2d_spot_heatmapA.mp4")

anim = @animate for t in range(0, 6, length=81)
    H = reshape(sol_hpc(t)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightH = [(flipc(H))[2:end-1, :]; H]
    Hfull = [(flipr(RightH))[:, 2:end-1] RightH]
    
    F = reshape(sol_f(t)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightF = [(flipc(F))[2:end-1, :]; F]
    Ffull= [(flipr(RightF))[:, 2:end-1] RightF]
    
    heatmap(X, Y, I0 * ((-exp.(-Φ * Ffull .* Hfull)) .+ 1) ./ ((Ffull .^ 2) .+ 1), color=range(RGB(0, 0, 0), RGB(0, 1, 0)), clims=(0, 1), dpi=200,
        title=@sprintf("I, t=%.3f", t),
    )

end
mp4(anim, "2d_spot_heatmapI_A.mp4")

I0 * ((-exp.(-Φ * Ffull .* Hfull)) .+ 1) ./ ((Ffull .^ 2) .+ 1)
    

# change V_max to 2,5,25 
A2 = surface(x,y,H',xlabel="x",ylabel="y",title="h");
B2 = surface(x,y,P',xlabel="x",ylabel="y",title="p");
C2 = surface(x,y,C',xlabel="x",ylabel="y",title="c");
D2 = surface(x,y,F',xlabel="x",ylabel="y",title="f");

layout = @layout [a b; c d]
plot(A2, B2, C2, D2; layout)