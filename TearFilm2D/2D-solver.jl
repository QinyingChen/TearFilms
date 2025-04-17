using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using Krylov
using LinearSolve
using Plots
using LaTeXStrings,Printf
using JLD2

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

# solve hpc

function fourier(m, n,
    xw,
    yw,
    ak,
    vb,
    Pc,
    invPec;
    tspan=(0.0,3.5),
    solver=QNDF(linsolve=KrylovJL_GMRES()),
    tol=1e-6
    )
    hx = 2π / m
    x = (-π .+ hx*(1:m))
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    y = (-π .+ hy*(1:n))
    entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

    
    
    function TF2d(du,u,params,t) 
       
        h, p, c = unpack(u, (m,n)) 
       # hump1 = [exp(-((x)/xw)^2/2)*exp(-((y)/yw)^2/2)  for x in x, y in y ]
       # hump2 = [exp(-((x+0.8)/xw)^2/2)*exp(-((y)/yw)^2/2)  for x in x, y in y ]
       hump1 = [exp(-((x)/xw)^2/2)*exp(-((y)/(yw))^2/2)  for x in x, y in y ]
       #hump2 = [exp(-((x-0.6)/0.5)^2/2)*exp(-((y)/0.5)^2/2)  for x in x, y in y ]
  
       Jval = vb.+(ak-vb)*hump1
        #Jval = vb.+(1-vb)*hump1 .+ (1-vb)*hump2
    
        ubar = (-h.^2/12).*(Dx*p)
        vbar = (-h.^2/12).*(p*Dy')
        c_x = Dx*c 
        c_y = c*Dy'
       
        osmo = Pc*(c .- 1)
        h_lap = Dxx*h + h*Dyy'
        
        tmp = Dx*(h.*ubar)+(h.*vbar)*Dy'
        dh = @. osmo - tmp - Jval
    
        #@timeit to "J" Jval = vb.+(1-vb)*hump 
            
        #@timeit to "h_lap" h_lap = Dxx*h + h*Dyy'
        #@timeit to "tmp" Dx*(h.*ubar)+(h.*vbar)*Dy'
        
        dp = @. -h_lap  - p
        tmp = Dx*(h.*c_x) + (h.*c_y)*Dy'
        dc = @. (invPec*tmp - osmo*c + Jval*c)/h - (ubar*c_x + vbar*c_y)
        pack!(du, dh, dp, dc)
        update!(prog, tspan[2]-t)
    end
    
    M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
    dudt = ODEFunction(TF2d, mass_matrix=M)
   # hump1 = [exp(-((x)/xw)^2/2)*exp(-((y)/yw)^2/2)  for x in x, y in y ]
   # hump2 = [exp(-((x+0.8)/xw)^2/2)*exp(-((y)/yw)^2/2)  for x in x, y in y ]
  # hump1 = [exp(-((x+0.6)/0.5)^2/2)*exp(-((y)/(0.5))^2/2)  for x in x, y in y ]
  # hump2 = [exp(-((x-0.6)/0.5)^2/2)*exp(-((y)/0.5)^2/2)  for x in x, y in y ]
   #hump1 = [exp(-((x-0.6)/0.5)^2/2)*exp(-((y)/0.5)^2/2) for x in x, y in y ]
  # hump2 = [exp(-((x+0.6)/0.5)^2/2)*exp(-((y)/0.5)^2/2) for x in x, y in y ]
  hump1 = [exp(-((x)/xw)^2/2)*exp(-((y)/(yw))^2/2)  for x in x, y in y ]
        Jval = vb.+(ak-vb)*hump1
        #Jval = vb.+(1-vb)*hump1 .+ (1-vb)*hump2
    center = argmax(Jval)
   
    
    u0 = pack(ones(m,n), Dxx*ones(m,n)+ones(m,n)*Dyy', ones(m,n))
    prob_hpc = ODEProblem(dudt, u0, tspan)
    
    condition(u,t,integrator)=(reshape(u[1:m*n],(m,n)))[center] < 1/4.5
    #condition(u,t,integrator)=(reshape(u[1:m*n],(m,n)))[Q[1],Q[2]] < 1/4.5
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
p0 = 12.1e-6; vw = 1.8e-5;c0 = 300;sigma_0 = 0.045;mu = 1.3e-3; d = 4.5e-6;Df = 0.39e-9;D0 = 1.6e-9;
vmax = 10
v_max = vmax*1e-6/60
ell = (sigma_0/mu/v_max)^(1/4)*d
eps = d/ell
Pc = (p0*vw*c0)/v_max
Pec = (v_max*ell)/(eps*D0)
Pecf = (v_max*ell)/(eps*Df)
invPec = 1/Pec
invPecf = 1/Pecf
vb = 1/vmax

m=40;n=40;xw=0.5;yw=0.5;

ak  = 1;
 x, y, sol_hpc= fourier(m,n,xw,yw,ak,vb,Pc,invPec);


 # solve for f
 
 function fourierf(m, n;
    tspan=(0.0,2.4),
    tol=1e-6
    )
    hx = 2π / m
        
        x = (-π .+ hx*(1:m))
        entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
        Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
        entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
        Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
        
        hy = 2π / n
       
        y = (-π .+ hy*(1:n))
        entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
        Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
        entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
        Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

   
function TF2df(df, f, params, t)
   
    f_x = Dx*f 
    f_y = f*Dy'
    h, p, c = unpack(sol_hpc(t), (m,n)) 
   #hump = [0.5*exp(-((x-0.6)/xw)^2/2)*exp(-(y/yw)^2/2)+exp(-((x+0.6)/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ] 
   # h,p,c=unpack([Bh2*pod_sol2(t)[1:sh];Bp2*pod_sol2(t)[sh+1:sh+sp];Bc2*pod_sol2(t)[sh+sp+1:sh+sp+sc]],(m,n)) 
    #h,p,c=unpack([Bh3*pod_sol3(t)[1:2*sh];Bp3*pod_sol3(t)[2*sh+1:2*(sh+sp)];Bc3*pod_sol3(t)[2*(sh+sp)+1:2*(sh+sp+sc)]],(m,n)) 
 # h,p,c=unpack([Bh*pod_sol(t)[1:sh];Bp*pod_sol(t)[sh+1:sh+sp];Bc*pod_sol(t)[sh+sp+1:sh+sp+sc]],(m,n)) 
  #h,p,c=unpack([Bh3*pod_sol3(t)[1:sh];Bp3*pod_sol3(t)[sh+1:sh+sp];Bc3*pod_sol3(t)[sh+sp+1:sh+sp+sc]],(m,n)) 
    hump1 = [ exp(-(x/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ]
    #hump = [ exp(-(x/0.35)^2/2)*exp(-(y/(0.25/0.35))^2/2) for x in x, y in y ]
   # hump = [ exp(-(x/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ]
   # hump1 = [exp(-((x-0.6)/0.5)^2/2)*exp(-((y)/(0.5))^2/2)  for x in x, y in y ]
   # hump2 = [exp(-((x+0.6)/0.5)^2/2)*exp(-((y)/(0.5))^2/2)  for x in x, y in y ]
   #hump1 = [exp(-((x-0.6)/0.5)^2/2)*exp(-((y)/0.5)^2/2) for x in x, y in y ]
   #hump2 = [exp(-((x+0.6)/0.5)^2/2)*exp(-((y)/0.5)^2/2) for x in x, y in y ]
   # Jval = vb.+(1-vb)*hump1 .+ (1-vb)*hump2
   
    Jval=vb.+(ak-vb)*hump1
    ubar = (-h.^2/12) .* (Dx*p)
    vbar = (-h.^2/12) .* (p*Dy')
    osmo = Pc*(c .- 1)
    tmp = Dx*(h.*f_x) + (h.*f_y)*Dy'
  
    @. df = (invPecf*tmp - osmo*f + Jval*f)/h - (ubar*f_x + vbar*f_y)
    update!(prog, tspan[2]-t)
end

dfdt = ODEFunction(TF2df)
f0 = ones(m,n)
prob_f = ODEProblem(dfdt, f0, tspan)
prog = ProgressThresh(0.0, 0.5)
sol = solve(prob_f, reltol=tol, abstol=tol)
update!(prog, 0.0)
return sol
end

sol_f = fourierf(m,n)





# Plotting

using DataInterpolations

t0 = 0.2;
t1 = 0.5;
t2 = 0.8;
t3 = 1.1;

t0 = 0.2;
t1 = 0.8;
t2 = 1.1;
t3 = 1.7;

t0 = 0.2;
t1 = 1.1;
t2 = 1.7;
t3 = 2.2;

interp_x = range(x[20],x[end],300)
H0= reshape(sol_hpc(t0)[1:m*n],(m,n))
H1= reshape(sol_hpc(t1)[1:m*n],(m,n))
H2= reshape(sol_hpc(t2)[1:m*n],(m,n))
H3= reshape(sol_hpc(t3)[1:m*n],(m,n))

H02= reshape(sol_hpc2(t0)[1:m*n],(m,n))
H12= reshape(sol_hpc2(t1)[1:m*n],(m,n))
H22= reshape(sol_hpc2(t2)[1:m*n],(m,n))
H32= reshape(sol_hpc2(t3)[1:m*n],(m,n))

H03 = reshape(sol_hpc3(t0)[1:m*n],(m,n))
H13 = reshape(sol_hpc3(t1)[1:m*n],(m,n))
H23 = reshape(sol_hpc3(t2)[1:m*n],(m,n))
H33 = reshape(sol_hpc3(t3)[1:m*n],(m,n))



plot(interp_x,QuadraticInterpolation(H0[20:end,20],x[20:end])(interp_x),lw=2,label="t=0.2");
plot!(interp_x,QuadraticInterpolation(H1[20:end,20],x[20:end])(interp_x),lw=2,label="t=0.5");
plot!(interp_x,QuadraticInterpolation(H2[20:end,20],x[20:end])(interp_x),lw=2,label="t=0.8");
A = plot!(interp_x,QuadraticInterpolation(H3[20:end,20],x[20:end])(interp_x),xlabel=L"r",ylabel=L"h",lw=2,legend=false,label="t=1.1");



C0=reshape(sol_hpc(t0)[2*m*n+1:3*m*n],(m,n));
C1=reshape(sol_hpc(t1)[2*m*n+1:3*m*n],(m,n));
C2=reshape(sol_hpc(t2)[2*m*n+1:3*m*n],(m,n));
C3=reshape(sol_hpc(t3)[2*m*n+1:3*m*n],(m,n));

C02 = reshape(sol_hpc2(t0)[2*m*n+1:3*m*n],(m,n));
C12=reshape(sol_hpc2(t1)[2*m*n+1:3*m*n],(m,n));
C22=reshape(sol_hpc2(t2)[2*m*n+1:3*m*n],(m,n));
C32=reshape(sol_hpc2(t3)[2*m*n+1:3*m*n],(m,n));

C03 = reshape(sol_hpc3(t0)[2*m*n+1:3*m*n],(m,n));
C13=reshape(sol_hpc3(t1)[2*m*n+1:3*m*n],(m,n));
C23=reshape(sol_hpc3(t2)[2*m*n+1:3*m*n],(m,n));
C33=reshape(sol_hpc3(t3)[2*m*n+1:3*m*n],(m,n));


plot(interp_x,QuadraticInterpolation(C0[20:end,20],x[20:end])(interp_x),lw=2,label=false);
plot!(interp_x,QuadraticInterpolation(C1[20:end,20],x[20:end])(interp_x),lw=2,label=false);
plot!(interp_x,QuadraticInterpolation(C2[20:end,20],x[20:end])(interp_x),lw=2,label=false);
CC1 = plot!(interp_x,QuadraticInterpolation(C3[20:end,20],x[20:end])(interp_x),xlabel="",ylabel="",lw=2,legendfontsize=10,guidefontsize = 14,label=false,title=L"c");


plot(interp_x,QuadraticInterpolation(C02[20:end,20],x[20:end])(interp_x),lw=2,label=false);
plot!(interp_x,QuadraticInterpolation(C12[20:end,20],x[20:end])(interp_x),lw=2,label=false);
plot!(interp_x,QuadraticInterpolation(C22[20:end,20],x[20:end])(interp_x),lw=2,label=false);
CC2 = plot!(interp_x,QuadraticInterpolation(C32[20:end,20],x[20:end])(interp_x),xlabel="",ylabel="",lw=2,label=false);

plot(interp_x,QuadraticInterpolation(C03[20:end,20],x[20:end])(interp_x),lw=2,label=false);
plot!(interp_x,QuadraticInterpolation(C13[20:end,20],x[20:end])(interp_x),lw=2,label=false);
plot!(interp_x,QuadraticInterpolation(C23[20:end,20],x[20:end])(interp_x),lw=2,label=false);
CC3 = plot!(interp_x,QuadraticInterpolation(C33[20:end,20],x[20:end])(interp_x),xlabel=L"r",ylabel="",lw=2,guidefontsize = 14,label=false);

F0= reshape(sol_f(t0)[1:m*n],(m,n));
F1= reshape(sol_f(t1)[1:m*n],(m,n));
F2= reshape(sol_f(t2)[1:m*n],(m,n));
F3= reshape(sol_f(t3)[1:m*n],(m,n));

F02 = reshape(sol_f2(t0)[1:m*n],(m,n));
F12= reshape(sol_f2(t1)[1:m*n],(m,n));
F22= reshape(sol_f2(t2)[1:m*n],(m,n));
F32= reshape(sol_f2(t3)[1:m*n],(m,n));

F03 = reshape(sol_f3(t0)[1:m*n],(m,n));
F13= reshape(sol_f3(t1)[1:m*n],(m,n));
F23= reshape(sol_f3(t2)[1:m*n],(m,n));
F33= reshape(sol_f3(t3)[1:m*n],(m,n));



plot(interp_x,QuadraticInterpolation(F0[20:end,20],x[20:end])(interp_x),lw=2,label="t=0");
plot!(interp_x,QuadraticInterpolation(F1[20:end,20],x[20:end])(interp_x),lw=2,label="t=0.5");
plot!(interp_x,QuadraticInterpolation(F2[20:end,20],x[20:end])(interp_x),lw=2,label="t=1.5");
D = plot!(interp_x,QuadraticInterpolation(F3[20:end,20],x[20:end])(interp_x),xlabel=L"r",ylabel=L"f",lw=2,legend=false,label="t=2.4");


d = 4.5e-6
ϵf = 1.75e+7  #Napierian extinction coefficient
fcr = 0.0053   #critical fluorescein concentration
Φ = ϵf * fcr * d
FI = ((-exp.(-Φ * ones(41,41) .* ones(41,41))) .+ 1) ./ ((ones(41,41) .^ 2) .+ 1)
FI[1]
I0 = 1 / FI[1]

I01 = I0 * ((-exp.(-Φ * F0 .* H0)) .+ 1) ./ ((F0 .^ 2) .+ 1);
I1 = I0 * ((-exp.(-Φ * F1 .* H1)) .+ 1) ./ ((F1 .^ 2) .+ 1);
I2 = I0 * ((-exp.(-Φ * F2 .* H2)) .+ 1) ./ ((F2 .^ 2) .+ 1);
I3 = I0 * ((-exp.(-Φ * F3 .* H3)) .+ 1) ./ ((F3 .^ 2) .+ 1);

I02 = I0 * ((-exp.(-Φ * F02 .* H02)) .+ 1) ./ ((F02 .^ 2) .+ 1);
I12 = I0 * ((-exp.(-Φ * F12 .* H12)) .+ 1) ./ ((F12 .^ 2) .+ 1);
I22 = I0 * ((-exp.(-Φ * F22 .* H22)) .+ 1) ./ ((F22 .^ 2) .+ 1);
I32 = I0 * ((-exp.(-Φ * F32 .* H32)) .+ 1) ./ ((F32 .^ 2) .+ 1);


I03 = I0 * ((-exp.(-Φ * F03 .* H03)) .+ 1) ./ ((F03 .^ 2) .+ 1);
I13 = I0 * ((-exp.(-Φ * F13 .* H13)) .+ 1) ./ ((F13 .^ 2) .+ 1);
I23 = I0 * ((-exp.(-Φ * F23 .* H23)) .+ 1) ./ ((F23 .^ 2) .+ 1);
I33 = I0 * ((-exp.(-Φ * F33 .* H33)) .+ 1) ./ ((F33 .^ 2) .+ 1);



plot(interp_x,QuadraticInterpolation(I01[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=0.2");
plot!(interp_x,QuadraticInterpolation(I1[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=0.5");
plot!(interp_x,QuadraticInterpolation(I2[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=0.8");
II1 = plot!(interp_x,QuadraticInterpolation(I3[20:end,20],x[20:end])(interp_x),xlabel="",ylabel=L"(a)",lw=2,legendfontsize=10,guidefontsize = 14,label=L"t=1.1",title=L"I");

plot(interp_x,QuadraticInterpolation(I02[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=0.2");
plot!(interp_x,QuadraticInterpolation(I12[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=0.8");
plot!(interp_x,QuadraticInterpolation(I22[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=1.1");
II2 = plot!(interp_x,QuadraticInterpolation(I32[20:end,20],x[20:end])(interp_x),xlabel="",ylabel=L"(b)",lw=2,guidefontsize = 14,legendfontsize=10,label=L"t=1.7");

plot(interp_x,QuadraticInterpolation(I03[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=0.2");
plot!(interp_x,QuadraticInterpolation(I13[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=1.1");
plot!(interp_x,QuadraticInterpolation(I23[20:end,20],x[20:end])(interp_x),lw=2,label=L"t=1.7");
II3 = plot!(interp_x,QuadraticInterpolation(I33[20:end,20],x[20:end])(interp_x),xlabel=L"r",ylabel=L"(c)",lw=2,legendfontsize=10,guidefontsize = 14,label=L"t=2.2");

plot(interp_x,QuadraticInterpolation(I23[20:end,20],x[20:end])(interp_x))

plot(interp_x,QuadraticInterpolation(I0[20:end,20],x[20:end])(interp_x),lw=2,label="t=0");
plot!(interp_x,QuadraticInterpolation(I1[20:end,20],x[20:end])(interp_x),lw=2,label="t=0.5");
plot!(interp_x,QuadraticInterpolation(I2[20:end,20],x[20:end])(interp_x),lw=2,label="t=1.5");
E = plot!(interp_x,QuadraticInterpolation(I3[20:end,20],x[20:end])(interp_x),xlabel=L"r",ylabel=L"I",lw=2,legend=false,label="t=2.4");

layout = @layout [a b ;c d ;e f]
S = plot(II1,CC1,II2,CC2,II3,CC3; layout,legendfontsize=10, size = (680,580))



savefig(S,"three_spot_over_r.pdf")

hump1 = [ exp(-(x/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ]

new_Jval=vb.+(ak-vb)*hump1

new_Jval = reshape(result.u,(m,n))

hx = 2π / m
x = (-π .+ hx*(1:m))
entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]

hy = 2π / n
y = (-π .+ hy*(1:n))