function fourierf(m, n;
    tspan=(0.0,1.0),
    tol=1e-8
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
    vmax = result[2]+result[4]
    v_max = vmax*1e-6/60
    ell = (sigma_0/mu/v_max)^(1/4)*d
    eps = d/ell
    Pc = (p0*vw*c0)/v_max
    Pec = (v_max*ell)/(eps*D0)
    Pecf = (v_max*ell)/(eps*Df)
    invPec = 1/Pec
    invPecf = 1/Pecf
    vb = 1/vmax
    f_x = Dx*f 
    f_y = f*Dy'
   # h, p, c = unpack(sol_hpc(t), (m,n)) 
   #hump = [0.5*exp(-((x-0.6)/xw)^2/2)*exp(-(y/yw)^2/2)+exp(-((x+0.6)/xw)^2/2)*exp(-(y/yw)^2/2) for x in x, y in y ] 
    h,p,c=unpack([Qh*pod_sol(t)[1:sh*2];Qp*pod_sol(t)[2*sh+1:2*(sh+sp)];Qc*pod_sol(t)[2*(sh+sp)+1:2*(sh+sp+sc)]],(m,n)) 
    hump1 = [ 0.5*exp(-((x-0.6)/result[1])^2/2)*exp(-(y/result[1])^2/2) for x in x, y in y ]
    hump2 = [ exp(-((x+0.6)/result[3])^2/2)*exp(-(y/result[3])^2/2) for x in x, y in y ]
    Jval = vb.+(1-vb)*(result[5]*hump1 + result[6]*hump2)
   
   # Jval=vb.+(1-vb)*hump
    ubar = (-h.^2/12) .* (Dx*p)
    vbar = (-h.^2/12) .* (p*Dy')
    osmo = Pc*(c .- 1)
    tmp = Dxx*(h.*f_x) + (h.*f_y)*Dy'
  
    @. df = (invPecf*tmp - osmo*f + Jval*f)/h - (ubar*f_x + vbar*f_y)
end

dfdt = ODEFunction(TF2df)
f0 = ones(m,n)
prob_f = ODEProblem(dfdt, f0, tspan)
return solve(prob_f, reltol=tol, abstol=tol)

end
m=40;n=40;xw=0.5;yw=0.5;vb=0.1;Pc=0.392;invPec=1/6.76;sh=20;sp=30;sc=20;invPecf=1/27.7;
sol_f2 = fourierf(m,n)

d = 4.5e-6
ϵf = 1.75e+7  #Napierian extinction coefficient
fcr = 0.0053   #critical fluorescein concentration
Φ = ϵf * fcr * d
FI = ((-exp.(-Φ * ones(41,41) .* ones(41,41))) .+ 1) ./ ((ones(41,41) .^ 2) .+ 1)
FI[1]
I0 = 1 / FI[1]
t = 0.95;
I1 = I0 * ((-exp.(-Φ * vec(sol_f(t)) .* sol_hpc(t)[1:m*n])) .+ 1) ./ ((vec(sol_f(t)) .^ 2) .+ 1);
I2 = I0 * ((-exp.(-Φ * vec(sol_f2(t)) .* [Qh*pod_sol(t)[1:sh*2];Qp*pod_sol(t)[2*sh+1:2*(sh+sp)];Qc*pod_sol(t)[2*(sh+sp)+1:2*(sh+sp+sc)]][1:m*n])) .+ 1) ./ ((vec(sol_f2(t)) .^ 2) .+ 1);
    
FI1=heatmap(x, y, I1', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none,title="Original");
FI2=heatmap(x, y, I2', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none,title="Optimized");
FI11=heatmap(x, y, I1', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none);
FI22=heatmap(x, y, I2', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none);
FI111=heatmap(x, y, I1', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none);
FI222=heatmap(x, y, I2', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none);
[Qh*pod_sol(t)[1:sh*2];Qp*pod_sol(t)[2*sh+1:2*(sh+sp)];Qc*pod_sol(t)[2*(sh+sp)+1:2*(sh+sp+sc)]][1:m*n]

S = plot(A,B,C,FI1,A2,B2,C2,FI2,A3,B3,C3,FI3,colorh,colorp,colorc,colorI;layout=lay,size=(680,400))
savefig(PP,"twoss.png")

xh = range(0,1,200);
xp=range(-4,4,200);
xc=range(1,3,200);
y10 = [0];
#lay = @layout [a{0.9h}; b]
lay = @layout [a b c; d e f;g h i;j{0.05h} k{0.05h} l{0.05h}]
colorh = heatmap(xh,y10,[x for x in xh, y in y10],colormap=:viridis,colorbar=false,yticks=[]);
colorp = heatmap(xp,y10,[x for x in xp, y in y10],colormap=:redsblues,colorbar=false,yticks=[]);
colorc = heatmap(xc,y10,[x for x in xc, y in y10],colormap=:viridis,colorbar=false,yticks=[]);
colorI = heatmap(xh,y10,[x for x in xh, y in y10],color=range(RGB(0,0,0),RGB(0,1,0)),colorbar=false,yticks=[]);

layout = @layout [a b ;c d; e f;g{0.05h}]
PP = plot(FI1,FI2,FI11,FI22,FI111,FI222,colorI;layout,size=(680,680))

anim = @animate for t in range(0, 1, length=81)
    I1 = I0 * ((-exp.(-Φ * vec(sol_f(t)) .* sol_hpc4(t)[1:m*n])) .+ 1) ./ ((vec(sol_f(t)) .^ 2) .+ 1);
    heatmap(x, y, I1', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none,title="I");
   
end
gif(anim,"I1.gif")


anim = @animate for t in range(0, 1, length=81)
    I1 = I0 * ((-exp.(-Φ * vec(sol_f(t)) .* sol_hpc(t)[1:m*n])) .+ 1) ./ ((vec(sol_f(t)) .^ 2) .+ 1);
    A = heatmap(x, y, I1', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none,title="I");
    I2 = I0 * ((-exp.(-Φ * vec(sol_f2(t)) .* [Qh*pod_sol(t)[1:sh*2];Qp*pod_sol(t)[2*sh+1:2*(sh+sp)];Qc*pod_sol(t)[2*(sh+sp)+1:2*(sh+sp+sc)]][1:m*n])) .+ 1) ./ ((vec(sol_f2(t)) .^ 2) .+ 1);
    B = heatmap(x, y, I2', color=range(RGB(0,0,0),RGB(0,1,0)),aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none);
    layout = @layout [a b]
    plot(A, B; layout)
end
gif(anim,"II.gif")






