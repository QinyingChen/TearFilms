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

diff1 = zeros(250,1);
diff2 = zeros(250,1);
diff3 = zeros(250,1);
evap1 = zeros(250,1);
osmo1 = zeros(250,1);
evap2 = zeros(250,1);
osmo2 = zeros(250,1);
evap3 = zeros(250,1);
osmo3 = zeros(250,1);

alph=4.06e-2;
#vb=0.1;
#Pc=0.392;
#invPec=1/6.76;
#invPecf=1/27.7;
vb=0.05;
Pc=0.196;
invPec=1/13.52;
invPecf=1/55.4;
hump1 = [ exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y ];
hump2 = [ exp(-(x/0.5)^2/2)*exp(-(y/1)^2/2) for x in x, y in y ];
hump3 = [ exp(-(x/0.5)^2/2)*exp(-(y/4)^2/2) for x in x, y in y ];

for i = eachindex(t1)
    H = reshape(sol_hpc(t1[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    P=reshape(sol_hpc(t1[i])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc(t1[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    F = reshape(sol_f(t1[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    J = vb .+ (1-vb)*hump1 .+ alph*P
    c_x = Dx*C
    c_y = C*Dy'
    ubar = (-H.^2/12).*(Dx*P)
    vbar = (-H.^2/12).*(P*Dy')
    #advection = ubar.*c_x + vbar.*c_y
    diffusion = (invPec*(Dx_2*(H.*c_x) + (H.*c_y)*Dy_2'))./H
    evap = (J.*C)./H
    osmo = ((Pc*(C.-1)).*C)./H
    evap1[i] = evap[1,1]
    osmo1[i] = osmo[1,1]
    diff1[i] = diffusion[1,1]
  
end
for i = eachindex(t2)
    H = reshape(sol_hpc2(t2[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    P=reshape(sol_hpc2(t2[i])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc2(t2[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    F = reshape(sol_f2(t2[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    J = vb .+ (1-vb)*hump2 .+ alph*P
    c_x = Dx*C
    c_y = C*Dy'
    ubar = (-H.^2/12).*(Dx*P)
    vbar = (-H.^2/12).*(P*Dy')   
    diffusion = (invPec*(Dx_2*(H.*c_x) + (H.*c_y)*Dy_2'))./H
    evap = (J.*C)./H
    osmo = ((Pc*(C.-1)).*C)./H
    evap2[i] = evap[1,1]
    osmo2[i] = osmo[1,1]

    diff2[i] = diffusion[1,1]
   
  
end
for i = eachindex(t3)
    H = reshape(sol_hpc3(t3[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    P=reshape(sol_hpc3(t3[i])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc3(t3[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    F = reshape(sol_f3(t3[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    J = vb .+ (1-vb)*hump3 .+ alph*P
    c_x = Dx*C
    c_y = C*Dy'
    ubar = (-H.^2/12).*(Dx*P)
    vbar = (-H.^2/12).*(P*Dy')
    diffusion = (invPec*(Dx_2*(H.*c_x) + (H.*c_y)*Dy_2'))./H
    evap = (J.*C)./H
    osmo = ((Pc*(C.-1)).*C)./H
    evap3[i] = evap[1,1]
    osmo3[i] = osmo[1,1]
    diff3[i] = diffusion[1,1]
   
  
end


plot(t1,diff1,label=L"yw=0.5");
plot!(t2,diff2,label=L"yw=1");
A = plot!(t3,diff3,title="diffusion",label=L"yw=4",legend=:bottomleft)


plot(t1,evap1,label=false);
plot!(t2,evap2,label=false);
B = plot!(t3,evap3,title="evaporation",label=false)


plot(t1,evap1-osmo1,label=false);
plot!(t2,evap2-osmo2,label=false);
C=plot!(t3,evap3-osmo3,xlabel=L"t",title="evap-osmo",label=false)


plot(t1,evap1-osmo1,xlabel="t",ylabel="evap-osmo",label="rw=0.5");
plot!(t2,evap2-osmo2,xlabel="t",ylabel="evap-osmo",label="rw=1");
plot!(t3,evap3-osmo3,xlabel="t",ylabel="evap-osmo",label="rw=4");
E = plot!(t4,evap4-osmo4,xlabel="t",ylabel="evap-osmo",label="rw=1e4")

layout = @layout [a ;b; c ]
S=plot(A,B,C; layout,size = (330, 500))

savefig(S,"diffevaposmo.png")

H = reshape(sol_hpc(t1[end])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
P=reshape(sol_hpc(t1[end-3])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
C=reshape(sol_hpc(t1[end-3])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
F = reshape(sol_f(t1[end])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
c_x = Dx*C
c_y = C*Dy'
diffusion = (invPec*(Dx_2*(H.*c_x) + (H.*c_y)*Dy_2'))./H

ubar = (-H.^2/12).*(Dx*P)
vbar = (-H.^2/12).*(P*Dy')
advection = ubar*c_x + vbar*c_y
surface(x,y,diffusion')
contour(x,y,advection')

J = vb .+ (1-vb)*hump1 .+ alpha*P
osmo = @. (Pc*(C.-1)).*C