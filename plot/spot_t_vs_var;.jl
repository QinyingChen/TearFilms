m = 80;
nx = Int(m/2);
n = 80;
ny = Int(n/2);
t0 = 0;
t1 = 2;
t2 = 4;
t3 = 6;
t4 = 8;
H0= reshape(sol_hpc(t0)[1:(nx+1)*(ny+1)],(nx+1,ny+1))
H2= reshape(sol_hpc(t1)[1:(nx+1)*(ny+1)],(nx+1,ny+1))
H4= reshape(sol_hpc(t2)[1:(nx+1)*(ny+1)],(nx+1,ny+1))
H5= reshape(sol_hpc(t3)[1:(nx+1)*(ny+1)],(nx+1,ny+1))
H6= reshape(sol_hpc(t4)[1:(nx+1)*(ny+1)],(nx+1,ny+1))
plot(x,diag(H0),label="t=0");
plot!(x,diag(H2),label="t=2");
plot!(x,diag(H4),label="t=4");
plot!(x,diag(H5),label="t=5");
A =plot!(x,diag(H6),xlabel=L"r",ylabel=L"h",legend=false)


P0=reshape(sol_hpc3(t0)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1));
P2=reshape(sol_hpc3(t1)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1));
P4=reshape(sol_hpc3(t2)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1));
P5=reshape(sol_hpc3(t3)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1));
P6=reshape(sol_hpc3(t4)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1));
plot(x,diag(P0),label="t=0");
plot!(x,diag(P2),label="t=2");
plot!(x,diag(P4),label="t=4");
plot!(x,diag(P5),label="t=5");
B =plot!(x,diag(P6),xlabel=L"r",ylabel=L"p",legend=false)

C0=reshape(sol_hpc(t0)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1));
C2=reshape(sol_hpc(t1)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1));
C4=reshape(sol_hpc(t2)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1));
C5=reshape(sol_hpc(t3)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1));
C6=reshape(sol_hpc(t4)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1));
plot(x,diag(C0),label="t=0");
plot!(x,diag(C2),label="t=2");
plot!(x,diag(C4),label="t=4");
plot!(x,diag(C5),label="t=5");
C =plot!(x,diag(C6),xlabel=L"r",ylabel=L"c",legend=false)

F0= reshape(sol_f(t0)[1:(nx+1)*(ny+1)],(nx+1,ny+1));
F2= reshape(sol_f(t1)[1:(nx+1)*(ny+1)],(nx+1,ny+1));
F4= reshape(sol_f(t2)[1:(nx+1)*(ny+1)],(nx+1,ny+1));
F5= reshape(sol_f(t3)[1:(nx+1)*(ny+1)],(nx+1,ny+1));
F6= reshape(sol_f(t4)[1:(nx+1)*(ny+1)],(nx+1,ny+1));
plot(x,diag(F0),label="t=0");
plot!(x,diag(F2),label="t=2");
plot!(x,diag(F4),label="t=4");
plot!(x,diag(F5),label="t=5");
D =plot!(x,diag(F6),xlabel=L"r",ylabel=L"f",legend=false)


d = 4.5e-6
ϵf = 1.75e+7  #Napierian extinction coefficient
fcr = 0.0053   #critical fluorescein concentration
Φ = ϵf * fcr * d
FI = ((-exp.(-Φ * ones(41,41) .* ones(41,41))) .+ 1) ./ ((ones(41,41) .^ 2) .+ 1)
FI[1]
I0 = 1 / FI[1]

I_0 = I0 * ((-exp.(-Φ * F0 .* H0)) .+ 1) ./ ((F0 .^ 2) .+ 1);
I2 = I0 * ((-exp.(-Φ * F2 .* H2)) .+ 1) ./ ((F2 .^ 2) .+ 1);
I4 = I0 * ((-exp.(-Φ * F4 .* H4)) .+ 1) ./ ((F4 .^ 2) .+ 1);
I5 = I0 * ((-exp.(-Φ * F5 .* H5)) .+ 1) ./ ((F5 .^ 2) .+ 1);
I6 = I0 * ((-exp.(-Φ * F6 .* H6)) .+ 1) ./ ((F6 .^ 2) .+ 1);
plot(x,diag(I_0),label=L"t=0");
plot!(x,diag(I2),label=L"t=2");
plot!(x,diag(I4),label=L"t=3");
plot!(x,diag(I5),label=L"t=4");
E =plot!(x,diag(I6),xlabel=L"r",ylabel=L"I",label=L"t=5",legend=:outertopright)

layout = @layout [a b c; d e]
plot(A, D, C, B,E; layout)






