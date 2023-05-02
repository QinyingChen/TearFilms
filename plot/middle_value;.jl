m = 80;
nx = Int(m/2);
n = 80;
ny = Int(n/2);
hx = 2π / m
#t1 = range(0,5.5,250);
#t2 = range(0,6.8,250);
#t3 = range(0,7.8,250);


t1 = range(0,5,250);
t2 = range(0,6,250);
t3 = range(0,7,250);


h = zeros(250,1);
h2 = zeros(250,1);
h3 = zeros(250,1);

p = zeros(250,1);
p2 = zeros(250,1);
p3 = zeros(250,1);

c = zeros(250,1);
c2 = zeros(250,1);
c3 = zeros(250,1);

f = zeros(250,1);
f2 = zeros(250,1);
f3 = zeros(250,1);

I1 = zeros(250,1);
I2 = zeros(250,1);
I3 = zeros(250,1);

#reshape(sol_hpc(5.6)[1:(nx+1)*(ny+1)],(nx+1,ny+1))

for i = eachindex(t1)
    H = reshape(sol_hpc(t1[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    P=reshape(sol_hpc(t1[i])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc(t1[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    F = reshape(sol_f(t1[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    I = I0 * ((-exp.(-Φ * F .* H)) .+ 1) ./ ((F .^ 2) .+ 1)
    I1[i] = I[1,1]
    h[i] = H[1,1]
    p[i] = P[1,1]
    c[i] = C[1,1]
    f[i] = F[1,1]
end
for i = eachindex(t2)
    H = reshape(sol_hpc2(t2[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    P=reshape(sol_hpc2(t2[i])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc2(t2[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    F = reshape(sol_f2(t2[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    I = I0 * ((-exp.(-Φ * F .* H)) .+ 1) ./ ((F .^ 2) .+ 1)
    I2[i] = I[1,1]
    h2[i] = H[1,1]
    p2[i] = P[1,1]
    c2[i] = C[1,1]
    f2[i] = F[1,1]
end
for i = eachindex(t3)
    H = reshape(sol_hpc3(t3[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    P=reshape(sol_hpc3(t3[i])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc3(t3[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    F = reshape(sol_f3(t3[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    I = I0 * ((-exp.(-Φ * F .* H)) .+ 1) ./ ((F .^ 2) .+ 1)
    I3[i] = I[1,1]
    h3[i] = H[1,1]
    p3[i] = P[1,1]
    c3[i] = C[1,1]
    f3[i] = F[1,1]
end
for i = eachindex(t4)
    H = reshape(sol_hpc4(t4[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    P=reshape(sol_hpc4(t4[i])[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc4(t4[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    F = reshape(sol_f4(t4[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    I = I0 * ((-exp.(-Φ * F .* H)) .+ 1) ./ ((F .^ 2) .+ 1)
    I4[i] = I[1,1]
    h4[i] = H[1,1]
    p4[i] = P[1,1]
    c4[i] = C[1,1]
    f4[i] = F[1,1]
end


plot(t1,h,ylabel=L"h",label=L"yw=0.5");
plot!(t2,h2,ylabel=L"h",label=L"yw=1");
A1 = plot!(t3,h3,ylabel=L"h",label=L"yw=4");


plot(t1,p,ylabel=L"p",legend=false);
plot!(t2,p2,ylabel=L"p",legend=false);
B1 = plot!(t3,p3,ylabel=L"p",legend=false);


plot(t1,c,ylabel=L"c");
plot!(t2,c2,ylabel=L"c");
C1=plot!(t3,c3,ylabel=L"c",legend=false);


plot(t1,f,xlabel=L"t",ylabel=L"f",legend=false);
plot!(t2,f2,xlabel=L"t",ylabel=L"f",legend=false);
D1 = plot!(t3,f3,xlabel=L"t",ylabel=L"f",legend=false)


plot(t1,I1,xlabel=L"t",ylabel=L"I");
plot!(t2,I2,xlabel=L"t",ylabel=L"I");
E1 = plot!(t3,I3,xlabel=L"t",ylabel=L"I",legend=false)


layout = @layout[a;b;c;d]
plot(A1, B1,C1,D1;layout)

K = plot(A1,B1,C1,D1;layout,size = (300, 500))
