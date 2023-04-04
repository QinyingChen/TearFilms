using Plots   
m = 40;
nx = Int(m/2);
n = 40;
ny = Int(n/2);
h = 2pi/n; 
T = sol_hpc.t;        # need to load solutions
T2 = sol_hpc2.t;
T3 = sol_hpc3.t;
T4 = sol_hpc4.t;
T5= sol_hpc5.t;
ch = zeros(length(T),1);
ch2 = zeros(length(T2),1);
ch3 = zeros(length(T3),1);
ch4 = zeros(length(T4),1);
ch5 = zeros(length(T5),1);

 
# c*h
for i = eachindex(T)
    H = reshape(sol_hpc(T[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc(T[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    ch[i] = C[1,1]*H[1,1]
end
for i = eachindex(T2)
    H = reshape(sol_hpc2(T2[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc2(T2[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    ch2[i] = C[1,1]*H[1,1]
end
for i = eachindex(T3)
    H = reshape(sol_hpc3(T3[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc3(T3[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    ch3[i] = C[1,1]*H[1,1]
end
for i = eachindex(T4)
    H = reshape(sol_hpc4(T4[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc4(T4[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    ch4[i] = C[1,1]*H[1,1]
end
for i = eachindex(T5)
    H = reshape(sol_hpc5(T5[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc5(T5[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    ch5[i] = C[1,1]*H[1,1]
end
#plot(T5,ch5,xlabel="t",ylabel="ch",label="larger spot",title="N=90")
plot(T,ch,xlabel="t",ylabel="ch",label="spot",title="N=40");
plot!(T2,ch2,xlabel="t",ylabel="ch",label="ellipse1");
plot!(T3,ch3,xlabel="t",ylabel="ch",label="ellipse2");
A_40 = plot!(T4,ch4,xlabel="t",ylabel="ch",label="streak")


# H

h = zeros(length(T),1);
h2 = zeros(length(T2),1);
h3 = zeros(length(T3),1);
h4 = zeros(length(T4),1);
for i = eachindex(T)
    H = reshape(sol_hpc(T[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    h[i] = H[1,1]
end
for i = eachindex(T2)
    H = reshape(sol_hpc2(T2[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    h2[i] = H[1,1]
end
for i = eachindex(T3)
    H = reshape(sol_hpc3(T3[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    h3[i] = H[1,1]
end
for i = eachindex(T4)
    H = reshape(sol_hpc4(T4[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    h4[i] = H[1,1]
end


plot(T,h,xlabel="t",ylabel="h",label="spot",title="N=40");
plot!(T2,h2,xlabel="t",ylabel="h",label="ellipse1");
plot!(T3,h3,xlabel="t",ylabel="h",label="ellipse2");
Ah_40 = plot!(T4,h4,xlabel="t",ylabel="h",label="streak")

 # c
c = zeros(length(T),1);
c2 = zeros(length(T2),1);
c3 = zeros(length(T3),1);
c4 = zeros(length(T4),1);
for i = eachindex(T)
    C=reshape(sol_hpc(T[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    c[i] = C[1,1]
end
for i = eachindex(T2)
    C=reshape(sol_hpc2(T2[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    c2[i] = C[1,1]
end
for i = eachindex(T3)
    C=reshape(sol_hpc3(T3[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    c3[i] = C[1,1]
end
for i = eachindex(T4)
    C=reshape(sol_hpc4(T4[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    c4[i] = C[1,1]
end


plot(T,c,xlabel="t",ylabel="c",label="spot",title="N=40");
plot!(T2,c2,xlabel="t",ylabel="c",label="ellipse1");
plot!(T3,c3,xlabel="t",ylabel="c",label="ellipse2");
Ac_40 = plot!(T4,c4,xlabel="t",ylabel="c",label="streak")


# Mass

#integral = (h^2/4)*(Cfull[1,1]*Hfull[1,1] + Cfull[1,n]*Hfull[1,n] + Cfull[n,1]*Hfull[n,1] + Cfull[n,n]*Hfull[n,n] + 2*(sum(Cfull[2:n-1,1].*Hfull[2:n-1,1]) + sum(Cfull[2:n-1,n].*Hfull[2:n-1,n]) + sum(Cfull[1,2:n-1].*Hfull[1,2:n-1]) + sum(Cfull[n,2:n-1].*Hfull[n,2:n-1])) + 4*sum(Cfull[2:n-1,2:n-1].*Hfull[2:n-1,2:n-1]))
M = zeros(length(T),1);
M2 = zeros(length(T2),1);
M3 = zeros(length(T3),1);
M4 = zeros(length(T4),1);
for i = eachindex(T)
    H = reshape(sol_hpc(T[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc(T[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    RightH = [(flipc(H))[2:end-1,:];H]
Hfull = [(flipr(RightH))[:,2:end-1] RightH]
RightC = [(flipc(C))[2:end-1,:];C]
Cfull = [(flipr(RightC))[:,2:end-1] RightC]
    M[i]=(h^2/4)*(Cfull[1,1]*Hfull[1,1] + Cfull[1,n]*Hfull[1,n] + Cfull[n,1]*Hfull[n,1] + Cfull[n,n]*Hfull[n,n] + 2*(sum(Cfull[2:n-1,1].*Hfull[2:n-1,1]) + sum(Cfull[2:n-1,n].*Hfull[2:n-1,n]) + sum(Cfull[1,2:n-1].*Hfull[1,2:n-1]) + sum(Cfull[n,2:n-1].*Hfull[n,2:n-1])) + 4*sum(Cfull[2:n-1,2:n-1].*Hfull[2:n-1,2:n-1]))
end
for i = eachindex(T2)
    H = reshape(sol_hpc2(T2[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc2(T2[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    RightH = [(flipc(H))[2:end-1,:];H]
Hfull = [(flipr(RightH))[:,2:end-1] RightH]
RightC = [(flipc(C))[2:end-1,:];C]
Cfull = [(flipr(RightC))[:,2:end-1] RightC]
    M2[i]=(h^2/4)*(Cfull[1,1]*Hfull[1,1] + Cfull[1,n]*Hfull[1,n] + Cfull[n,1]*Hfull[n,1] + Cfull[n,n]*Hfull[n,n] + 2*(sum(Cfull[2:n-1,1].*Hfull[2:n-1,1]) + sum(Cfull[2:n-1,n].*Hfull[2:n-1,n]) + sum(Cfull[1,2:n-1].*Hfull[1,2:n-1]) + sum(Cfull[n,2:n-1].*Hfull[n,2:n-1])) + 4*sum(Cfull[2:n-1,2:n-1].*Hfull[2:n-1,2:n-1]))
end
for i = eachindex(T3)
    H = reshape(sol_hpc3(T3[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc3(T3[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    RightH = [(flipc(H))[2:end-1,:];H]
Hfull = [(flipr(RightH))[:,2:end-1] RightH]
RightC = [(flipc(C))[2:end-1,:];C]
Cfull = [(flipr(RightC))[:,2:end-1] RightC]
    M3[i]=(h^2/4)*(Cfull[1,1]*Hfull[1,1] + Cfull[1,n]*Hfull[1,n] + Cfull[n,1]*Hfull[n,1] + Cfull[n,n]*Hfull[n,n] + 2*(sum(Cfull[2:n-1,1].*Hfull[2:n-1,1]) + sum(Cfull[2:n-1,n].*Hfull[2:n-1,n]) + sum(Cfull[1,2:n-1].*Hfull[1,2:n-1]) + sum(Cfull[n,2:n-1].*Hfull[n,2:n-1])) + 4*sum(Cfull[2:n-1,2:n-1].*Hfull[2:n-1,2:n-1]))
end
for i = eachindex(T4)
    H = reshape(sol_hpc4(T4[i])[1:(nx+1)*(ny+1)],(nx+1,ny+1))
    C=reshape(sol_hpc4(T4[i])[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)],(nx+1,ny+1))
    RightH = [(flipc(H))[2:end-1,:];H]
Hfull = [(flipr(RightH))[:,2:end-1] RightH]
RightC = [(flipc(C))[2:end-1,:];C]
Cfull = [(flipr(RightC))[:,2:end-1] RightC]
    M4[i]=(h^2/4)*(Cfull[1,1]*Hfull[1,1] + Cfull[1,n]*Hfull[1,n] + Cfull[n,1]*Hfull[n,1] + Cfull[n,n]*Hfull[n,n] + 2*(sum(Cfull[2:n-1,1].*Hfull[2:n-1,1]) + sum(Cfull[2:n-1,n].*Hfull[2:n-1,n]) + sum(Cfull[1,2:n-1].*Hfull[1,2:n-1]) + sum(Cfull[n,2:n-1].*Hfull[n,2:n-1])) + 4*sum(Cfull[2:n-1,2:n-1].*Hfull[2:n-1,2:n-1]))
end

plot(T,M.-M[1],xlabel="t",ylabel="M(t)-M(0)",title="N=40",label="spot");
plot!(T2,M2.-M2[1],xlabel="t",ylabel="M(t)-M(0)",label="ellipse1");
plot!(T3,M3.-M3[1],xlabel="t",ylabel="M(t)-M(0)",label="ellipse2");
M_40 = plot!(T4,M4.-M4[1],xlabel="t",ylabel="M(t)-M(0)",label="streak")


plot(M_20,M_40,M_60,M_90;legend = :left,layout)
plot(A_20,A_40,A_60,A_90;legend = :left,layout)
plot(Ac_20,Ac_40,Ac_60,Ac_90;legend = :bottomright,layout)
plot(Ah_20,Ah_40,Ah_60,Ah_90;legend = :topright,layout)
