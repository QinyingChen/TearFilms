using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf

plot(t1,h,xlabel="",ylabel=L"h",title="",label=false,lw = 2,ls=:solid)
plot!(t2,h2,label=false,lw = 2,ls=:dash)
plot!(t3,h3,label=false,lw = 2,ls=:dot)
A = plot!(t4,h4,label=false,lw = 2,guidefontsize=12,ls=:dashdot)

plot(t1,c,xlabel="",ylabel=L"c",title="",label=false,lw = 2,ls=:solid)
plot!(t2,c2,label=false,lw = 2,ls=:dash)
plot!(t3,c3,label=false,lw = 2,ls=:dot)
B = plot!(t4,c4,label=false,lw = 2,guidefontsize=12,ls=:dashdot)

plot(t1,f,xlabel="",ylabel=L"f",title="",label=false,lw = 2,ls=:solid)
plot!(t2,f2,label=false,lw = 2,ls=:dash)
plot!(t3,f3,label=false,lw = 2,ls=:dot)
C = plot!(t4,f4,label=false,lw = 2,guidefontsize=12,ls=:dashdot)

plot(t1,I1,xlabel=L"t",ylabel=L"I",label=L"y_w=0.5",title="",lw = 2,ls=:solid)
plot!(t2,I2,label=L"y_w=1",lw = 2,ls=:dash)
plot!(t3,I3,label=L"y_w=4",lw = 2,ls=:dot)
D = plot!(t4,I4,label="streak",lw = 2,guidefontsize=12,legendfontsize=10,ls=:dashdot)

# Need: a 2*2 plot and a 4*2 plot

# 2*2 plot

layout = @layout [a b;c d ]
Q1=plot(A,B,C,D; layout, size = (550,550),legendfontsize=7)


# 4*2 plot
# first column has 4 linestyles
# second column has 3 linestyles

##
plot(t1,h,xlabel="",ylabel=L"h",title="",label=false,lw = 2,ls=:solid)
plot!(t2,h2,label=false,lw = 2,ls=:dash)
AA = plot!(t3,h3,label=false,lw = 2,guidefontsize=12,ls=:dot)


plot(t1,c,xlabel="",ylabel=L"c",title="",label=false,lw = 2,ls=:solid)
plot!(t2,c2,label=false,lw = 2,ls=:dash)
BB = plot!(t3,c3,label=false,lw = 2,guidefontsize=12,ls=:dot)


plot(t1,f,xlabel="",ylabel=L"f",title="",label=false,lw = 2,ls=:solid)
plot!(t2,f2,label=false,lw = 2,ls=:dash)
CC = plot!(t3,f3,label=false,lw = 2,guidefontsize=12,ls=:dot)


plot(t1,I1,xlabel=L"t",ylabel=L"I",title="",label=L"x_k=1.5",lw = 2,ls=:solid)
plot!(t2,I2,label=L"x_k=0.8",lw = 2,ls=:dash)
DD = plot!(t3,I3,label=L"x_k=0.6",guidefontsize=12,legendfontsize=11,lw = 2,ls=:dot)

layout = @layout [a{0.25h} b;c{0.25h} d;e{0.25h} f; g{0.25h} h ]
Q2=plot(A,AA,B,BB,C,CC,D,DD; layout,size = (800, 800))

savefig(Q2,"sptosk_fixed_xw_hcfi_center_and_connected_spots_center.pdf")




