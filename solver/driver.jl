include("TF2Dsolver.jl")
include("1dstreak.jl")
 # m, n, xw, yw, vb, Pc, invPec
 m=80;n=80;xw=0.5;yw=1;vb=0.1;Pc=0.392;invPec=1/6.76;
x6, y6, sol_hpc6 = fourier(m,n,xw,yw,vb,Pc,invPec);
params=(;m,n,xw,yw,vb,Pc,invPec)

# 1d_streak
m=80;xw=0.5;vb=0.1;Pc=0.392;invPec=1/6.76;invPecf=1/27.7;
params=(;m,xw,vb,Pc,invPec,invPecf)



@load "clean_e1vm10.jld2" x2 y2 sol_hpc2 params
@load "clean_spotvm10.jld2" x y sol_hpc params
@load "1dstreak_vm10.jld2" x sol params

