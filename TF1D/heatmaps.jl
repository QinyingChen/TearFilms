anim = @animate for t in range(0, 3, length=81)
    H1 = reshape(sol_hpc(t)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightH1 = [(flipc(H1))[2:end-1, :]; H1]
    Hfull1 = [(flipr(RightH1))[:, 2:end-1] RightH1]
    P1 = reshape(sol_hpc(t)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightP1 = [(flipc(P1))[2:end-1, :]; P1]
    Pfull1 = [(flipr(RightP1))[:, 2:end-1] RightP1]
    C1 = reshape(sol_hpc(t)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightC1 = [(flipc(C1))[2:end-1, :]; C1]
    Cfull1 = [(flipr(RightC1))[:, 2:end-1] RightC1]
   
    A = heatmap(X, Y, Hfull1, color=range(RGB(0, 0, 0), RGB(0, 1, 0)), clims=(0, 1), dpi=200,
        title=@sprintf("h, t=%.3f", t),
    )
    B = heatmap(X, Y, Pfull1, color=range(RGB(0, 0, 0), RGB(0, 1, 0)), clims=(0, 1), dpi=200,
        title=@sprintf("P, t=%.3f", t),
    )
    O = heatmap(X, Y, Cfull1, color=range(RGB(0, 0, 0), RGB(0, 1, 0)), clims=(0, 1), dpi=200,
        title=@sprintf("C, t=%.3f", t),
    )
  



    layout = @layout [a b c]
    plot(A, B, O, D; layout)



end

mp4(anim, "heatmap3.mp4")




H = reshape(sol_hpc(1)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightH= [(flipc(H))[2:end-1, :]; H]
    Hfull= [(flipr(RightH))[:, 2:end-1] RightH]
    A = heatmap(X, Y, Hfull, color=:viridis,aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none,title="h")

    H2 = reshape(sol_hpc(2)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightH2= [(flipc(H2))[2:end-1, :]; H2]
    Hfull2= [(flipr(RightH2))[:, 2:end-1] RightH2]
    A2 = heatmap(X, Y, Hfull2, color=:viridis,aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none)

    H3 = reshape(sol_hpc(3)[1:(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightH3= [(flipc(H3))[2:end-1, :]; H3]
    Hfull3= [(flipr(RightH3))[:, 2:end-1] RightH3]
    A3 = heatmap(X, Y, Hfull3, color=:viridis,aspect_ratio=1, clims=(0, 1),dpi=200,colorbar=:none)


    P = reshape(sol_hpc(1)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightP = [(flipc(P))[2:end-1, :]; P]
    Pfull= [(flipr(RightP))[:, 2:end-1] RightP]

    P2 = reshape(sol_hpc(2)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightP2 = [(flipc(P2))[2:end-1, :]; P2]
    Pfull2= [(flipr(RightP2))[:, 2:end-1] RightP2]
    P3 = reshape(sol_hpc(3)[(nx+1)*(ny+1)+1:2*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightP3 = [(flipc(P3))[2:end-1, :]; P3]
    Pfull3= [(flipr(RightP3))[:, 2:end-1] RightP3]

    C1 = reshape(sol_hpc(1)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightC1 = [(flipc(C1))[2:end-1, :]; C1]
    Cfull1 = [(flipr(RightC1))[:, 2:end-1] RightC1]

    C2 = reshape(sol_hpc(2)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightC2 = [(flipc(C2))[2:end-1, :]; C2]
    Cfull2 = [(flipr(RightC2))[:, 2:end-1] RightC2]

    C3 = reshape(sol_hpc(3)[2*(nx+1)*(ny+1)+1:3*(nx+1)*(ny+1)], (nx + 1, ny + 1))
    RightC3= [(flipc(C3))[2:end-1, :]; C3]
    Cfull3 = [(flipr(RightC3))[:, 2:end-1] RightC3]

    B = heatmap(X, Y, Pfull, color=:redsblues,aspect_ratio=1,clims=(-4,4),dpi=200,colorbar=:none,title="p");
    B2 = heatmap(X, Y, Pfull2, color=:redsblues,aspect_ratio=1,clims=(-4,4),dpi=200,colorbar=:none);
    B3 = heatmap(X, Y, Pfull3, color=:redsblues,aspect_ratio=1,clims=(-4,4),dpi=200,colorbar=:none);
    C = heatmap(X, Y, Cfull1, color=:viridis,aspect_ratio=1,clims=(1, 3),dpi=200,colorbar=:none,title="c");
    C2 = heatmap(X, Y, Cfull2, color=:viridis,aspect_ratio=1,clims=(1, 3),dpi=200,colorbar=:none);
    C3 = heatmap(X, Y, Cfull3, color=:viridis,aspect_ratio=1,clims=(1, 3),dpi=200,colorbar=:none);

    layout = @layout [a b c;d e f;g h i]
    plot(A, B, C, A2,B2,C2 , A3,B3,C3; layout)


    xh = range(0,1,200);
    xp=range(-4,4,200);
    xc=range(1,3,200);
    y10 = [0];
    #lay = @layout [a{0.9h}; b]
    lay = @layout [a b c; d e f;g h i;j{0.05h} k{0.05h} l{0.05h}]
    colorh = heatmap(xh,y10,[x for x in xh, y in y10],colormap=:viridis,colorbar=false,yticks=[]);
    colorp = heatmap(xp,y10,[x for x in xp, y in y10],colormap=:redsblues,colorbar=false,yticks=[]);
    colorc = heatmap(xc,y10,[x for x in xc, y in y10],colormap=:viridis,colorbar=false,yticks=[]);
    S = plot(A,B,C,A2,B2,C2,A3,B3,C3,colorh,colorp,colorc;layout=lay)
    savefig(S,"heatmap_spot.png")
    




    
  