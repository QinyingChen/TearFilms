using Optim
using NLopt
using Optimization, OptimizationNLopt
using Interpolations
using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using Krylov
using LinearSolve
using Plots
using LaTeXStrings,Printf
using JLD2
using TimerOutputs
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

@load "t01_xy.jld2" x y sol_hpc


m2=40;n2=m2;
theo_data = zeros(3*m2*n2,50)
true_data = zeros(3*m2*n2,50)
interp_h = zeros(m2*n2,50)
interp_p = zeros(m2*n2,50)
interp_c = zeros(m2*n2,50)
interp_h2 = zeros(m2*n2,50)
interp_p2 = zeros(m2*n2,50)
interp_c2 = zeros(m2*n2,50)
sh=20;sp=30;sc=20;
p0 = 12.1e-6; vw = 1.8e-5;c0 = 300;sigma_0 = 0.045;mu = 1.3e-3; d = 4.5e-6;Df = 0.39e-9;D0 = 1.6e-9;
function obj(params)

    # 1d spot solution   
    m=160;  
    n=Int(m/2);
    hx = 2π / m
    x = -π .+ hx*(1:m)
    x = x[n:end]  
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    Dx = Dx/sqrt(2)
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    Dxx = Dxx/(2*sqrt(2))
    hx = 2*sqrt(2)π / m
    x = -sqrt(2)π .+ hx*(1:m)
    x = x[n:end]  

# Modify Dx to T*Dx*K
K=zeros(m,n+1)
T = zeros(n+1,m)
for i = 1:Int(m/2)
    K[m+1-i,n+2-i]=1
    K[i,n+1-i]=1    
end
for i = 1:n+1
    T[n+2-i,m+1-i]=1
end
K1 = zeros(m,n+1)
for i = 1:n
    K1[m+1-i,n+2-i]=-1
    K1[m,n+1]=1
    K1[i,n+1-i]=1   
end
Dx_2 = -(T*Dx*K1)
Dx = T*Dx*K
Dxx = T*Dxx*K
  


function TF1d(dv,v,p,t)  
      v_max = p[2]*1e-6/60
        ell = (sigma_0/mu/v_max)^(1/4)*d
        eps = d/ell
        Pc = (p0*vw*c0)/v_max
        Pec = (v_max*ell)/(eps*D0)
        invPec = 1/Pec
        vb = 1/p[2]
    h=v[1:n+1]
    P=v[n+2:2(n+1)]
    c=v[2(n+1)+1:3(n+1)]
 
    c_x=Dx*c

    
    hump = [ exp(-(x/p[1])^2/2) for x in x]
  
    Jval=vb.+(1-vb)*hump
    ubar=((-h.^2/12).*(Dx*P))

    h_xx = [(Dxx*h)[1];(Dx_2*(x.*(Dx*h)))[2:end]./x[2:end]]
    tmp = [(Dx_2*(h.*ubar))[1];(Dx_2*(x.*(h.*ubar)))[2:end]./x[2:end]]
    tmp2 = [(Dx_2*(h.*c_x))[1];(Dx_2*(x.*(h.*c_x)))[2:end]./x[2:end]]
    

    dv[1:n+1].= Pc*(c.-1)-tmp.-Jval
   
    dv[n+2:2(n+1)].=-h_xx - P
 
    dv[2(n+1)+1:3(n+1)].=(invPec*tmp2-Pc*(c.-1).*c+Jval.*c)./h-(ubar.*c_x)

  
end

M=Diagonal([ones(n+1);zeros(n+1);ones(n+1)])
u0=[ones(n+1);[(Dxx*ones(n+1))[1];(Dx_2*(x.*(Dx*ones(n+1))))[2:end]./x[2:end]];ones(n+1)]
tspan=(0.0,0.5)

sol = ODESolution[]
      #params = [0.5,10,0.6,10]
      for i = 1:2
      solution = solve(ODEProblem(ODEFunction(TF1d, mass_matrix=M), u0, tspan,params[2*i-1:2*i]),reltol=1e-6, abstol=1e-6)
      push!(sol,solution)
      end

sol1 = sol[1]
sol2 = sol[2]


 
        # Interpolate to the rectangular grid
        tdata = range(tspan[1],tspan[2],50)
        for T in eachindex(tdata)
            h = sol1(tdata[T])[1:n+1]
            p = sol1(tdata[T])[n+2:2(n+1)]
            c = sol1(tdata[T])[2(n+1)+1:3(n+1)]
            h2 = sol2(tdata[T])[1:n+1]
            p2 = sol2(tdata[T])[n+2:2(n+1)]
            c2 = sol2(tdata[T])[2(n+1)+1:3(n+1)]
            itp_h = interpolate((x,),h, Gridded(Linear()))
            itp_p = interpolate((x,),p, Gridded(Linear()))
            itp_c = interpolate((x,),c, Gridded(Linear()))
            itp_h2 = interpolate((x,),h2, Gridded(Linear()))
            itp_p2 = interpolate((x,),p2, Gridded(Linear()))
            itp_c2 = interpolate((x,),c2, Gridded(Linear()))
            yy = range(-pi,pi,m2)
            xx = range(-pi,pi,n2)
            h = [itp_h(sqrt.(xx.^2 .+ yy.^2)) for xx in xx, yy in yy]
            p = [itp_p(sqrt.(xx.^2 .+ yy.^2)) for xx in xx, yy in yy]
            c = [itp_c(sqrt.(xx.^2 .+ yy.^2)) for xx in xx, yy in yy]
            interp_h[:,T] = vec(h)
            interp_p[:,T] = vec(p)
            interp_c[:,T] = vec(c)
            h2 = [itp_h2(sqrt.(xx.^2 .+ yy.^2)) for xx in xx, yy in yy]
            p2 = [itp_p2(sqrt.(xx.^2 .+ yy.^2)) for xx in xx, yy in yy]
            c2 = [itp_c2(sqrt.(xx.^2 .+ yy.^2)) for xx in xx, yy in yy]
            interp_h2[:,T] = vec(h2)
            interp_p2[:,T] = vec(p2)
            interp_c2[:,T] = vec(c2)
        end
        # shift spots
        yy = range(-pi,pi,m2)
        xx = range(-pi,pi,n2)
        index1 = findfirst(x->x > pi-0.6,xx)
        index2 = findfirst(x->x > -pi + 0.6,xx)
        for i in 1:50
            interp_h[:,i] = [interp_h[:,i][index1:end];interp_h[:,i][1:index1-1]]
            interp_p[:,i] = [interp_p[:,i][index1:end];interp_p[:,i][1:index1-1]]
            interp_c[:,i] = [interp_c[:,i][index1:end];interp_c[:,i][1:index1-1]]
            interp_h2[:,i] = [interp_h2[:,i][index2:end];interp_h2[:,i][1:index2-1]]
            interp_p2[:,i] = [interp_p2[:,i][index2:end];interp_p2[:,i][1:index2-1]]
            interp_c2[:,i] = [interp_c2[:,i][index2:end];interp_c2[:,i][1:index2-1]]
        
        end
     
        # POD

        U,σ,V = svd(interp_h)
        Bh=U[:,1:sh]
        
        U,σ,V = svd(interp_p)
        Bp=U[:,1:sp]
        
        U,σ,V = svd(interp_c)
        Bc=U[:,1:sc]

        U,σ,V = svd(interp_h2)
        Bh2=U[:,1:sh]
        
        U,σ,V = svd(interp_p2)
        Bp2=U[:,1:sp]
        
        U,σ,V = svd(interp_c2)
        Bc2=U[:,1:sc]
        Bh_hat = hcat(Bh,Bh2)
        Bp_hat = hcat(Bp,Bp2)
        Bc_hat = hcat(Bc,Bc2)
        Q1,R1 = qr(Bh_hat)
        Q2,R2 = qr(Bp_hat)
        Q3,R3 = qr(Bc_hat)
        Qh = Matrix(Q1)
        Qp = Matrix(Q2)
        Qc = Matrix(Q3)

        hx = 2π / m2
    
    x = (-π .+ hx*(1:m2))
    q(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ q(mod(i-j,m2)) for i in 1:m2, j in 1:m2 ]
    q2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ q2(mod(i-j,m2)) for i in 1:m2, j in 1:m2 ]
    
    hy = 2π / n2
   
    y = (-π .+ hy*(1:n2))
    q3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ q3(mod(i-j,n2)) for i in 1:n2, j in 1:n2 ]
    q4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ q4(mod(i-j,n2)) for i in 1:n2, j in 1:n2 ]
        function PODfun(v, params, t)
            vmax = params[2]+params[4]
            v_max = vmax*1e-6/60
            ell = (sigma_0/mu/v_max)^(1/4)*d
            eps = d/ell
            Pc = (p0*vw*c0)/v_max
            Pec = (v_max*ell)/(eps*D0)
            invPec = 1/Pec
            vb = 1/vmax
            h,p,c = unpack([Qh*v[1:sh*2];Qp*v[sh*2+1:2*(sh+sp)];Qc*v[2*(sh+sp)+1:2*(sh+sp+sc)]], (m2,n2))
            hump1 = [ exp(-((x-0.6)/params[1])^2/2)*exp(-(y/params[1])^2/2) for x in x, y in y ]
            hump2 = [ exp(-((x+0.6)/params[3])^2/2)*exp(-(y/params[3])^2/2) for x in x, y in y ]
            Jval = vb.+(1-vb)*(params[5]*hump1 + params[6]*hump2)
            ubar = (-h.^2/12).*(Dx*p)
            vbar = (-h.^2/12).*(p*Dy')
            c_x = Dx*c 
            c_y = c*Dy'
        
            osmo = Pc*(c .- 1)
            h_lap = Dxx*h + h*Dyy'
            
            tmp = Dx*(h.*ubar)+(h.*vbar)*Dy'
            dh = @. osmo - tmp - Jval
            
            dp = @. -h_lap  - p
            tmp = Dx*(h.*c_x) + (h.*c_y)*Dy'
            dc = @. (invPec*tmp - osmo*c + Jval*c)/h - (ubar*c_x + vbar*c_y)
            return [Qh'*vec(dh);Qp'*vec(dp);Qc'*vec(dc)]
            
        end
        
        u0 = pack(ones(m2,n2), Dxx*ones(m2,n2)+ones(m2,n2)*Dyy', ones(m2,n2))
        U1 = Qh'*u0[1:m2*n2]
        U2 = Qp'*u0[m2*n2+1:2*m2*n2]
        U3 = Qc'*u0[2*m2*n2+1:end]
        newu = [U1;U2;U3] 
        M2 = cat(Qh',Qp',Qc'; dims=(1,2))
        M3 = cat(Qh,Qp,Qc; dims=(1,2))
        
        M = Diagonal([ones(m2*n2); zeros(m2*n2); ones(m2*n2)])
        M_new = M2*M*M3
        POD = ODEFunction(PODfun, mass_matrix=M_new)
        reducedprob = ODEProblem(POD, newu, (0,1),params)
    
        pod_sol = solve(reducedprob, reltol=1e-6, abstol=1e-6)
        for T in eachindex(tdata)
            theo_data[:,T] = [Qh*pod_sol(tdata[T])[1:sh*2];Qp*pod_sol(tdata[T])[2*sh+1:2*(sh+sp)];Qc*pod_sol(tdata[T])[2*(sh+sp)+1:2*(sh+sp+sc)]]
            true_data[:,T] = sol_hpc(tdata[T])
            end
        return sum((theo_data - true_data).^2)
    end
    
    initial_guess = [0.5,10.0,0.5,10.0,0.5,1.0]
  

    obj2(x,p) =  obj(x)

    ## 
    f = OptimizationFunction(obj2)
    prob = Optimization.OptimizationProblem(f, initial_guess, lb = [0.1,0.1,0.1,0.1,0.1,0.1], ub = [2.0,15.0,2.0,15.0,2.0,2.0],xtol_abs=1e-2)
    result = solve(prob, NLopt.LN_NELDERMEAD(),maxtime=4800.0,callback = my_callback)

    
    function my_callback(params,obj)
        println("Objective value: ", obj)
        println("Current x: ", params)
        return false
    end

