using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings
using FFTW
using Krylov
using LinearSolve
using BenchmarkTools

m=50;
unvec = z ->reshape(z,m,m);

hx = 2π / m
    x=-π.+hx*(1:m)
    y=-π.+hx*(1:m)
    #entry(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hx / 2)
    #Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    
    #entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    #Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
   u0 = (x,y) -> sin(4*x)*cos(2*y)
      mtx = f-> [f(x,y) for x in x, y in y]
        u_0=mtx(u0)

    
    
    function fft1(v)
        N = length(v)
        v̂ = fft(v)
        ik = 1im * [0:N/2-1; 0; -N/2+1:-1]
        ŵ = ik .* v̂
        return real(ifft(ŵ))
    end

    function fft2(v)
        N = length(v)
        v̂ = fft(v)
        ik = 1im * [0:N/2;-N/2+1:-1]
        ik2 = ik.*ik
        ŵ = ik2 .* v̂
        return real(ifft(ŵ))
    end

  

    function dudt(du,u,p,t)
     U = unvec(u)
     u_x = zeros(m,m);
     u_xx = zeros(m,m);
     u_yy = zeros(m,m);
    for j=1:m
       u_x[:,j]=fft1(U[:,j])
       u_xx[:,j]=fft2(U[:,j])
       u_yy[j,:]=fft2(U[j,:])
     end
      du.= vec(U.*u_x + p*(u_xx + u_yy))
     
      
   end

   function f_jac(J,u,p,t)
    U = unvec(u)
    u_x = zeros(m,m);
    u_xx = zeros(m,m);
    u_yy = zeros(m,m);
    for j=1:m
      u_x[:,j]=fft1(U[:,j])
      u_xx[:,j]=fft2(U[:,j])
      u_yy[j,:]=fft2(U[j,:])
    end
   J .= Diagonal(vec(U))*kron(I(m),u_x) + kron(u_0,u_x) + p*(kron(I(m),u_xx) + kron(u_yy,I(m)))
   # print("hi")
  end
  function f_jvp(Jv,v,u,p,t)
    U = unvec(u)
    u_x = zeros(m,m);
    u_xx = zeros(m,m);
    u_yy = zeros(m,m);
    for j=1:m
      u_x[:,j]=fft1(U[:,j])
      u_xx[:,j]=fft2(U[:,j])
      u_yy[j,:]=fft2(U[j,:])
    end
    Jv .= (Diagonal(vec(U))*kron(I(m),u_x) + kron(u_0,u_x) + p*(kron(I(m),u_xx) + kron(u_yy,I(m))))*vec(v)
    print("hi")
  end

   tspan=(0.0,3.0)  
        f = ODEFunction(dudt) 
        f = ODEFunction(dudt;jac=f_jac) 
        f = ODEFunction(dudt;jvp=f_jvp) 

     
        prob_mm = ODEProblem(f,vec(u_0),tspan,0.1)
       


        @btime sol2 = solve(prob_mm,reltol=1e-8,abstol=1e-8)

#time = @elapsed sol1 = solve(prob_mm,QNDF(linsolve=LinSolveCG(),autodiff=false),reltol=1e-8,abstol=1e-8)

time = @elapsed sol2 = solve(prob_mm,TRBDF2(linsolve=LinSolveFactorize(qr!),autodiff=false),reltol=1e-8,abstol=1e-8)

time = @elapsed sol3 = solve(prob_mm,TRBDF2(linsolve=LinSolveGMRES(),autodiff=false),reltol=1e-8,abstol=1e-8)

time1 = @btime sol = solve(prob_mm,Rosenbrock23(linsolve=LinSolveGMRES(),autodiff=false),reltol=1e-8,abstol=1e-8)
time2 = @btime sol1 = solve(prob_mm,TRBDF2(linsolve=LinSolveGMRES(),autodiff=false),reltol=1e-8,abstol=1e-8)

t=range(0,3,21);
surface(x,y,sol2[5],
xlabel=L"x",ylabel=L"y",zaxis=((-5,5),L"u(x,y)"),
color=:viridis,alpha=0.66,clims=(-4,4),colorbar=:none,
title="graph",dpi=100 )








