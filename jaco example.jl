# u_t = u_xx - u*u_x
using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings
using FFTW
m=16;

hx = 2π / m
    x=-π.+hx*(0:m-1)
    entry(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    u0 = @. exp(-100(x-1)^2) 
    du = @. ones(m)
    eps = 1e-8
    Dx*u0
    Diagonal(u0)*Dx

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
 
    function dudt(u)
        return Dxx*u - u.*(Dx*u)
    end
   
    function jaco(u,v)
       # return (Dxx - 2*Diagonal(u)*Dx)*v
        return fft2(v)-2*Diagonal(u)*fft1(v)
    end
    A=jaco(u0,du)
    C=Dxx*du - 2*u0.*(Dx*du)
    B=(dudt(u0+eps*du)-dudt(u0))*1e6
    dudt(u0+(eps*du))
    dudt(u0)
    (dudt(u0+eps*du)-dudt(u0))

    function heat(du,u,p,t)
       u_x  = Dx*u
       u_xx = Dxx*u
       
       du[1:m]=u_xx - u.*u_x
      
    end

    function f_jac(J,u,p,t)
        v=ones(m)
        J=jaco(u,v)
       
        nothing
      end
       
        tspan=(0.0,3.0)  
        f = ODEFunction(heat;jac=f_jac)    
       
        prob_mm = ODEProblem(f,u0,tspan)

sol = solve(prob_mm,QNDF(linsolve=LinSolveGMRES()),reltol=1e-8,abstol=1e-8)

t=range(0,3,21)
c=hcat([sol(t) for t in range(0,1,21)]...)

label=hcat(["t=$t" for t in t]...)
plot(x,c,label=label)