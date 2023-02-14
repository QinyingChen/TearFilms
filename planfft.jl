using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings,Printf
using FFTW
using ComponentArrays


m=16;
n=16;
function plan_fderiv(N::Integer, order::Integer=1)
    @assert iseven(N)
    v = ones(N)
    F = plan_rfft(v)
    F⁻¹ = plan_irfft(F*v, N)
    sgn = iseven(order÷2) ? 1 : -1
    if isodd(order)
        mult = sgn * 1im * [(0:N/2-1).^order; 0] 
    else
        mult = sgn * (0:N/2).^order
    end
    return F, F⁻¹, mult
end

    


   function xderiv(v,r) 
    F, F⁻¹, mult = plan_fderiv(m,r)
    k=zeros(m,n)
    for j=1:m
        k[:,j]=F⁻¹ * ( mult .* (F*v[:,j]))
       
      end  
      return k 
    end

 function yderiv(v,r) 
        F, F⁻¹, mult = plan_fderiv(n,r)
        k=zeros(m,n)
        for j=1:n
            k[j,:]=F⁻¹ * ( mult .* (F*v[j,:]))
           
          end  
          return k 
        end

        
hx = 2π / m
    x=-π.+hx*(1:m)
   
    
    hy = 2π / n
    y=-π.+hy*(1:n)
hump=[exp(-(x/0.5)^2/2)*exp(-(y/0.5)^2/2) for x in x, y in y]

function fourier(m,n)
    

    function TF2d(du,u,params,t)
       
        J=(p)->constants.vb.+(1-constants.vb)*hump.+ constants.alpha*p 
         
        Jval = J(u.p)
            
        ubar=(-u.h.^2/12).*xderiv(u.p,1)
        vbar = (-u.h.^2/12).*yderiv(u.p,1)
        c_x = xderiv(u.c,1) 
        c_y = yderiv(u.c,1)
       
        
        osmo = params.Pc*(u.c .- 1)
       
        u_xx = xderiv(u.h,2) + yderiv(u.h,2)
        
        tmp = xderiv(u.h.*ubar,1) + yderiv(u.h.*vbar,1)
        @. du.h = osmo - tmp - Jval
    
        @. du.p = -u_xx - params.A/u.h^3 - u.p
        
        tmp2 = xderiv(u.h.*c_x,1) + yderiv(u.h.*c_y,1)
        @. du.c = (params.invPec*tmp2 - osmo*u.c + Jval*u.c)/u.h - (ubar*c_x + vbar*c_y)
    
    end
        
    M = Diagonal([ones(m*n);zeros(m*n);ones(m*n)])
    f = ODEFunction(TF2d,mass_matrix=M)


    constants = (vb=0.1,alpha=4.06e-2,A=5.5e-3,Pc=0.392,invPec=1/6.76)
    u0 = ComponentArray(h=ones(m,n),p=constants.A*ones(m,n),c=ones(m,n))
    
    tspan=(0.0,3.0)
   
    prob_mm = ODEProblem(f,u0,tspan,constants)
    
    return x,y, solve(prob_mm,QNDF(linsolve=LinSolveGMRES(),autodiff=false),reltol=1e-8,abstol=1e-8)
 
end


@elapsed x,y,sol=fourier(m,n)