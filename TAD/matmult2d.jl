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

function matmult(m, n, tspan=(0.,3.) )
    hx = 2π / m
    x=-π.+hx*(0:m-1)
    entry(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    y=-π.+hy*(0:n-1)
    entry3(k) = k==0 ? 0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]
    
    function TF2d(du,u,allparams,t)
        h, p, c = unpack(u, (m,n))
        params, J = allparams
        Jval = J(x, y, p)
            
        ubar=(-h.^2/12) .* (Dx*p)
        vbar = (-h.^2/12) .* (p*Dy')
        c_x = Dx*c 
        c_y = c*Dy'
        
        osmo = params.Pc*(c .- 1)
        u_xx = Dxx*h + h*Dyy'
        
        tmp = Dx*(h.*ubar) + (h.*vbar)*Dy'
        dh = osmo - tmp - Jval
        dp = @. -u_xx - params.A/h^3 - p
        tmp = Dx*(h.*c_x) + (h.*c_y)*Dy'
        dc = @. (params.invPec*tmp - osmo*c + Jval*c)/h - (ubar*c_x + vbar*c_y)
        pack!(du, dh, dp, dc)

        return du
    end
        
    M = Diagonal([ones(m*n);zeros(m*n);ones(m*n)])
    f = ODEFunction(TF2d,mass_matrix=M)

    hump(x,y)=[exp(-(x/constants.xw)^2/2)*exp(-(y/constants.xw)^2/2) for x in x, y in y]
   
    constants = (vb=0.1, xw=0.5/L, alpha=4.06e-2, A=5.5e-3, Pc=0.392, invPec=1/6.76)
    J = (x,y,p) -> constants.vb .+ (1-constants.vb)*hump(x,y) .+ constants.alpha*p

    u0 = pack(ones(m,n), fill(constants.A,(m,n)), ones(m,n)) 
    prob_mm = ODEProblem( f, u0, tspan, (constants, J) )
    
    return x, y, prob_mm 
end
