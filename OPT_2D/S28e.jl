using FFTW
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
using Statistics
using Images
using Unitful
using Optimization, OptimizationMOI, OptimizationOptimJL

cd("C:\\Users\\cqyzxy\\Downloads\\jld\\jld")

function registration(X::AbstractArray{<:Number, 3}, center, radius, search_radius=div(radius, 3))
    nframes = size(X, 3)
    w = -search_radius:search_radius
    δ = CartesianIndices((w, w))
    window = CartesianIndices((-radius:radius, -radius:radius))  # viewing window

    centers = fill(Tuple(center), nframes)
    indices = fill(window .+ CartesianIndex(Tuple(center)), nframes)

    # fixed reference = last frame
    ref_center = CartesianIndex(center)
    Xref = view(X, window .+ ref_center, nframes)

    for frame in nframes-1:-1:1
        center_now = CartesianIndex(centers[frame + 1])

        dif(shift) = norm(view(X, window .+ (ref_center + shift), frame) - Xref)
       
        stay_diff = norm(view(X, window .+ center_now, frame) - view(X, window .+ center_now, frame + 1))
        min_diff, idx = findmin( dif(shift) for shift in δ )
        updated_diff = norm(view(X, window .+ (ref_center + δ[idx]), frame) - view(X, window .+ center_now, frame + 1))
        
        if abs(updated_diff - stay_diff) / stay_diff > 0.1 && updated_diff > stay_diff
            centers[frame] = Tuple(center_now)
            indices[frame] = window .+ center_now
            
        else
            centers[frame] = Tuple(ref_center + δ[idx])
            indices[frame] = window .+ (ref_center + δ[idx])

        end
       
        
    end

    return indices, centers, X
end

# Wrappers
registration(X::AbstractArray{<:Gray}, args...) = registration(Float16.(X), args...)
registration(X::AbstractArray{<:RGB}, args...) = registration(Gray.(X), args...)
registration(X::AbstractVector{<:AbstractMatrix}, args...) = registration(cat(X..., dims=3), args...)

# MANUALLY SELECT TBU REGION...

D = load("28_visit1_t3_box1.jld2")
X = load.(readdir("/Users/cqyzxy/extracted-data/28_/visit1/t3/extracted", join=true));
#x1 = 40; x2 = 182; radius = 35; search_radius = 10;
x1 = 50; x2 = 200; radius = 35; search_radius = 10;
X = X[x1:x2];
data, centers, Y = registration(X, (360,170), radius, search_radius); 

m1 = size(X[end][data[end]])[1]

# smoothing and interpolation
new_G = zeros(size(X[end][data[end]])[1],size(X[end][data[end]])[1],size(Y)[3])


for i = x1:x2
    new_G[:,:,i-(x1-1)] = Y[data[i-(x1-1)],i-(x1-1)]
end



new_GG = zeros(40,40,size(new_G)[3])

x = y = range(-pi,pi,length=size(new_G[:,:,1])[1]);

    
a = -2.6;b = 2.6; k =5;

I2 =@. [((tanh(k*(x-a)) - tanh(k*(x-b)))/2) * ((tanh(k*(y-a)) - tanh(k*(y-b)))/2) for x in x, y in y];

for i = 1:size(new_G)[3]
    filtered_2 = imfilter(new_G[:,:,i], Kernel.gaussian(2))
    new_opt = mean(filtered_2) .+ I2.*(filtered_2 .- mean(filtered_2))
    itp2 = cubic_spline_interpolation((x, y), new_opt)
    hx = 2π / 40
    y_1 = x_1 = (-π .+ hx*(1:40))
    new_GG[:,:,i] = itp2(x_1, y_1)
end

x = y = range(-pi,pi,length=size(new_GG[:,:,1])[1]);


#scale
scale_max = findmax(new_GG[:,:,1])[1]
new_GG = new_GG/scale_max

###
m=40;n=40;sh=20;sp=30;sc=20;


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

p0 = 12.1e-6; vw = 1.8e-5;c0 = 300;sigma_0 = 0.045;mu = 1.3e-3; d = 4.44793e-6;Df = 0.39e-9;D0 = 1.6e-9;


ϵf = 1.75e+7  #Napierian extinction coefficient
fcr = 0.0053   #critical fluorescein concentration unit: M
f0 = 0.20953
f_0 = f0/0.2 #initial fluorescein concentration, 0.2% is critial FL concentration 
Φ = ϵf * fcr * d

# MANUALLY SET THE INITIAL GUESS

params = [ 0.1  # initial guess
  1
  0.1
  1
  0.1
  0
  0.8
  1
]

m=40;n=40;sh=20;sp=30;sc=20;

  function gaussian_ellipse(x, y; xk, yk, fx, fy, e)
    # Focal vector and angle
    f = [fx, fy]
    θ = atan(fy, fx)
    c = norm(f)  # distance from center to focus
    
    # Semi-major axis
    a = c / e
    b = a * sqrt(1 - e^2)

    # Rotation matrix
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    Λ_inv = Diagonal([1 / a^2, 1 / b^2])
    A = R * Λ_inv * R'

  
    Δx = x .- xk
    Δy = y .- yk
    Q = A[1,1]*Δx.^2 .+ 2*A[1,2]*Δx.*Δy .+ A[2,2]*Δy.^2
    return @. exp(-0.5 * Q)
end

function obj(params)
    beta = params[8]
    vmax = 1/params[1]
    v_max = vmax*1e-6/60
    ell = (sigma_0/mu/v_max)^(1/4)*d
    eps = d/ell
    Pc = (p0*vw*c0)/v_max
    Pec = (v_max*ell)/(eps*D0)
    invPec = 1/Pec
    Pecf = (v_max*ell)/(eps*Df)
    invPecf = 1/Pecf
   
    hx = 2π / m
    x = (-π .+ hx*(1:m))
    entry(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hx / 2)
    Dx = [ entry(mod(i-j,m)) for i in 1:m, j in 1:m ]
    entry2(k) = k==0 ? -π^2/3hx^2-1/6 : -(-1)^k/(2*(sin(k*hx/2))^2)
    Dxx = [ entry2(mod(i-j,m)) for i in 1:m, j in 1:m ]
    
    hy = 2π / n
    y = (-π .+ hy*(1:n))
    entry3(k) = k==0 ? 0.0 : 0.5 * (-1)^k * cot(k * hy / 2)
    Dy = [ entry3(mod(i-j,n)) for i in 1:n, j in 1:n ]
    entry4(k) = k==0 ? -π^2/3hy^2-1/6 : -(-1)^k/(2*(sin(k*hy/2))^2)
    Dyy = [ entry4(mod(i-j,n)) for i in 1:n, j in 1:n ]

    X, Y = [xv for xv in x, yv in y], [yv for xv in x, yv in y]
   hump1 = gaussian_ellipse(X, Y; xk=params[5], yk=params[6], fx=params[2], fy=params[3], e=params[7])
   row_diff = hump1[1, :] .- hump1[end, :]
col_diff = hump1[:, 1] .- hump1[:, end]

sse_row = sum(abs2, row_diff)
sse_col = sum(abs2, col_diff)

threshold = 1e-2

if sse_row + sse_col > threshold
    return 1e6
end

if beta*params[1] > params[4]
    return 1e6
end
tdata = 0:0.2:(size(new_GG,3)-1)*0.2

scale_t = (d*1e6/(1/params[1]))*60  # in seconds
tdata = tdata/scale_t

    # full solver on t = [0,0.5]
  
    function TF2d(du,u,params,t) 

       
        h, p, c = unpack(u, (m,n)) 

        Jval = beta*params[1].+(params[4]-beta*params[1])*hump1 

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
        pack!(du, dh, dp, dc)
        
    end
    
    M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
    dudt = ODEFunction(TF2d, mass_matrix=M)
  
    u0 = pack(ones(m,n), Dxx*ones(m,n)+ones(m,n)*Dyy', ones(m,n))
   
    Jval = beta*params[1].+(params[4]-beta*params[1])*hump1 
    condition1(u,t,integrator)=(reshape(u[1:m*n],(m,n)))[argmax(Jval)] < 0.5/(d*1e6)
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition1,affect!)
    if tdata[end] > 0.5
        prob_hpc = ODEProblem(dudt, u0, (0,0.5),params,callback=cb)
    else
        prob_hpc = ODEProblem(dudt, u0, (0,tdata[end]),params,callback=cb)
    end
    sol_hpc2 = solve(prob_hpc, QNDF(linsolve=KrylovJL_GMRES()),reltol=1e-6, abstol=1e-6)

t = range(0,sol_hpc2.t[end],100)
        W1 = zeros(m*n, length(t))
        for j in eachindex(t)
        W1[:,j] = sol_hpc2(t[j])[1:m*n]
        end
        U,σ,V = svd(W1)
        Bh=U[:,1:sh]
        
        for j in eachindex(t)
        W1[:,j] = sol_hpc2(t[j])[m*n+1:m*n*2]
        end

        U,σ,V = svd(W1)
        Bp=U[:,1:sp]
    
        for j in eachindex(t)
        W1[:,j] = sol_hpc2(t[j])[2*m*n+1:3*m*n]
        end
        U,σ,V = svd(W1)
        Bc=U[:,1:sc]

     # POD solver 
function PODfun(v, params, t)
    h,p,c = unpack([Bh*v[1:sh];Bp*v[sh+1:sh+sp];Bc*v[sh+sp+1:sh+sp+sc]], (m,n)) 
   
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
    return [Bh'*vec(dh);Bp'*vec(dp);Bc'*vec(dc)]
    
end
Jval = beta*params[1].+(params[4]-beta*params[1])*hump1 
condition2(v,t,integrator)=(reshape(Bh*v[1:sh],(m,n)))[argmax(Jval)] < 0.5/(d*1e6)
affect2!(integrator) = terminate!(integrator)
cb2 = DiscreteCallback(condition2,affect2!)
u0 = pack(ones(m,n), Dxx*ones(m,n)+ones(m,n)*Dyy', ones(m,n))
U1 = Bh'*u0[1:m*n]
U2 = Bp'*u0[m*n+1:2*m*n]
U3 = Bc'*u0[2*m*n+1:end]
newu = [U1;U2;U3] 
M2 = cat(Bh',Bp',Bc'; dims=(1,2))
M3 = cat(Bh,Bp,Bc; dims=(1,2))

M = Diagonal([ones(m*n); zeros(m*n); ones(m*n)])
M_new = M2*M*M3
POD = ODEFunction(PODfun, mass_matrix=M_new)
reducedprob = ODEProblem(POD, newu, (0,5),params,callback=cb2)

pod_sol = solve(reducedprob,reltol=1e-6, abstol=1e-6)

if pod_sol.t[end] < tdata[end]
    return 1e6
end

# solver for f

function TF2df(df, f, params, t)
   
    f_x = Dx*f 
    f_y = f*Dy'

  h,p,c=unpack([Bh*pod_sol(t)[1:sh];Bp*pod_sol(t)[sh+1:sh+sp];Bc*pod_sol(t)[sh+sp+1:sh+sp+sc]],(m,n)) 
  
  Jval = beta*params[1].+(params[4]-beta*params[1])*hump1 
   
    ubar = (-h.^2/12) .* (Dx*p)
    vbar = (-h.^2/12) .* (p*Dy')
    osmo = Pc*(c .- 1)
    tmp = Dx*(h.*f_x) + (h.*f_y)*Dy'
  
    @. df = (invPecf*tmp - osmo*f + Jval*f)/h - (ubar*f_x + vbar*f_y)
 
end

dfdt = ODEFunction(TF2df)
f0 = f_0*ones(m,n)
prob_f = ODEProblem(dfdt, f0, (0,maximum(pod_sol.t)),params)

solf = solve(prob_f,reltol=1e-6, abstol=1e-6)


theo_data = zeros(m*n,length(tdata))
true_data = zeros(m*n,length(tdata))
initial_value = ((-exp.(-Φ * vec(solf(0)) .* (Bh*pod_sol(0)[1:sh]))) .+ 1) ./ ((vec(solf(0)) .^ 2) .+ 1)
I_0 = 1 / initial_value[1]

    for T in eachindex(tdata)
        true_data[:,T]  = vec(new_GG[:,:,T]);
        theo_data[:,T] = I_0*((-exp.(-Φ * vec(solf(tdata[T])) .* (Bh*pod_sol(tdata[T])[1:sh]))) .+ 1) ./ ((vec(solf(tdata[T])) .^ 2) .+ 1);
    end
    return sum((theo_data[:,1:end] - true_data[:,1:end]).^2)
#
end

obj2(x,p) =  obj(x)

    ## 
    f = OptimizationFunction(obj2)
    #f = OptimizationFunction(obj2, Optimization.AutoFiniteDiff())
    prob = Optimization.OptimizationProblem(f, initial_guess,lb=[0.01,-2.5,-3.5,0.1,-3.0,-3.0,0.1,0.1],ub=[1.0,3.5,4.0,3.0,3.0,3.0,0.99,10]) #v_b,f_x,f_y,a_1,x_k,y_k,e,beta single ellipse
   
    result = solve(prob,NLopt.LN_PRAXIS(),callback = my_callback)

    result = solve(prob,NLopt.LN_NELDERMEAD(),callback = my_callback)

function my_callback(params,obj)
        println("Objective value: ", obj)
        println("Current x: ", params)
        return false
    end


# Optimized parameters
params = [ 0.20694066098747754
  1.007424477777116
  0.17203020180454767
  1.6169082913575206
  0.0711191799577547
 -0.22043687351035351
  0.7794533698746218
  2.4691128670219236
]






