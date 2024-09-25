# 1D Heat equation 
# Optimize D
@parameters x t
@variables u(..)
Dxx = Differential(x)^2
Dt = Differential(t)

@parameters D

eq = Dt(u(x,t)) ~ D * Dxx(u(x,t))


u0(x) = sin(pi*x)
bcs = [u(x, 0) ~ u0(x),
       u(0, t) ~ 0.0,
       u(1, t) ~ 0.0]

domains = [x ∈ Interval(0.0, 1.0),
           t ∈ Interval(0.0, 1.0)]



u_chain = Lux.Chain(Dense(2, 8, Lux.σ), Dense(8, 8, Lux.σ), Dense(8, 1))  


# D = 1 Exact solution
u_exact(x, t) = exp(-1 * pi^2 * t) * sin(pi * x)

ts = range(0,1,100)
xs = range(0,1,100)

[u_exact(x,t) for t in ts for x in xs]

    function additional_loss(phi, θ, p)
        u_true = [u_exact(x,t) for t in ts for x in xs]
    
      #  u_predict = [first((phi[1]([t, x], θ[:u]))) for t in ts for x in xs]
       
        return sum(abs2, first((phi[1]([t,x], θ[:u]))) - u_exact(x,t) for t in ts for x in xs) 
    end 


discretization = NeuralPDE.PhysicsInformedNN([u_chain],
    NeuralPDE.QuadratureTraining(; abstol = 1e-6, reltol = 1e-6, batch = 200), param_estim = true,
    additional_loss = additional_loss)
@named pde_system = PDESystem(eq, bcs, domains, [t,x], [u(x,t)], [D],defaults = Dict([p .=> 1.0 for p in [D]]))
prob = NeuralPDE.discretize(pde_system, discretization)

callback = function (p, l)
    println("Current loss is: $l")
    return false
end
res = Optimization.solve(prob, BFGS(linesearch = BackTracking()); maxiters = 1000, time_limit = 60.0, callback = callback)



