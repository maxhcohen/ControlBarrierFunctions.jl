"Abstract Lyapunov function type"
abstract type LyapunovFunction end

"""
    ControlLyapunovFunction

Control Lyapunov function V for a control affine system
# Fields
- `V`: function V(x) that represents the CLF
- `∇V`: function ∇V(x) defining the gradient of the CLF
- `α`: class K function α(V(x))
"""
struct ControlLyapunovFunction <: LyapunovFunction
    V
    ∇V
    α
end

"""
    ControlLyapunovFunction(V; α)

Construct a CLF given V(x) and class K function α(V(x)).
Extended classK function α() defaults to the identity i.e., α(V)=V.
"""
function ControlLyapunovFunction(V; α=r->r)
    ∇V(x) = ForwardDiff.gradient(V,x)'
    return ControlLyapunovFunction(V, ∇V, α)
end

"""
    control(x, Σ::ControlAffineSystem, clf::ControlLyapunovFunction)

Solve standard CLF-QP to get stabilizing control input.
---------------------------------------
"""
function control(x, Σ::ControlAffineSystem, CLF::ControlLyapunovFunction)
    u = Convex.Variable(Σ.m)
    problem = minimize(
        0.5sumsquares(u),
        [CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x))]
        )
    Convex.solve!(problem, ECOS.Optimizer, silent_solver=true)
    if Σ.m == 1
        return u.value
    else
        return vec(u.value)
    end
end

"""
    run_sim(t0, tf, dt, x, Σ::ControlAffineSystem, clf::ControlLyapunovFunction)

Simulate an open-loop system under the influence of a CLF-QP based controller.
---------------------------------------
"""
function run_sim(
	t0::Float64, 
	tf::Float64, 
	dt::Float64, 
	x0::Vector{Float64}, 
	Σ::ControlAffineSystem, 
	CLF::ControlLyapunovFunction
)
    # Allocate data for system trajectory
	ts = t0:dt:tf
    xs = zeros(Σ.n, length(ts))
    xs[:,1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[:,i]
        u = control(x, Σ, CLF)
        xs[:,i+1] = step(x, u, t, t+dt, Σ)
    end

    return ts, xs
end

function run_sim(
	t0::Float64, 
	tf::Float64, 
	dt::Float64, 
	x0::Float64, 
	Σ::ControlAffineSystem, 
	CLF::ControlLyapunovFunction
)
    # Allocate data for system trajectory
	ts = t0:dt:tf
    xs = zeros(Σ.n, length(ts))
    xs[1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[i]
        u = control(x, Σ, CLF)
        xs[i+1] = step(x, u, t, t+dt, Σ)
    end

    return ts, xs
end