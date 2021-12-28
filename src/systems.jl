"Abstract type used to derive special classes of control systems"
abstract type ControlSystem end

"""
    ControlAffineSystem

Models a nonlinear control affine system of the form 
    x' = f(x) + g(x)u
# Fields
- `n::Int`: state dimension
- `m::Int`: control dimension
- `f`: function f(x) that models the drift dynamics
- `g`: function g(x) that models the control directions
"""
struct ControlAffineSystem <: ControlSystem
    n::Int
    m::Int
    f
    g
end

"""
    step(x, u, t0, tf, Σ::ControlAffineSystem)
    
Integrate the dynamics of a ControlSystemAffine starting at state x under control u over 
from t0 to tf.
"""
function step(x, u, t0, tf, Σ::ControlAffineSystem)
	rhs(x, u, t) = Σ.f(x) + Σ.g(x)u
	sol = DifferentialEquations.solve(ODEProblem(rhs, x, (t0, tf), u))
	
	return sol[:,end]
end

"""
    run_sim(t0, tf, dt, x, k, Σ)

Simulate system under a nominal feedback control policy k(x)
"""
function run_sim(
	t0::Float64, 
	tf::Float64, 
	dt::Float64, 
	x0::Vector{Float64}, 
	k, 
	Σ::ControlAffineSystem
)
    # Allocate data for system trajectories
	ts = t0:dt:tf
    xs = zeros(Σ.n, length(ts))
    xs[:,1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[:,i]
        u = k(x)
		xs[:,i+1] = step(x, u, t, t+dt, Σ)
    end

    return ts, xs
end

function run_sim(
	t0::Float64, 
	tf::Float64, 
	dt::Float64, 
	x0::Float64, 
	k, 
	Σ::ControlAffineSystem
)
    # Allocate data for system trajectories
	ts = t0:dt:tf
    xs = zeros(length(ts))
    xs[1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[i]
        u = k(x)
		xs[i+1] = step(x, u, t, t+dt, Σ)
    end

    return ts, xs
end