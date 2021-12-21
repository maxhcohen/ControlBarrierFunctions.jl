# Base functionality for Control Barrier Functions (CBFs)
abstract type BarrierFunction end

"""
    ControlBarrierFunction

Control barrier function h(x) for a safe set C.
# Fields
- `h`: function h(x) that represents the CBF
- `∇h`: function ∇h(x) defining the gradient of the CBF
- `α`: extended class K function α(h(x))
"""
struct ControlBarrierFunction <: BarrierFunction
    h
    ∇h
    α
end

"""
    ControlBarrierFunction(h; α)

Construct a CBF given h(x) and extended class K function α(h(x)).
Extended classK function α() defaults to the identity i.e., α(h)=h.
"""
function ControlBarrierFunction(h; α=r->r)
    ∇h(x) = ForwardDiff.gradient(h,x)'

    return ControlBarrierFunction(h, ∇h, α)
end

"""
    control(x, Σ::ControlAffineSystem, cbf::ControlBarrierFunction)

Solve standard CBF-QP to get safe control input.
---------------------------------------
"""
function control(x, Σ::ControlAffineSystem, CBF::ControlBarrierFunction)
    u = Convex.Variable(Σ.m)
    problem = minimize(
        0.5sumsquares(u),
        [CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF.h(x))]
        )
    Convex.solve!(problem, ECOS.Optimizer, silent_solver=true)
    if Σ.m == 1
        return u.value
    else
        return vec(u.value)
    end
end

"""
    control(x, kd, Σ::ControlAffineSystem, CBF::ControlBarrierFunction)

Filter nominal policy to get safe input by solving CBF-QP.
---------------------------------------
"""
function control(x, kd, Σ::ControlAffineSystem, CBF::ControlBarrierFunction)
    u = Convex.Variable(Σ.m)
    problem = minimize(
        0.5sumsquares(u) - kd'u,
        [CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF.h(x))]
        )
    Convex.solve!(problem, ECOS.Optimizer, silent_solver=true)
    if Σ.m == 1
        return u.value
    else
        return vec(u.value)
    end
end

"""
    control(
        x, 
        Σ::ControlAffineSystem, 
        CBF::ControlBarrierFunction, 
        CLF::ControlLyapunovFunction,
        p::Float64
    )

Solve standard CBF-CLF-QP to get safe control input, where p is the relaxation penalty.
---------------------------------------
"""
function control(
    x, 
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction,
    CLF::ControlLyapunovFunction,
    p::Float64,
)
    u = Convex.Variable(Σ.m)
    δ = Convex.Variable()
    problem = minimize(
        0.5sumsquares(u) + p*sumsquares(δ),
        [CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)),
        CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ]
        )
    Convex.solve!(problem, ECOS.Optimizer, silent_solver=true)
    if Σ.m == 1
        return u.value
    else
        return vec(u.value)
    end
end

"""
    run_sim(t0, tf, dt, x, Σ::ControlAffineSystem, CBF::ControlBarrierFunction)

Simulate an open-loop system under the influence of a CBF-QP based controller.
---------------------------------------
"""
function run_sim(
	t0::Float64, 
	tf::Float64, 
	dt::Float64, 
	x0::Vector{Float64}, 
	Σ::ControlAffineSystem, 
	CBF::ControlBarrierFunction
)
    # Allocate data for system trajectory
	ts = t0:dt:tf
    xs = zeros(Σ.n, length(ts))
    xs[:,1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[:,i]
        u = control(x, Σ, CBF)
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
	CBF::ControlBarrierFunction
)
    # Allocate data for system trajectory
	ts = t0:dt:tf
    xs = zeros(length(ts))
    xs[1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[i]
        u = control(x, Σ, CBF)
        xs[i+1] = step(x, u, t, t+dt, Σ)
    end

    return ts, xs
end

"""
    run_sim(t0, tf, dt, x, k, Σ::ControlAffineSystem, CBF::ControlBarrierFunction)

Simulate system with nominal policy k(x) filtered through a CBF-QP.
---------------------------------------
"""
function run_sim(
	t0::Float64, 
	tf::Float64, 
	dt::Float64, 
	x0::Vector{Float64}, 
	k, 
	Σ::ControlAffineSystem, 
	CBF::ControlBarrierFunction
)
    # Allocate data for system trajectory
	ts = t0:dt:tf
    xs = zeros(Σ.n, length(ts))
    xs[:,1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[:,i]
        u = control(x, k(x), Σ, CBF)
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
	Σ::ControlAffineSystem, 
	CBF::ControlBarrierFunction
)
    # Allocate data for system trajectory
	ts = t0:dt:tf
    xs = zeros(length(ts))
    xs[1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[i]
        u = control(x, k(x), Σ, CBF)
        xs[i+1] = step(x, u, t, t+dt, Σ)
    end

    return ts, xs
end

"""
    run_sim(
        t0, 
        tf, 
        dt, 
        x, 
        Σ::ControlAffineSystem, 
        cbf::ControlBarrierFunction,
        clf::ControlLyapunovFunction,
        p
    )

Simulate system under the influence of a CBF-CLF-QP based controller.
---------------------------------------
"""
function run_sim(
    t0::Float64, 
    tf::Float64, 
    dt::Float64, 
    x0::Vector{Float64}, 
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction,
    CLF::ControlLyapunovFunction,
    p::Float64
)
    # Allocate data for system trajectory
	ts = t0:dt:tf
    xs = zeros(Σ.n, length(ts))
    xs[:,1] = x0

    # Run simulation
    for i in 1:length(ts)-1
        t = ts[i]
        x = xs[:,i]
        u = control(x, Σ, CBF, CLF, p)
        xs[:,i+1] = step(x, u, t, t+dt, Σ)
    end

    return ts, xs
end