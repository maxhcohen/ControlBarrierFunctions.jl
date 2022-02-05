"""
	FeedbackPolicy

Struct representing a feedback control policy.

# Fields
- `control`: function u = u = control(x) describing the control law.
"""
struct FeedbackPolicy <: Policy
	control
end

"""
	(k::FeedbackPolicy)(x)

Evaluate FeedbackPolicy at state x.
"""
function (k::FeedbackPolicy)(x)
	return k.control(x)
end

"""
    (sim::Simulation)(Σ::ControlAffineSystem, k::FeedbackPolicy, x0)

Simulate the closed-loop trajectory of a control affine system under a FeedbackPolicy from
initial condition x0.
"""
function (sim::Simulation)(Σ::ControlAffineSystem, k::FeedbackPolicy, x0)
	xs = Vector{typeof(x0)}(undef, length(sim.ts))
	xs[1] = x0
	for i in 1:length(sim.ts)-1
		t = sim.ts[i]
		x = xs[i]
		u = k(x)
		xs[i+1] = step(Σ, x, u, t, t + sim.dt)
	end

	return vec2mat(xs)
end

"""
	TimeVaryingFeedbackPolicy

Struct representing a non-stationary feedback control policy.

# Fields
- `control`: function u = control(x, t) describing the control law.
"""
struct TimeVaryingFeedbackPolicy <: Policy
	control
end

"""
	(k::TimeVaryingFeedbackPolicy)(x)

Evaluate TimeVaryingFeedbackPolicy at state x and time t.
"""
function (k::TimeVaryingFeedbackPolicy)(x, t)
	return k.control(x, t)
end

"""
    (sim::Simulation)(Σ::ControlAffineSystem, k::FeedbackPolicy, x0)

Simulate the closed-loop trajectory of a control affine system under a FeedbackPolicy from
initial condition x0.
"""
function (sim::Simulation)(Σ::ControlAffineSystem, k::TimeVaryingFeedbackPolicy, x0)
	xs = Vector{typeof(x0)}(undef, length(sim.ts))
	xs[1] = x0
	for i in 1:length(sim.ts)-1
		t = sim.ts[i]
		x = xs[i]
		u = k(x, t)
		xs[i+1] = step(Σ, x, u, t, t + sim.dt)
	end

	return vec2mat(xs)
end
