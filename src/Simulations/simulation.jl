"""
    Trajectory

Struct containing time, state, and control trajectory from a simulation
"""
struct Trajectory
    t
    x
    u
end

"""
    Simulation

Type used to run various simulations.

# Fields
- `t0::Float64`: initial time
- `tf::Float64`: final time
- `dt::Float64`: time-step
- `ts`: t0:dt:tf
"""
struct Simulation
    t0::Float64
    tf::Float64
    dt::Float64
    ts
end

"""
    Simulation(tf::Float64)
    Simulation(tf::Float64, dt::Float64)
    Simulation(t0::Float64, tf::Float64, dt::Float64)

Construct a simulation object from an initial time, final time, and time-step.
"""
Simulation(tf::Float64) = Simulation(0.0, tf, 0.01, 0.0:0.01:tf)
Simulation(tf::Float64, dt::Float64) = Simulation(0.0, tf, dt, 0.0:dt:tf)
Simulation(t0::Float64, tf::Float64, dt::Float64) = Simulation(t0, tf, dt, t0:dt:tf)

"""
    Base.length(sim::Simulation)

Get number of timesteps in simulation.
"""
Base.length(sim::Simulation) = length(sim.ts)

"""
    (sim::Simulation)(Σ::ControlAffineSystem, x0)
    (sim::Simulation)(Σ::ControlAffineSystem, k::FeedbackPolicy, x0)
    (sim::Simulation)(Σ::ControlAffineSystem, k::TimeVaryingFeedbackPolicy, x0)
    (sim::Simulation)(Σ::ControlAffineSystem, k::CLFQP, x0)
    (sim::Simulation)((Σ::ControlAffineSystem, k::CBFQP, x0)

Simulate a trajectory of a control affine system.
"""
function (sim::Simulation)(Σ::ControlAffineSystem, x0)
    xs = Vector{typeof(x0)}(undef, length(sim))
    xs[1] = x0
    for i in 1:(length(sim) - 1)
        t = sim.ts[i]
        x = xs[i]
        u = Σ.m == 1 ? 0.0 : zeros(Σ.m)
        xs[i + 1] = step(Σ, x, u, t, t + sim.dt)
    end

    return Trajectory(sim.ts, vec2mat(xs), missing)
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::FeedbackPolicy, x0)
    xs = Vector{typeof(x0)}(undef, length(sim))
    us = Σ.m == 1 ? zeros(length(sim)) : zeros(Σ.m, length(sim))
    xs[1] = x0
    for i in 1:(length(sim) - 1)
        t = sim.ts[i]
        x = xs[i]
        u = k(x)
        xs[i + 1] = step(Σ, x, u, t, t + sim.dt)
        Σ.m == 1 ? us[i] = u : us[:, i] = u
    end
    Σ.m == 1 ? us[end] = k(xs[end]) : us[:, end] = k(xs[end])

    return Trajectory(sim.ts, vec2mat(xs), us)
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::TimeVaryingFeedbackPolicy, x0)
    xs = Vector{typeof(x0)}(undef, length(sim))
    us = Σ.m == 1 ? zeros(length(sim)) : zeros(Σ.m, length(sim))
    xs[1] = x0
    for i in 1:(length(sim) - 1)
        t = sim.ts[i]
        x = xs[i]
        u = k(x, t)
        xs[i + 1] = step(Σ, x, u, t, t + sim.dt)
        Σ.m == 1 ? us[i] = u : us[:, i] = u
    end
    Σ.m == 1 ? us[end] = k(xs[end], sim.tf) : us[:, end] = k(xs[end], sim.tf)

    return Trajectory(sim.ts, vec2mat(xs), us)
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::CLFQP, x0)
    xs = Vector{typeof(x0)}(undef, length(sim))
    us = Σ.m == 1 ? zeros(length(sim)) : zeros(Σ.m, length(sim))
    xs[1] = x0
    for i in 1:(length(sim) - 1)
        t = sim.ts[i]
        x = xs[i]
        u = k(x)
        xs[i + 1] = step(Σ, x, u, t, t + sim.dt)
        Σ.m == 1 ? us[i] = u : us[:, i] = u
    end
    Σ.m == 1 ? us[end] = k(xs[end]) : us[:, end] = k(xs[end])

    return Trajectory(sim.ts, vec2mat(xs), us)
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::CBFQP, x0)
    xs = Vector{typeof(x0)}(undef, length(sim))
    us = Σ.m == 1 ? zeros(length(sim)) : zeros(Σ.m, length(sim))
    xs[1] = x0
    for i in 1:(length(sim) - 1)
        t = sim.ts[i]
        x = xs[i]
        u = k(x)
        xs[i + 1] = step(Σ, x, u, t, t + sim.dt)
        Σ.m == 1 ? us[i] = u : us[:, i] = u
    end
    Σ.m == 1 ? us[end] = k(xs[end]) : us[:, end] = k(xs[end])

    return Trajectory(sim.ts, vec2mat(xs), us)
end
