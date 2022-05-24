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
    (sim::Simulation)(Σ::ControlAffineSystem)
    (sim::Simulation)(Σ::ControlAffineSystem, k::FeedbackPolicy)
    (sim::Simulation)(Σ::ControlAffineSystem, k::CLFQP)
    (sim::Simulation)((Σ::ControlAffineSystem, k::CBFQP)

Simulate a trajectory of a control affine system.
"""
function (sim::Simulation)(Σ::ControlAffineSystem)
    initialize!(Σ)
    for i in 1:(length(sim) - 1)
        u = Σ.m == 1 ? 0.0 : zeros(Σ.m)
        step!(Σ, u, sim.dt)
    end
    vec2mat!(Σ)

    return Σ
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::FeedbackPolicy)
    initialize!(Σ)
    for i in 1:(length(sim) - 1)
        u = k(Σ.x)
        step!(Σ, u, sim.dt)
    end
    vec2mat!(Σ)

    return Σ
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::CLFQP)
    initialize!(Σ)
    for i in 1:(length(sim) - 1)
        u = k(Σ.x)
        step!(Σ, u, sim.dt)
    end
    vec2mat!(Σ)

    return Σ
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::CBFQP)
    initialize!(Σ)
    for i in 1:(length(sim) - 1)
        u = k(Σ.x)
        step!(Σ, u, sim.dt)
    end
    vec2mat!(Σ)

    return Σ
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::TimeVaryingCLFQP)
    initialize!(Σ)
    for i in 1:(length(sim) - 1)
        t = sim.ts[i]
        u = k(Σ.x, t)
        step!(Σ, u, sim.dt)
    end
    vec2mat!(Σ)

    return Σ
end

function (sim::Simulation)(Σ::ControlAffineSystem, k::TimeVaryingCBFQP)
    initialize!(Σ)
    for i in 1:(length(sim) - 1)
        t = sim.ts[i]
        u = k(Σ.x, t)
        step!(Σ, u, sim.dt)
    end
    vec2mat!(Σ)

    return Σ
end

