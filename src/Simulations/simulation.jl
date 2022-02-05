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
    Simulation(t0::Float64, tf::Float64, dt::Float64)

Construct a simulation object from an initial time, final time, and time-step.
"""
function Simulation(t0::Float64, tf::Float64, dt::Float64)
    ts = t0:dt:tf
    return Simulation(t0, tf, dt, ts)
end
