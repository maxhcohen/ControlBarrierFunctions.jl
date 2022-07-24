"""
    Simulation

Type representing a simulation object, which can be used to simulate the trajectory of a
dynamical system from an initial condition under a specified control policy.

# Fields
- `t0::Float64`: Initial simulation time - defaults to zero.
- `tf::Float64`: Ending simulation time.
"""
struct Simulation
    t0::Float64
    tf::Float64
end

# Simulation constructor from simulation end time
Simulation(T::Float64) = Simulation(0.0, T)
