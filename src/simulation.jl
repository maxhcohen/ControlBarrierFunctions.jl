"""
    Simulation

Type representing a simulation object, which can be used to simulate the trajectory of a
dynamical system from an initial condition under a specified control policy.

# Fields
- `t0::Real`: Initial simulation time - defaults to zero.
- `tf::Real`: Ending simulation time.
"""
struct Simulation
    t0::Real
    tf::Real
end

# Simulation constructor from simulation end time
Simulation(T::Real) = Simulation(0.0, T)