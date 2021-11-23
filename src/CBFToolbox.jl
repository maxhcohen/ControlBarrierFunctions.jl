module CBFToolbox

# Julia packages required for this module
using DifferentialEquations
using ForwardDiff
using Convex
using ECOS

# Export types
export ControlAffineSystem
export ControlLyapunovFunction
export ControlBarrierFunction

# Export core functions
export control
export run_sim

# Export various cases of ControlAffineSystem
export single_integrator
export adaptive_cruise_control

# Export utility functions 
export circle_shape

# Source code
include("systems.jl")
include("system_library.jl")
include("control_lyapunov_functions.jl")
include("control_barrier_functions.jl")
include("utils.jl")

end
