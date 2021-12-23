module CBFToolbox

# Julia packages required for this module
using DifferentialEquations
using Zygote
using LinearAlgebra
using Convex
using ECOS
using Plots

# Export types
export ControlSystem
export ControlAffineSystem
export LyapunovFunction
export ControlLyapunovFunction
export BarrierFunction
export ControlBarrierFunction

# Export core functions
export control
export run_sim
export step

# Export various cases of ControlAffineSystem
export 
    single_integrator,
    adaptive_cruise_control

# Export various CBFs
export
    cbf_obstacle

# Export plot utility functions 
export 
	circle_shape,
	custom_plots,
	custom_colors

# Source code
include("systems.jl")
include("system_library.jl")
include("control_lyapunov_functions.jl")
include("control_barrier_functions.jl")
include("cbf_library.jl")
include("plot_utils.jl")

end
