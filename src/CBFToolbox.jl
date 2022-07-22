module CBFToolbox

# Julia packages required for this module
using LinearAlgebra
using Zygote
using JuMP
using OSQP
using Plots
using PlotThemes
import DifferentialEquations: ODEProblem, solve

# Abstract types
abstract type System end
abstract type Policy end
abstract type CertificateFunction end

# Export types and functions related to control affine system
export ControlAffineSystem

# Export base simulation functionality
export Simulation

# System definitions
include("control_affine_system.jl")

# Simulator
include("simulation.jl")

# # Export Simulation type
# export Simulation
# export Trajectory

# # Export systems
# export System
# export ControlAffineSystem
# export step!
# export initialize!

# # Export Certificate functions
# export LyapunovFunction
# export ControlLyapunovFunction
# export TimeVaryingCLF
# export BarrierFunction
# export ControlBarrierFunction
# export HighOrderCBF

# # Export control policies
# export Policy
# export FeedbackPolicy
# export CLFQP
# export TimeVaryingCLFQP
# export CBFQP
# export TimeVaryingCBFQP

# # Export utility functions and types
# export CircularObstacle
# export custom_plots
# export custom_theme
# export circle_shape
# export vec2mat

# # System definitions
# include("Systems/control_affine_system.jl")
# include("Systems/feedback_policy.jl")
# include("Systems/system_library.jl")

# # Lyapunov functions
# include("LyapunovFunctions/control_lyapunov_functions.jl")
# include("LyapunovFunctions/time_varying_clfs.jl")
# include("LyapunovFunctions/clf_quadratic_programs.jl")
# include("LyapunovFunctions/time_varying_clfqp.jl")

# # Barrier functions
# include("BarrierFunctions/control_barrier_functions.jl")
# include("BarrierFunctions/cbf_quadratic_programs.jl")
# include("BarrierFunctions/high_order_cbfs.jl")
# include("BarrierFunctions/hocbf_quadratic_programs.jl")
# include("BarrierFunctions/time_varying_cbfqp.jl")

# # Base simulation type
# include("Simulations/simulation.jl")

# # Utilities
# include("Utils/data_utils.jl")
# include("Utils/cbf_utils.jl")
# include("Utils/plot_utils.jl")

end
