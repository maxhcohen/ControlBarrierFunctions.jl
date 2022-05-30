module CBFToolbox

# Julia packages required for this module
using LinearAlgebra
using OrdinaryDiffEq
using Zygote
using ForwardDiff
using JuMP
using OSQP
using Plots
using PlotThemes

# Import various utilities from other packages
import ControlSystems: lqr, are, Continuous

# Base Abstract types
abstract type System end
abstract type Controller end
abstract type FeedbackController <: Controller end
abstract type CertificateFunction end

# Export Systems
export ControlAffineSystem
export SingleIntegrator
export DoubleIntegrator
export InvertedPendulum

# Export functions for ControlAffineSystems
export state_dim
export control_dim
export drift
export actuation
export _f
export _g
export dynamics
export linearize
export integrate

# Export outputs
export ConfigurationError

# Export certificate functions
export ControlLyapunovFunction
export ControlBarrierFunction

# Export simulator
export Simulator

# Export controllers
export LQRController
export FLController
export CLFController
export CBFController
export CBFCLFController

# Export Simulation type
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

# System definitions
include("Systems/control_affine_system.jl")
include("Systems/single_integrator.jl")
include("Systems/double_integrator.jl")
include("Systems/inverted_pendulum.jl")

# Outputs
include("Outputs/output.jl")
include("Outputs/configuration_error.jl")

# Controllers
include("Controllers/lqr_controller.jl")
include("Controllers/fl_controller.jl")
include("Controllers/control_lyapunov_function.jl")
include("Controllers/clf_controller.jl")
include("Controllers/control_barrier_function.jl")
include("Controllers/cbf_controller.jl")
include("Controllers/cbf_clf_controller.jl")

# Simulator
include("simulator.jl")

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
