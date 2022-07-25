module CBFToolbox

# Julia packages required for this module
using LinearAlgebra
using ForwardDiff
using JuMP
using OSQP
import DifferentialEquations: ODEProblem, solve

# Import various utilities from other packages
import ControlSystems: lqr, are, Continuous

# Abstract types
abstract type System end
abstract type Controller end
abstract type FeedbackController <: Controller end
abstract type CertificateFunction end
abstract type Disturbance end

# Export concrete types
export Simulation
export ControlAffineSystem
export StateFeedbackController
export LQRController
export ControlLyapunovFunction
export CLFSontag
export CLFQuadProg
export ControlBarrierFunction
export CBFQuadProg
export MatchedDisturbance

# Source code
include("simulation.jl")
include("control_affine_system.jl")
include("state_feedback_controller.jl")
include("lqr_controller.jl")
include("control_lyapunov_function.jl")
include("clf_sontag.jl")
include("clf_quad_prog.jl")
include("control_barrier_function.jl")
include("cbf_quad_prog.jl")
include("matched_disturbance.jl")

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
