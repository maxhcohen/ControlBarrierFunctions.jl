module CBFToolbox

# Julia packages required for this module
using LinearAlgebra
using Zygote
using JuMP
using OSQP
using Plots
import DifferentialEquations: ODEProblem, solve

# Abstract types
abstract type System end
abstract type Policy end
abstract type CertificateFunction end

# Export Simulation type
export Simulation

# Export systems
export System
export ControlAffineSystem

# Export Certificate functions
export LyapunovFunction
export ControlLyapunovFunction
export BarrierFunction
export ControlBarrierFunction
export HighOrderCBF

# Export control policies
export Policy
export FeedbackPolicy
export CLFQP
export CBFQP

# Export utility functions and types
export CircularObstacle
export latexify_plots
export circle_shape
export vec2mat

# Base simulation type
include("Simulations/simulation.jl")

# System definitions
include("Systems/control_affine_system.jl")
include("Systems/feedback_policy.jl")
include("Systems/system_library.jl")

# Lyapunov functions
include("LyapunovFunctions/control_lyapunov_functions.jl")
include("LyapunovFunctions/clf_quadratic_programs.jl")

# Barrier functions
include("BarrierFunctions/control_barrier_functions.jl")
include("BarrierFunctions/cbf_quadratic_programs.jl")
include("BarrierFunctions/high_order_cbfs.jl")
include("BarrierFunctions/hocbf_quadratic_programs.jl")

# Utilities
include("Utils/data_utils.jl")
include("Utils/cbf_utils.jl")
include("Utils/plot_utils.jl")

end
