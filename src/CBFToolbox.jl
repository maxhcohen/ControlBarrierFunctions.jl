module CBFToolbox

# Julia packages required for this module
using LinearAlgebra
using OrdinaryDiffEq
using ForwardDiff
using JuMP
using OSQP
using Plots
using VectorFieldPlots

# Abstract types
abstract type System end
abstract type Controller end
abstract type FeedbackController <: Controller end
abstract type CertificateFunction end
abstract type HighOrderCBF <: CertificateFunction end
abstract type Disturbance end

# Export concrete types
export Simulation
export ControlAffineSystem
export StateFeedbackController
export ControlLyapunovFunction
export CLFSontag
export CLFQuadProg
export ISSCLFQuadProg
export ControlBarrierFunction
export SecondOrderCBF
export CBFQuadProg
export ISSfCBFQuadProg
export MatchedDisturbance

# Plotting functions
export plot_vector_field
export plot_vector_field!
export plot_phase_portrait
export plot_phase_portrait!
export plot_circle
export plot_circle!
export plot_safe_set
export plot_safe_set!
export plot_constraint_set
export plot_constraint_set!

# Source code
include("Systems/simulation.jl")
include("Systems/control_affine_system.jl")
include("Systems/state_feedback_controller.jl")
include("Lyapunov/control_lyapunov_function.jl")
include("Lyapunov/clf_sontag.jl")
include("Lyapunov/clf_quad_prog.jl")
include("Lyapunov/iss_clf_quad_prog.jl")
include("Barrier/control_barrier_function.jl")
include("Barrier/second_order_cbf.jl")
include("Barrier/cbf_quad_prog.jl")
include("Barrier/hocbf_quad_prog.jl")
include("Barrier/issf_cbf_quad_prog.jl")
include("Barrier/hoissfcbf_quad_prog.jl")
include("Systems/matched_disturbance.jl")
include("Plotting/plots.jl")

end
