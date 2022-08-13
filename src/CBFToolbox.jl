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
# export LQRController
export ControlLyapunovFunction
export CLFSontag
export CLFQuadProg
export ControlBarrierFunction
export SecondOrderCBF
export CBFQuadProg
export MatchedDisturbance
export InputToStateCLF

# Plotting functions
export plot_vector_field
export plot_vector_field!
export plot_phase_portrait
export plot_phase_portrait!
export plot_circle
export plot_circle!

# Source code
include("simulation.jl")
include("control_affine_system.jl")
include("state_feedback_controller.jl")
# include("lqr_controller.jl")
include("control_lyapunov_function.jl")
include("clf_sontag.jl")
include("clf_quad_prog.jl")
include("control_barrier_function.jl")
include("second_order_cbf.jl")
include("cbf_quad_prog.jl")
include("hocbf_quad_prog.jl")
include("matched_disturbance.jl")
include("input_to_state_clf.jl")
include("plots.jl")

end
