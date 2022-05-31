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

end
