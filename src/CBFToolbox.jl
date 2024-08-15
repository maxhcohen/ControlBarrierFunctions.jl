module CBFToolbox

# Julia packages required for this module
using LinearAlgebra
using ForwardDiff
using JuMP
using OSQP
using DifferentialEquations

# Export types
export ControlAffineSystem
export ControlBarrierFunction
export ExplicitSafetyFilter
export QPSafetyFilter
export TunableQPSafetyFilter
export SmoothSafetyFilter
export ISSfSmoothSafetyFilter

# Export core functions
export simulate

# Source code
include("control_affine_system.jl")
include("control_barrier_function.jl")
include("safety_filter.jl")
include("explicit_safety_filter.jl")
include("qp_safety_filter.jl")
include("tunable_qp_safety_filter.jl")
include("smooth_safety_filter.jl")
include("issf_smooth_safety_filter.jl")

end
