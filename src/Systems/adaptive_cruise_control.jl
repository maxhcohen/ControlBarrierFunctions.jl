"""
    AdaptiveCruiseControl <: ControlAffineSystem

ControlAffineSystem representing the adaptive cruise control (ACC) problem commonly used in
the CBF literature. The defauly parameters used here are taken from

A. D. Ames, J. W. Grizzle, P. Tabuada, "Control barrier function based quadratic programs
with application to adaptive cruise control," Proceedings of the IEEE Conference on 
Decision and Control, 2014.
"""
struct AdaptiveCruiseControl <: ControlAffineSystem
    m
    f0
    f1
    f2
    vd
    v0
end

function AdaptiveCruiseControl()
    m = 1650.0
    f0 = 0.1
    f1 = 5.0
    f2 = 0.25
    vd = 24.0
    v0 = 13.89

    return AdaptiveCruiseControl(m, f0, f1, f2, vd, v0)
end

state_dim(::AdaptiveCruiseControl) = 2
control_dim(::AdaptiveCruiseControl) = 1
degrees_of_freedom(::AdaptiveCruiseControl) = 1
euler_lagrange(::AdaptiveCruiseControl) = true

function drift(Σ::AdaptiveCruiseControl, x)
    # Pull out parameters
    m = Σ.m
    f0 = Σ.f0
    f1 = Σ.f1
    f2 = Σ.f2
    v0 = Σ.v0

    # Pull out states
    v = x[1] 

    # Friction term
    Fr = f0 + f1*v + f2*v^2

    return [-Fr/m, v0 - v]
end

function actuation(Σ::AdaptiveCruiseControl, x)
    m = Σ.m

    return [1/m, 0.0]
end