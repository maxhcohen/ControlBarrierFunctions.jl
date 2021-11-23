# Collection of common systems used throughout CBFToolbox.jl

"""
    single_integrator()
    
Construct a two-dimensional single integrator ControlAffineSystem of the form x' = u
"""
function single_integrator()
    n = 2
    m = 2
    f(x) = [0.0, 0.0]
    g(x) = [1.0 0.0; 0.0 1.0]
    return ControlAffineSystem(n, m, f, g)
end

"""
    adaptive_cruise_control()
    
Adaptive cruise control example from

A. D. Ames, X. Xu, J. W. Grizzle, and P. Tabuada, "Control barrier function based quadratic
programs for safety critical systems," IEEE Transactions on Automatic Control, vol. 62,
no. 8, 2017.

Keyword args:
 - `M`: mass of ego vehicle 
 - `f0`: friction coefficient
 - `f1`: friction coefficient
 - `f2`: friction coefficient
 - `vl`: velocity of lead vehicle
"""
function adaptive_cruise_control(; M=1650, f0=0.1, f1=5, f2=0.25, vl=13.89)
    n = 2
    m = 1
    Fr(x) = f0 + f1*x[1] + f2*x[2]^2
    f(x) = [-Fr(x)/M, vl - x[1]]
    g(x) = [1/M, 0.0]
    return ControlAffineSystem(n, m, f, g)
end