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