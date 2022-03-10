"""
	ControlAffineSystem(system::String)

Construct a control affine system from a string.
"""

function ControlAffineSystem(system::String)
    if system == "single_integrator2D"
        return single_integrator2D()
    elseif system == "single_integrator3D"
        return single_integrator3D()
    elseif system == "double_integrator1D"
        return double_integrator1D()
    elseif system == "double_integrator2D"
        return double_integrator2D()
    elseif system == "double_integrator3D"
        return double_integrator3D()
    elseif system == "unicycle_3D"
        return unicycle_3D()
    elseif system == "unicycle_4D"
        return unicycle_4D()
    end
end

"""
	single_integrator2D()

Construct a ControlAffineSystem with 2D single integrator dynamics.
"""
function single_integrator2D()
    n = 2
    m = 2
    f(x) = zeros(n)
    g(x) = Matrix(1.0I, m, m)

    return ControlAffineSystem(n, m, f, g)
end

"""
	single_integrator3D()

Construct a ControlAffineSystem with 3D single integrator dynamics.
"""
function single_integrator3D()
    n = 3
    m = 3
    f(x) = zeros(n)
    g(x) = Matrix(1.0I, m, m)

    return ControlAffineSystem(n, m, f, g)
end

"""
	double_integrator1D()

Construct a ControlAffineSystem with 1D double integrator dynamics.
"""
function double_integrator1D()
    n = 2
    m = 1
    f(x) = [x[2], 0.0]
    g(x) = [0.0, 1.0]

    return ControlAffineSystem(n, m, f, g)
end

"""
	double_integrator2D()

Construct a ControlAffineSystem with 2D double integrator dynamics.
"""
function double_integrator2D()
    n = 4
    m = 2
    f(x) = vcat(x[(m + 1):end], zeros(m))
    g(x) = vcat(zeros(m, m), Matrix(1.0I, m, m))

    return ControlAffineSystem(n, m, f, g)
end

"""
	double_integrator3D()

Construct a ControlAffineSystem with 3D double integrator dynamics.
"""
function double_integrator3D()
    n = 6
    m = 3
    f(x) = vcat(x[(m + 1):end], zeros(m))
    g(x) = vcat(zeros(m, m), Matrix(1.0I, m, m))

    return ControlAffineSystem(n, m, f, g)
end

"""
    unicycle_3D()

Construct a ControlAffineSystem with 3D unicycle dynamics.
"""
function unicycle_3D()
    n = 3
    m = 2
    f(x) = zeros(n)
    g(x) = [cos(x[3]) 0.0; sin(x[3]) 0.0; 0.0 1.0]

    return ControlAffineSystem(n, m, f, g)
end

"""
    unicycle_4D()

Construct a ControlAffineSystem with 4D unicycle dynamics.
"""
function unicycle_4D()
    n = 4
    m = 2
    f(x) = [x[4] * cos(x[3]), x[4] * sin(x[3]), 0.0, 0.0]
    g(x) = [0.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 1.0]

    return ControlAffineSystem(n, m, f, g)
end
