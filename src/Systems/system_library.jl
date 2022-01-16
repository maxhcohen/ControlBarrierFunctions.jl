"""
	ControlAffineSystem(system::Symbol)

Construct a control affine system from a system symbol. Current supported symbols are:
[:single_integrator2D, :single_integrator3D, :double_integrator1D, :double_integrator2D, 
:double_integrator3D, :adaptive_cruise_control, :simplified_acc]
"""

function ControlAffineSystem(system::Symbol)
	if system == :single_integrator2D 
		return single_integrator2D()
	elseif system == :single_integrator3D 
		return single_integrator3D()
	elseif system == :double_integrator1D
		return double_integrator1D()
	elseif system == :double_integrator2D
		return double_integrator2D()
	elseif system == :double_integrator3D
		return double_integrator3D()
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
	f(x) = vcat(x[m+1:end], zeros(m))
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
	f(x) = vcat(x[m+1:end], zeros(m))
	g(x) = vcat(zeros(m, m), Matrix(1.0I, m, m))

	return ControlAffineSystem(n, m, f, g)
end	

"""
	adaptive_cruise_control(;v0=13.89, M=1650, f0=0.1, f1=5, f2=0.25)

Benchmark adaptive cruise control dynamics.
"""
function adaptive_cruise_control(;v0=13.89, M=1650, f0=0.1, f1=5, f2=0.25)
	# System dimensions
	n = 2
	m = 1

	# Friction
	F(x) = f0 + f1*x + f2*x^2

	# Dynamics
	f(x) = [v0 - x[2], -F(x[2])/M]
	g(x) = [0.0, 1/M]

	return ControlAffineSystem(n, m, f, g)
end

"""
	simplified_acc(v0=13.89)

Simplified version of the benchmark adaptive cruise control dynamics used in CBF papers.
"""
function simplified_acc(v0=13.89)
	n = 2
	m = 1
	f(x) = [v0 - x[2], 0.0]
	g(x) = [0.0, 1.0]

	return ControlAffineSystem(n, m, f, g)
end