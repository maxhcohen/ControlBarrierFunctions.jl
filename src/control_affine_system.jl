"""
    ControlAffineSystem

Control affine system described by the dynamics ``ẋ = f(x) + g(x)u`` where ``x ∈ ℝⁿ`` is the system state, ``u ∈ ℝᵐ`` is the control input and ``f : ℝⁿ → ℝⁿ`` and ``g : ℝⁿ → ℝⁿˣᵐ`` represent the dynamics.

# Fields
- `name::String` : name of system, e.g., "double integrator", "unicycle", "quadrotor", etc.
- `n::Int` : state dimension
- `m::Int` : control dimension
- `f::Function` : drift dynamics
- `g::Function` : control directions
"""
struct ControlAffineSystem
    name::String
    n::Int
    m::Int
    f::Function
    g::Function
end

"""
    ControlAffineSystem(n::Int, m::Int, f::Function, g::Function)

Construct a control affine system. Name defaults to `missing` if not provided.
"""
ControlAffineSystem(n::Int, m::Int, f::Function, g::Function) =
    ControlAffineSystem(missing, n, m, f, g)

"""
    dynamics(Σ::ControlAffineSystem, x)
    dynamics(Σ::ControlAffineSystem, x, u)

Compute the dynamics `ẋ` of a control affine system
"""
dynamics(Σ::ControlAffineSystem, x) = Σ.f(x)
dynamics(Σ::ControlAffineSystem, x, u) = Σ.f(x) + Σ.g(x) * u

"""
    simulate(Σ::ControlAffineSystem, x0, T::Float64)

Simulate a control affine system from initial condition `x0` for `T` seconds.
"""
function simulate(Σ::ControlAffineSystem, x0, T::Float64)
    function odefun(dx, x, p, t)
        return dx[1:(Σ.n)] .= dynamics(Σ, x)
    end

    return solve(ODEProblem(odefun, x0, (0, T)))
end

"""
    simulate(Σ::ControlAffineSystem, k::Function, x0, T::Float64)

Simulate a closed-loop control affine system from initial condition `x0` for `T` seconds. 
"""
function simulate(Σ::ControlAffineSystem, k::Function, x0, T::Float64)
    try
        k(Σ.n == 1 ? rand() : rand(Σ.n), 0.0)
    catch e
        if isa(e, MethodError)
            function odefun(dx, x, p, t)
                return dx[1:(Σ.n)] .= dynamics(Σ, x, k(x))
            end

            return solve(ODEProblem(odefun, x0, (0, T)))
        else
            return e
        end
    else
        function odefun(dx, x, p, t)
            return dx[1:(Σ.n)] .= dynamics(Σ, x, k(x, t))
        end

        return solve(ODEProblem(odefun, x0, (0, T)))
    end
end
