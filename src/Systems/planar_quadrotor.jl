"""
    PlanarQuadrotor <: ControlAffineSystem

ControlAffineSystem representing a quadrotor confined to a vertical plane.
"""
struct PlanarQuadrotor <: ControlAffineSystem
    mass
    inertia
    radius
end

PlanarQuadrotor() = PlanarQuadrotor(1.0, 1.0, 1.0)

state_dim(Σ::PlanarQuadrotor) = 6
control_dim(Σ::PlanarQuadrotor) = 2
degrees_of_freedom(Σ::PlanarQuadrotor) = 3

function drift(Σ::PlanarQuadrotor, x)
    return vcat(x[4:6], 0.0, -9.81, 0.0)
end

function actuation(Σ::PlanarQuadrotor, x)
    #  Pull out parameters
    m = Σ.mass
    inertia = Σ.inertia
    r = Σ.radius

    # Pull out states
    θ = x[3]

    # Columns of g
    g1 = [0.0, 0.0, 0.0, -sin(θ)/m, cos(θ)/m, r/inertia]
    g2 = [0.0, 0.0, 0.0, -sin(θ)/m, cos(θ)/m, -r/inertia]

    return hcat(g1, g2)
end