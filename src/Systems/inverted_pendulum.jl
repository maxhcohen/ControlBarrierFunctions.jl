struct InvertedPendulum <: ControlAffineSystem
    mass
    len
    friction
    gravity
end

InvertedPendulum() = InvertedPendulum(1.0, 1.0, 1.0, 9.81)
InvertedPendulum(mass, len, friction) = InvertedPendulum(mass, len, friction, 9.81)

state_dim(Σ::InvertedPendulum) = 2
control_dim(Σ::InvertedPendulum) = 1
degrees_of_freedom(Σ::InvertedPendulum) = 1

function drift(Σ::InvertedPendulum, x)
    # Pull out states
    q = x[1]
    q̇ = x[2]

    # Pull out parameters
    m = Σ.mass
    l = Σ.len
    γ = Σ.friction
    g = Σ.gravity

    # Drift dynamics
    f = vcat(q̇, -(γ/m)*q̇ + (g/l)*sin(q))

    return f
end

function actuation(Σ::InvertedPendulum, x)
    # Pull out parameters
    m = Σ.mass
    l = Σ.len

    # Control directions
    g = [0.0, 1/(m*l^2)]

    return g
end