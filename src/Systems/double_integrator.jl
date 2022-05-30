struct DoubleIntegrator <: ControlAffineSystem
    N
    mass
    friction
end

DoubleIntegrator(N::Int) = DoubleIntegrator(N, 1.0, 1.0)

state_dim(Σ::DoubleIntegrator) = 2*Σ.N
control_dim(Σ::DoubleIntegrator) = Σ.N
degrees_of_freedom(Σ::DoubleIntegrator) = Σ.N
euler_lagrange(Σ::DoubleIntegrator) = true

function drift(Σ::DoubleIntegrator, x)
    # Pull out states
    q̇ = Σ.N == 1 ? x[2] : x[Σ.N+1:end]

    # Pull out parameters
    m = Σ.mass
    γ = Σ.friction

    # Drift dynamics
    f = vcat(q̇, -(γ/m)*q̇)

    return f
end

function actuation(Σ::DoubleIntegrator, x)
    # Pull out parameters
    m = Σ.mass

    # Control directions
    g = Σ.N == 1 ? [0.0, 1/m] : vcat(zeros(Σ.N, Σ.N), Matrix((1/m)I, Σ.N, Σ.N))

    return g
end

function inertia_matrix(Σ::DoubleIntegrator, q)
    return Matrix(Σ.mass*I, Σ.N, Σ.N)
end

function coriolis_matrix(Σ::DoubleIntegrator, q, q̇)
    return zeros(Σ.N, Σ.N)
end

function viscous_friction(Σ::DoubleIntegrator, q̇)
    γ = Σ.friction
    if length(γ) == 1
        F = γ*q̇
    else
        F = [γ[i]*q̇[i] for i in 1:length(γ)]
    end

    return F
end

function actuation_matrix(Σ::DoubleIntegrator, q)
    return Matrix(1.0I, Σ.N, Σ.N)
end