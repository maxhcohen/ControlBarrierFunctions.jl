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
parameter_dim(Σ::DoubleIntegrator) = length(Σ.mass) + length(Σ.friction)

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

function regressor(Σ::DoubleIntegrator, x)
    # Pull out states
    N = Σ.N
    q̇ = N == 1 ? x[2] : x[N+1:end]

    # Pull out parameters
    m = Σ.mass

    # Regressor for uncertain parameters
    F = vcat(zeros(N,N), -(1/m)*Matrix(1.0I, N, N).*q̇)

    return F
end

function matched_regressor(Σ::DoubleIntegrator, x)
    # Pull out states
    N = Σ.N
    q̇ = N == 1 ? x[2] : x[N+1:end]

    # Matched regressor for uncertain parameters
    φ = -1.0*Matrix(1.0I, N, N).*q̇

    return φ
end

function nominal_drift(Σ::DoubleIntegrator, x)
    # Pull out states
    N = Σ.N
    q̇ = N == 1 ? x[2] : x[N+1:end]

    return vcat(q̇, zeros(N))
end

function parameters(Σ::DoubleIntegrator)
    return vcat(Σ.mass, Σ.friction)
end

function inertia_matrix(Σ::DoubleIntegrator, q)
    return Matrix(Σ.mass*I, Σ.N, Σ.N)
end

function coriolis_matrix(Σ::DoubleIntegrator, q, q̇)
    return zeros(Σ.N, Σ.N)
end

function damping_matrix(Σ::DoubleIntegrator)
    γ = Σ.friction
    if length(γ) == 1
        F = γ*Matrix(1.0I, Σ.N, Σ.N)
    else
        F = Matrix(1.0I, Σ.N, Σ.N)*γ
    end

    return F
end

function actuation_matrix(Σ::DoubleIntegrator, q)
    return Matrix(1.0I, Σ.N, Σ.N)
end