struct DoublePendulum <: ControlAffineSystem
    m1
    m2
    l1
    l2
    fd1
    fd2
end

DoublePendulum() = DoublePendulum(1.0, 1.0, 1.0, 1.0, 0.0, 0.0)

state_dim(::DoublePendulum) = 4
control_dim(::DoublePendulum) = 2
degrees_of_freedom(::DoublePendulum) = 2
euler_lagrange(::DoublePendulum) = true

function inertia_matrix(Σ::DoublePendulum, q)
    # Pull out states
    c2 = cos(q[2])

    # Pull out parameters
    m1 = Σ.m1
    m2 = Σ.m2
    l1 = Σ.l1
    l2 = Σ.l2

    # Entries of inertia matrix
    m11 = (m1 + m2)*l1^2 + m2*l2^2 + 2*m2*l1*l2*c2
    m12 = m2*l2^2 + m2*l1*l2*c2
    m21 = m12
    m22 = m2*l2^2

    return [m11 m12; m21 m22]
end

function coriolis_matrix(Σ::DoublePendulum, q, q̇)
    # Pull out states
    s2 = sin(q[2])

    # Pull out parameters
    m1 = Σ.m1
    m2 = Σ.m2
    l1 = Σ.l1
    l2 = Σ.l2

    # Entries of coriolis matrix
    c11 = 0.0
    c12 = -m2*l1*l2*(2*q̇[1] + q̇[2])*s2
    c21 = 0.5*m2*l1*l2*(2*q̇[1] + q̇[2])*s2
    c22 = -0.5*m2*l1*l2*q̇[1]*s2

    return [c11 c12; c21 c22]
end

function gravity_vector(Σ::DoublePendulum, q)
     # Pull out states
     s1 = sin(q[1])
     s12 = sin(q[1] + q[2])

     # Pull out parameters
     m1 = Σ.m1
     m2 = Σ.m2
     l1 = Σ.l1
     l2 = Σ.l2

     # Entries of gravity vector
     g1 = (m1 + m2)*l1*s1 + m2*l2*s12
     g2 = m2*l2*s12

     return -9.81*[g1, g2]
end

function damping_matrix(Σ::DoublePendulum)
    # Pull out parameters
    fd1 = Σ.fd1
    fd2 = Σ.fd2

    return [fd1 0.0; 0.0 fd2]
end

function actuation_matrix(Σ::DoublePendulum, q)
    return Matrix(1.0I, 2, 2)
end

function drift(Σ::DoublePendulum, x)
    # Pull out states
    q = x[1:2]
    q̇ = x[3:4]

    # Dynamics
    M = inertia_matrix(Σ, q)
    C = coriolis_matrix(Σ, q, q̇)
    F = damping_matrix(Σ)
    τg = gravity_vector(Σ, q)

    return vcat(q̇, -M\(C*q̇ + F*q̇ - τg))
end

function actuation(Σ::DoublePendulum, x)
    # Pull out states
    q = x[1:2]

    # Dynamics
    M = inertia_matrix(Σ, q)
    B = actuation_matrix(Σ, q)

    return vcat(zeros(2,2), M/B)
end