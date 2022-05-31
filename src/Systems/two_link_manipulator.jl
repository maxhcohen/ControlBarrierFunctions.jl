struct TwoLinkManipulator <: ControlAffineSystem
    p1
    p2
    p3
    fd1
    fd2
end

state_dim(Σ::TwoLinkManipulator) = 4
control_dim(Σ::TwoLinkManipulator) = 2
degrees_of_freedom(Σ::TwoLinkManipulator) = 2
euler_lagrange(Σ::TwoLinkManipulator) = true

function inertia_matrix(Σ::TwoLinkManipulator, q)
    # Pull out states
    c2 = cos(q[2])

    # Pull out parameters
    p1 = Σ.p1
    p2 = Σ.p2
    p3 = Σ.p3

    # Entries of inertia matrix
    m11 = p1 + 2*p3*c2
    m12 = p2 + p3*c2
    m21 = m12
    m22 = p2

    return [m11 m12; m21 m22]
end

function coriolis_matrix(Σ::TwoLinkManipulator, q, q̇)
    # Pull out states
    s2 = sin(q[2])

    # Pull out parameters
    p3 = Σ.p3

    # Entries of coriolis matrix
    c11 = -p3*s2*q̇[2]
    c12 = -p3*s2*(q̇[1] + q̇[2])
    c21 = p3*s2*q̇[1]
    c22 = 0.0

    return [c11 c12; c21 c22]
end

function damping_matrix(Σ::TwoLinkManipulator)
    # Pull out parameters
    fd1 = Σ.fd1
    fd2 = Σ.fd2

    return [fd1 0.0; 0.0 fd2]
end

function actuation_matrix(Σ::TwoLinkManipulator, q)
    return Matrix(1.0I, 2, 2)
end

function drift(Σ::TwoLinkManipulator, x)
    # Pull out states
    q = x[1:2]
    q̇ = x[3:4]

    # Dynamics
    M = inertia_matrix(Σ, q)
    C = coriolis_matrix(Σ, q, q̇)
    F = damping_matrix(Σ)

    return vcat(q̇, -M\(C*q̇ + F*q̇))
end

function actuation(Σ::TwoLinkManipulator, x)
    # Pull out states
    q = x[1:2]

    # Dynamics
    M = inertia_matrix(Σ, q)
    B = actuation_matrix(Σ, q)

    return vcat(zeros(2,2), M/B)
end
