"""
    SingleIntegrator <: ControlAffineSystem

ControlAffineSystem representing a single integrator with N states.
"""
struct SingleIntegrator <: ControlAffineSystem
    N
end

state_dim(Σ::SingleIntegrator) = Σ.N
control_dim(Σ::SingleIntegrator) = Σ.N
drift(Σ::SingleIntegrator, x) = Σ.N == 1 ? 0.0 : zeros(Σ.N)
actuation(Σ::SingleIntegrator, x) = Σ.N == 1 ? 1.0 : Matrix(1.0I, Σ.N, Σ.N)