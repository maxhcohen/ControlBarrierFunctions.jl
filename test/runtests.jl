using Test
using CBFToolbox
using LinearAlgebra
using ForwardDiff

### Test functionality of ControlAffineSystem
n = 2
m = 2
f(x) = zeros(2)
g(x) = [1.0 0.0; 0.0 1.0]
Σ = ControlAffineSystem("single integrator", n, m, f, g)

# Make sure ControlAffineSystem was constructed correctly
@test Σ.name == "single integrator"
@test Σ.n == n
@test Σ.m == m
@test Σ.f(zeros(2)) == f(zeros(2))
@test Σ.g(zeros(2)) == g(zeros(2))

# Make sure we can simulate a simple system
x0 = [-1.0, 1.0]
T = 20.0
sol = simulate(Σ, x0, T)
@test sol(T) == x0

# Make sure we can simulate with a controller
k(x) = -x
sol = simulate(Σ, k, x0, T)
@test norm(sol(T)) ≤ 0.01

# Repeat with a time-varying controller
k(x, t) = -x
sol = simulate(Σ, k, x0, T)
@test norm(sol(T)) ≤ 0.01

### Test functionality of ControlBarrierFunction
xo = [-1.0, 1.0]
ro = 0.4
h(x) = norm(x - xo)^2 - ro^2
cbf = ControlBarrierFunction(h, Σ)

# Make sure CBF was constructed correctly
@test cbf(x0) == h(x0)
@test cbf.∇h(x0) == ForwardDiff.gradient(h, x0)
@test cbf.Lfh(x0) == ForwardDiff.gradient(h, x0)' * f(x0)
@test cbf.Lgh(x0) == ForwardDiff.gradient(h, x0)' * g(x0)
@test cbf.α(cbf(x0)) == cbf(x0)

### Test functionality for ExplicitSafetyFilter
kS = ExplicitSafetyFilter(cbf, Σ, k)

# Make sure CBF was positive along solution
x0 = [-2.1, 2.0]
sol = simulate(Σ, kS, x0, T)
ts = 0.0:0.01:T
@test minimum(h.(sol.(ts))) ≥ 0.0

### Test functionality for QPSafetyFilter
kd(x) = -x

# Test standard QP
kQP = QPSafetyFilter(cbf, Σ, kd)
sol = simulate(Σ, kQP, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

# Test with time-varying controller
kQP = QPSafetyFilter(cbf, Σ, (x, t) -> kd(x))
sol = simulate(Σ, kQP, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

### Test functionality for TunableQPSafetyFilter

# Test standard QP
ktQP = TunableQPSafetyFilter(cbf, Σ, kd)
sol = simulate(Σ, ktQP, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

# Test with time-varying controller
ktQP = TunableQPSafetyFilter(cbf, Σ, (x, t) -> kd(x))
sol = simulate(Σ, ktQP, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

### Test functionality for SmoothSafetyFilter
kHalfSontag = SmoothSafetyFilter(cbf, Σ, kd)
sol = simulate(Σ, kHalfSontag, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSontag = SmoothSafetyFilter(cbf, Σ, kd; formula="sontag")
sol = simulate(Σ, kSontag, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSoftplus = SmoothSafetyFilter(cbf, Σ, kd; formula="softplus")
sol = simulate(Σ, kSoftplus, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

# Make sure time-varying controllers work as well
kHalfSontag = SmoothSafetyFilter(cbf, Σ, (x, t) -> kd(x))
sol = simulate(Σ, kHalfSontag, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSontag = SmoothSafetyFilter(cbf, Σ, (x, t) -> kd(x); formula="sontag")
sol = simulate(Σ, kSontag, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSoftplus = SmoothSafetyFilter(cbf, Σ, (x, t) -> kd(x); formula="softplus")
sol = simulate(Σ, kSoftplus, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

### Test functionality for ISSfSmoothSafetyFilter
ε = 1.0
kHalfSontagISSf = ISSfSmoothSafetyFilter(cbf, Σ, kd, ε)
sol = simulate(Σ, kHalfSontagISSf, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSontagISSf = ISSfSmoothSafetyFilter(cbf, Σ, kd; formula="sontag", ε)
sol = simulate(Σ, kSontagISSf, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSoftplusISSf = ISSfSmoothSafetyFilter(cbf, Σ, kd; formula="softplus", ε)
sol = simulate(Σ, kSoftplusISSf, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

# Make sure time-varying controllers work as well
kHalfSontagISSf = ISSfSmoothSafetyFilter(cbf, Σ, (x, t) -> kd(x), ε)
sol = simulate(Σ, kHalfSontagISSf, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSontagISSf = ISSfSmoothSafetyFilter(cbf, Σ, (x, t) -> kd(x); formula="sontag", ε)
sol = simulate(Σ, kSontagISSf, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0

kSoftplusISSf = ISSfSmoothSafetyFilter(cbf, Σ, (x, t) -> kd(x); formula="softplus", ε)
sol = simulate(Σ, kSoftplusISSf, x0, T)
@test minimum(h.(sol.(ts))) ≥ 0.0
