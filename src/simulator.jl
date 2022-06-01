struct Simulator
    x0
    t0
    dt
    tf
    ts
end

Simulator(x0) = Simulator(x0, 0.0, 0.01, 10.0, 0.0:0.01:10.0)
Simulator(x0, tf) = Simulator(x0, 0.0, 0.01, tf, 0.0:0.01:tf)
Simulator(x0, dt, tf) = Simulator(x0, 0.0, dt, tf, 0.0:dt:tf)

function (sim::Simulator)(Σ::ControlAffineSystem)
    trajectory = simulate(Σ, sim.x0, [sim.t0, sim.tf])

    return trajectory
end

Base.length(sim::Simulator) = length(sim.ts)

function (sim::Simulator)(Σ::ControlAffineSystem, k::FeedbackController)
    rhs(x, p, t) = _f(Σ, x) + _g(Σ, x)*k(x)
    prob = ODEProblem(rhs, sim.x0, [sim.t0, sim.tf])
    sol = solve(prob, Tsit5())

    return sol
end

function (sim::Simulator)(Σ::ControlAffineSystem, k::Union{CLFController, CBFController, CBFCLFController, ISSCLFController, ISSCBFController, ISSCBFCLFController, TISSCBFController})
    rhs(x, p, t) = _f(Σ, x) + _g(Σ, x)*k(Σ, x)
    prob = ODEProblem(rhs, sim.x0, [sim.t0, sim.tf])
    sol = solve(prob, Tsit5())

    return sol
end

function (sim::Simulator)(Σ::ControlAffineSystem, k::Union{CBFController, ISSCBFController, TISSCBFController}, k0::Union{CLFController, ISSCLFController})
    rhs(x, p, t) = _f(Σ, x) + _g(Σ, x)*k(Σ, x, k0(Σ, x))
    prob = ODEProblem(rhs, sim.x0, [sim.t0, sim.tf])
    sol = solve(prob, Tsit5())

    return sol
end

function (sim::Simulator)(Σ::ControlAffineSystem, k0::Union{CLFController, ISSCLFController}, k::Union{CBFController, ISSCBFController, TISSCBFController})
    rhs(x, p, t) = _f(Σ, x) + _g(Σ, x)*k(Σ, x, k0(Σ, x))
    prob = ODEProblem(rhs, sim.x0, [sim.t0, sim.tf])
    sol = solve(prob, Tsit5())

    return sol
end

function (sim::Simulator)(Σ::ControlAffineSystem, k::FLController, y::ConfigurationError)
    rhs(x, p, t) = _f(Σ, x) + _g(Σ, x)*k(Σ, y, x)
    prob = ODEProblem(rhs, sim.x0, [sim.t0, sim.tf])
    sol = solve(prob, Tsit5())

    return sol
end

function (sim::Simulator)(Σ::ControlAffineSystem, k::CBFController, k0::FLController, y::ConfigurationError)
    rhs(x, p, t) = _f(Σ, x) + _g(Σ, x)*k(Σ, x, k0(Σ, y, x))
    prob = ODEProblem(rhs, sim.x0, [sim.t0, sim.tf])
    sol = solve(prob, Tsit5())

    return sol
end