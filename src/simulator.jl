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
    X = zeros(state_dim(Σ), length(sim.ts))
    X[:,1] = sim.x0
    for i in 1:length(sim)-1
        x = X[:, i]
        u = k(x)
        X[:,i+1] = integrate(Σ, x, u, [sim.ts[i], sim.ts[i+1]])
    end

    return X
end

function (sim::Simulator)(Σ::ControlAffineSystem, k::CLFController)
    X = zeros(state_dim(Σ), length(sim.ts))
    X[:,1] = sim.x0
    for i in 1:length(sim)-1
        x = X[:, i]
        u = k(Σ, x)
        X[:,i+1] = integrate(Σ, x, u, [sim.ts[i], sim.ts[i+1]])
    end

    return X
end

function (sim::Simulator)(Σ::ControlAffineSystem, k::CBFController)
    X = zeros(state_dim(Σ), length(sim.ts))
    X[:,1] = sim.x0
    for i in 1:length(sim)-1
        x = X[:, i]
        u = k(Σ, x)
        X[:,i+1] = integrate(Σ, x, u, [sim.ts[i], sim.ts[i+1]])
    end

    return X
end

function (sim::Simulator)(Σ::ControlAffineSystem, k::CBFController, k0::CLFController)
    X = zeros(state_dim(Σ), length(sim.ts))
    X[:,1] = sim.x0
    for i in 1:length(sim)-1
        x = X[:, i]
        ud = k0(Σ, x)
        u = k(Σ, x, ud)
        X[:,i+1] = integrate(Σ, x, u, [sim.ts[i], sim.ts[i+1]])
    end

    return X
end

function (sim::Simulator)(Σ::ControlAffineSystem, k0::CLFController, k::CBFController)
    X = zeros(state_dim(Σ), length(sim.ts))
    X[:,1] = sim.x0
    for i in 1:length(sim)-1
        x = X[:, i]
        ud = k0(Σ, x)
        u = k(Σ, x, ud)
        X[:,i+1] = integrate(Σ, x, u, [sim.ts[i], sim.ts[i+1]])
    end

    return X
end