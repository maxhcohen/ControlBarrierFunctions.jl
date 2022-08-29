"""
    plot_phase_portrait(X0::Vector{Vector{Float64}}, f::Function, T::Float64; lw=1)
    plot_phase_portrait(X0::Vector{Float64}, f::Function, T::Float64; lw=1)
    plot_phase_portrait(xs::Vector{Float64}, ys::Vector{Float64}, f::Function, T::Float64; lw=1)
    plot_phase_portrait(xs::StepRangeLen, ys::StepRangeLen, f::Function, T::Float64; lw=1)

Plot the phase portrait of a dynamical system defined by `ẋ = f(x)` over `t ∈ [0,T]`.

The initial conditions for the trajectories can be specified as a list of `x` and `y`
coordinates `xs` and `ys`, respectively, or by directly passing in a list of initial
conditions `X0`.
"""
function plot_phase_portrait(X0::Vector{Vector{Float64}}, f::Function, T::Float64; lw=1)
    # Simulate vector field from time 0 to T from each initial condition
    X = [simulate(f, x0, T) for x0 in X0]

    # Plot each trajectory on the (x, y) plane
    fig = plot()
    for x in X
        plot!(x, idxs=(1,2), label="", lw=lw)
    end

    return fig
end

function plot_phase_portrait(X0::Vector{Float64}, f::Function, T::Float64; lw=1)
    return plot_phase_portrait([X0], f, T, lw=lw)
end

function plot_phase_portrait(xs::Vector{Float64}, ys::Vector{Float64}, f::Function, T::Float64; lw=1)
    # Make meshgrid from (xs, ys)
    Xs, Ys = meshgrid(xs, ys)

    # Construct vector of initial conditions from meshgrid
    X0 = [[x, y] for (x, y) in zip(Xs, Ys)]

    return plot_phase_portrait(X0, f, T, lw=lw)
end

function plot_phase_portrait(xs::Union{StepRangeLen, LinRange}, ys::Union{StepRangeLen, LinRange}, f::Function, T::Float64; lw=1)
    return plot_phase_portrait(collect(xs), collect(ys), f, T, lw=lw)
end

"""
    plot_phase_portrait!(X0::Vector{Vector{Float64}}, f::Function, T::Float64; lw=1)
    plot_phase_portrait!(X0::Vector{Float64}, f::Function, T::Float64; lw=1)
    plot_phase_portrait!(xs::Vector{Float64}, ys::Vector{Float64}, f::Function, T::Float64; lw=1)
    plot_phase_portrait!(xs::StepRangeLen, ys::StepRangeLen, f::Function, T::Float64; lw=1)

Plot the phase portrait of a dynamical system defined by `ẋ = f(x)` over `t ∈ [0,T]`.

The initial conditions for the trajectories can be specified as a list of `x` and `y`
coordinates `xs` and `ys`, respectively, or by directly passing in a list of initial
conditions `X0`.
"""
function plot_phase_portrait!(X0::Vector{Vector{Float64}}, f::Function, T::Float64; lw=1)
    # Simulate vector field from time 0 to T from each initial condition
    X = [simulate(f, x0, T) for x0 in X0]

    # Get colors for each trajectory: TODO add this as a kwarg
    # colors = colormap("RdBu", length(X))

    # Plot each trajectory on the (x, y) plane
    fig = Plots.current()
    for (i, x) in enumerate(X)
        plot!(x, idxs=(1,2), label="", lw=lw)
    end

    return fig
end

function plot_phase_portrait!(X0::Vector{Float64}, f::Function, T::Float64; lw=1)
    return plot_phase_portrait!([X0], f, T, lw=lw)
end

function plot_phase_portrait!(xs::Vector{Float64}, ys::Vector{Float64}, f::Function, T::Float64; lw=1)
    # Make meshgrid from (xs, ys)
    Xs, Ys = meshgrid(xs, ys)

    # Construct vector of initial conditions from meshgrid
    X0 = [[x, y] for (x, y) in zip(Xs, Ys)]

    return plot_phase_portrait!(X0, f, T, lw=lw)
end

function plot_phase_portrait!(xs::Union{StepRangeLen, LinRange}, ys::Union{StepRangeLen, LinRange}, f::Function, T::Float64; lw=1)
    return plot_phase_portrait!(collect(xs), collect(ys), f, T, lw=lw)
end

"""
    simulate(f::Function, x, T::Float64)

Simulate the differential equation `ẋ = f(x)` from initial condition `x` over `[0,T]`
"""
function simulate(f::Function, x, T::Float64)
    rhs(x, p, t) = f(x[1], x[2])
    problem = ODEProblem(rhs, x, [0.0, T])

    return solve(problem, Tsit5())
end