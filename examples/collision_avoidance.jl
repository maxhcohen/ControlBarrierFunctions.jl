### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ fc0c7450-cae2-4e3c-a0ec-d985605b5e9d
begin
    import Pkg
    Pkg.activate(Base.current_project())

    using CBFToolbox, Plots, LaTeXStrings
end

# ╔═╡ a025cade-52f5-11ec-0ee0-53eeaa34a5e8
md"## Control barrier function tutorial

This notebook contains a quick tutorial on using control barrier functions to develop safety-critical controllers for dynamical systems. Formally, given a nonlinear control affine system of the form

``\dot{x} = f(x) + g(x)u,``

the objective is to find a control policy ``u=k(x)`` such that the system trajectory ``x(t)`` remains in a given safe set ``\mathcal{C}`` for all time.

### Step 1: Activate current environment
To get started, we first have to activate the current environment so that we can access CBFToolbox.jl
"

# ╔═╡ ed61e480-c8cd-4d70-9349-dd90b4b88459
md"### Step 2: Define system dynamics

Next we need to define the dynamics of our system. This is done by constructing a `ControlAffineSystem` by specifying four different properties: 
1. an integer `n::Int` representing the state dimension
2. an integer `m::Int` representing the control dimension
3. a function `f(x)` representing the drift dynamics
4. a function `g(x)` representing the actuation dynamics

For this example, we consider a simple two-dimensional single integrator $\dot{x}=u$ with the following properties:
"

# ╔═╡ 090088ee-78d4-45f1-b0eb-7db13518311d
n = 2

# ╔═╡ 49a748c4-1fbb-457f-a85b-b0f6add48daa
m = 2

# ╔═╡ 45ca320e-d66a-4ad3-b5f7-307ba76196af
f(x) = zeros(2)

# ╔═╡ c3ee9765-830b-4200-9b78-8f764ac6b759
g(x) = [1.0 0.0; 0.0 1.0]

# ╔═╡ 34fa1efa-9376-4ef2-a902-f56c46e5c0a7
Σ = ControlAffineSystem(n, m, f, g)

# ╔═╡ 4dfdb669-e4cf-4e1c-b831-4606879b0a39
md"### Step 3: Construct CBF
Next, we need to construct a CBF $h(x)$ that describes our safe set as

``\mathcal{C}=\{x\in\mathbb{R}^n\,|\,h(x)\geq0\}``.

We do this by constructing `ControlBarrierFunction` directly from a function `h(x)` that defines the safe set. For this simple example, we want to avoid a circular obstacle of radius $r\in\mathbb{R}_{>0}$ centered at $y\in\mathbb{R}^2$. The safe set associated with this obstacle can be defined as

``h(x)=(x_1 - y_1)^2 + (x_2 - y_2)^2 - r^2,``

where $(x_1,x_2)$ denotes the position of the integrator. We can then construct the CBF as follows:
"

# ╔═╡ 9dbd8773-fa7a-475a-b839-4f104d4d9f7c
y = [-1.0, 1.0]

# ╔═╡ 1cd40dea-7513-472d-843c-9475418a1058
r = 0.4

# ╔═╡ b6398139-0160-44e7-94ce-f9c146a5b5e5
h(x) = (x[1] - y[1])^2 + (x[2] - y[2])^2 - r^2

# ╔═╡ adc1ed79-3ef5-43c3-a34d-c44ab248ce5e
md"The last thing we need to do before constructing the CBF is to specify the *extended class* $\mathcal{K}$ *function* $\alpha(h(x))$ associated with the CBF that determines the rate at which the system is allowed to approach the boundary of the safe set. Some examples of these functions are

$\alpha(h(x))=\gamma h(x)$

$\alpha(h(x))=\gamma h^3(x)$

for some positive constant $\gamma$.
"

# ╔═╡ 35165558-fce1-41cb-b7e1-4a1ab2194030
γ = 1.0

# ╔═╡ 9a5d2cce-1ef9-4b62-a23c-11a12d6fadce
linear_classK(h) = γ*h

# ╔═╡ 335636dc-3f40-47a5-a904-c9cb481b7ccc
cubic_classK(h) = γ*h^3

# ╔═╡ d8dbe60e-d059-49bd-8945-eb55b0ff7024
md"Finally, we can construct our CBF by calling

`CBF = ControlBarrierFunction(h, α=classK)`

where `classK` is the extended class $\mathcal{K}$ function of choice. Note that `α` is a `kwarg` - it defaults to $\alpha(h(x))=h(x)$ if left unspecified.
"

# ╔═╡ 8776a1b2-dca6-4813-8255-6ea61a0e16f5
CBF = ControlBarrierFunction(h, α=linear_classK)

# ╔═╡ e148f29b-dd9f-4dfe-bb46-ca8b84301b05
md"### Step 5: Construct nominal feedback control policy
One of the primary advantages of CBFs is their ability to modify any nominal control policy $k(x)$ in a minimally invasive fashion so as to guarantee safety. Specifically, given a CBF $h(x)$ and a nominal control policy $k(x)$ one can solve the following optimization problem to find input closest to the nominal input that also keeps the system safe

``
\begin{equation}
\begin{aligned}
\min_{u\in\mathcal{U}} && \tfrac{1}{2}\|u - k(x)\|^2 \\ 
\text{subject to} && L_fh(x) + L_gh(x)u\geq -\alpha(h(x)).
\end{aligned}
\end{equation}
``

The above optimization problem is a quadratic program (QP) as the objective function is quadratic in the decision variable $u$ and the constraints are linear in $u$, and can thus be solved efficiently in real-time. For this simple example, we'll take the nominal policy as the simple state feedback control law $k(x)=-\kappa x$, where $\kappa$ is a control gain, which attempts to drive the system to the origin.
"

# ╔═╡ 6f8d8dac-5cf8-4607-99ea-0c7dd3fbc56a
κ = 1.0

# ╔═╡ b13803aa-3147-4e07-a2d6-2b0bc2b77b74
k(x) = -κ*x

# ╔═╡ 59fde591-2486-4a8c-840d-91f81c42bc14
md"### Step 6: Run a simulation
Given our system, CBF, and nominal control policy we can now run a simulation to see how the system responds to the given control policy. The only additional information we have to specify is
 - the initial time of the simulation `t0`
 - the final time of the simulation `tf`
 - the simulation step size `dt` used by the ODE solver for integrating the closed-loop system
 - the initial condition of our system `x0`

With all this in place we can simulate the system by using the `run_sim` function as

`t, x = run_sim(t0, tf, dt, x0, k, Σ, CBF)`

which returns

 - a vector representing the simulation time steps `t`
 - the resulting trajectory of the system `x`

"

# ╔═╡ 3df77cf5-6f46-4499-983f-36306a4776e7
t0 = 0.0

# ╔═╡ 3320c450-cf97-4d04-bfff-5fb528366f3b
tf = 10.0

# ╔═╡ d64fa132-0df2-45bf-beac-ce3c08ff53e6
dt = 0.01

# ╔═╡ ee675b63-da4d-4fa6-acc2-49573a8e4b58
x0 = [-2.1, 2.0]

# ╔═╡ 6e6d6dee-30bf-4339-b52f-873744359462
t, x = run_sim(t0, tf, dt, x0, k, Σ, CBF)

# ╔═╡ 7b9bab0b-c1b5-4ade-bf3a-23959f21846b
md"### Step 7: Plot the results"

# ╔═╡ 59f39a6c-f83b-4740-b1ce-110c9549cab7
custom_plots()

# ╔═╡ ea4ce2da-de68-4b5d-9bd5-def3e1ebacc0
mycolors = custom_colors()

# ╔═╡ aa0d16c8-5ffa-4b2e-9aa3-9ac37a7a0e09
begin
	plot(xlabel=L"x_1", ylabel=L"x_2")
	plot!(x[1,:], x[2,:], label=L"x(t)")
	plot!(circle_shape(y[1], y[2], r), seriestype=[:shape], fillcolor=mycolors[2], 
			fillalpha=0.5, linecolor=:black, edgecolor=:black, label="Obstacle")
end

# ╔═╡ f3366d63-25d9-4c72-a78c-3b49aec5906c
begin
	plot(xlabel=L"t", ylabel=L"x(t)")
	plot!(t, x', label=[L"x_1" L"x_2"])
end

# ╔═╡ 75409a90-ef0d-490b-8faf-89cbf0a88699
begin
	plot(xlabel=L"t", ylabel=L"h(x(t))")
	plot!(t, [h(x[:,i]) for i in 1:length(t)])
	hline!([0.0], ls=:dot, c=:black, label=L"h(x)=0")
end

# ╔═╡ Cell order:
# ╟─a025cade-52f5-11ec-0ee0-53eeaa34a5e8
# ╠═fc0c7450-cae2-4e3c-a0ec-d985605b5e9d
# ╟─ed61e480-c8cd-4d70-9349-dd90b4b88459
# ╠═090088ee-78d4-45f1-b0eb-7db13518311d
# ╠═49a748c4-1fbb-457f-a85b-b0f6add48daa
# ╠═45ca320e-d66a-4ad3-b5f7-307ba76196af
# ╠═c3ee9765-830b-4200-9b78-8f764ac6b759
# ╠═34fa1efa-9376-4ef2-a902-f56c46e5c0a7
# ╟─4dfdb669-e4cf-4e1c-b831-4606879b0a39
# ╠═9dbd8773-fa7a-475a-b839-4f104d4d9f7c
# ╠═1cd40dea-7513-472d-843c-9475418a1058
# ╠═b6398139-0160-44e7-94ce-f9c146a5b5e5
# ╟─adc1ed79-3ef5-43c3-a34d-c44ab248ce5e
# ╠═35165558-fce1-41cb-b7e1-4a1ab2194030
# ╠═9a5d2cce-1ef9-4b62-a23c-11a12d6fadce
# ╠═335636dc-3f40-47a5-a904-c9cb481b7ccc
# ╟─d8dbe60e-d059-49bd-8945-eb55b0ff7024
# ╠═8776a1b2-dca6-4813-8255-6ea61a0e16f5
# ╟─e148f29b-dd9f-4dfe-bb46-ca8b84301b05
# ╠═6f8d8dac-5cf8-4607-99ea-0c7dd3fbc56a
# ╠═b13803aa-3147-4e07-a2d6-2b0bc2b77b74
# ╟─59fde591-2486-4a8c-840d-91f81c42bc14
# ╠═3df77cf5-6f46-4499-983f-36306a4776e7
# ╠═3320c450-cf97-4d04-bfff-5fb528366f3b
# ╠═d64fa132-0df2-45bf-beac-ce3c08ff53e6
# ╠═ee675b63-da4d-4fa6-acc2-49573a8e4b58
# ╠═6e6d6dee-30bf-4339-b52f-873744359462
# ╟─7b9bab0b-c1b5-4ade-bf3a-23959f21846b
# ╠═59f39a6c-f83b-4740-b1ce-110c9549cab7
# ╠═ea4ce2da-de68-4b5d-9bd5-def3e1ebacc0
# ╠═aa0d16c8-5ffa-4b2e-9aa3-9ac37a7a0e09
# ╠═f3366d63-25d9-4c72-a78c-3b49aec5906c
# ╠═75409a90-ef0d-490b-8faf-89cbf0a88699
