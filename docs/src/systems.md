# System Definitions

## Background
The starting point for most the utilities in this package is the `ControlAffineSystem` type. This type represents a nonlinear control affine system described by:
```math
\dot{x} = f(x) + g(x)u,
```
where ``x\in\mathbb{R}^n`` is the system state, ``u\in\mathbb{R}^n`` is the control input, ``f\,:\,\mathbb{R}^n\rightarrow\mathbb{R}^n`` is a vector field characterizing the drift dynamics, and ``g\,:\,\mathbb{R}^n\rightarrow\mathbb{R}^{n\times m}`` is a matrix whose columns capture the control directions. A control affine system can be constructed with:
```julia
Σ = ControlAffineSystem(name, n, m, f, g)
```
where `name::String` is a string describing the system, `n::Int` is an integer denoting the state dimension, `m::Int` is an integer denoting the control dimension, and `f::Function` and `g::Function` are functions of the form `f(x)` and `g(x)` which take as input an n-dimensional vector `x` or a scalar if `n=1`. Alternatively, if you do not want to assign a name to system you can use:
```julia
Σ = ControlAffineSystem(n, m, f, g)
```
to construct a system without a name. 

## Implementation
```@docs
ControlAffineSystem
ControlAffineSystem(n::Int, m::Int, f::Function, g::Function)
```