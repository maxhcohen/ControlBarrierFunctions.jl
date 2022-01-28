# CBFToolbox.jl
Control barrier functions (CBFs) and control Lyapunov functions (CLFs) written in Julia.

## Overview
This toolbox provides utilities to construct nonlinear systems and control policies based on control barrier functions and control Lyapunov functions. The utilities in this toolbox make heavy use of Julia's multiple dispatch functionality and are intended to provide a lightweight base for more complex packages that leverage CBFs and CLFs.

The utilities in this package are built into four main types:
1. `ControlAffineSystem`
2. `ControlLyapunovFunction`
3. `ControlBarrierFunction`
4. `Policy`

The `ControlAffineSystem` type represents a nonlinear control affine system of the form

    ẋ = f(x) + g(x)u

where `x ∈ ℝⁿ` is the system state, `u ∈ ℝᵐ` is the control input, `f(x)` captures the drift dynamics,  and `g(x)` captures the control directions. A `ControlAffineSystem` can be constructed as

    Σ = ControlAffineSystem(n, m, f, g)

where `n` is an integer representing the state dimension, `m` is an integer representing the control dimension, `f(x)` is a function representing the drift dynamics, and `g(x)` is a function representing the control directions. The control input for a system can be computed by constructing a `Policy` from either a `ControlLyapunovFunction`, a `ControlBarrierFunction`, or both. A `ControlLyapunovFunction` can be constructed as 

    CLF = ControlLyapunovFunction(V, γ)

where `V(x)` is a function representing the Lyapunov function candidate and `γ(s)` is a class K function that specifies the rate of convergence of the Lyapunov function. A `ControlBarrierFunction` can be constructed in a similar fashion as

    CBF = ControlBarrierFunction(h, α)

where `h(x)` is a function defining the safe set `C  ` as `C = {x ∈ ℝⁿ | h(x) ≥ 0}` and `α` is an extended class K function that specifies how quickly the system may approach the boundary of the safe set. The `Policy` type makes heavy use of multiple dispatch. For example, the commands

    k = CLFQP(Σ, CLF)
    k = CLFQP(Σ, CLF, U)

construct a `Policy` represented as a CLF quadratic program without any control constraints and with control constraints specified by `U`, respectively. Similarly, the following commands

    k = CBFQP(Σ, CBF)
    k = CBFQP(Σ, CBF, CLF)

construct a `Policy` represented as a CBF quadratic program without and with an additional CLF constraint, respectively. Once a `Policy` is constructed, it can be evaluated by calling `u = k(x)`, which will solve the corresponding quadratic program for the control input.

## Potential changes to be made.
- Should we make control constraints a field of the `ControlAffineSystem` type rather than as an explicit input to the `Policy` constructor?
- Should we generalize how control constraints are specified using `LazySets.jl`? This would add an additional dependency to the package but could be helpful.
- Add more concrete constructions of systems and CBFs? There is some support for this currently, but could be greatly expanded.
- Should `CLFQP` and `CBFQP` types just be one `QPPolicy` type?

## Questions and Contributions
If you have any questions about the toolbox, have suggestions for improvements, or would like to make your own contribution to the toolbox feel free to reach out to the repo's owner at maxcohen@bu.edu.
 
