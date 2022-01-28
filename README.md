# CBFToolbox.jl
Control barrier functions (CBFs) and control Lyapunov functions (CLFs) written in Julia.

## Overview
This toolbox provides utilities to construct nonlinear systems and control policies based on control barrier functions and control Lyapunov functions. The utilities in this toolbox make heavy use of Julia's multiple dispatch functionality and are intended to provide a lightweight base for more complex packages that leverage CBFs and CLFs.

The utilities in this package are built into four main types:
1. `ControlAffineSystem`
2. `ControlLyapunovFunction`
3. `ControlBarrierFunction`
4. `Policy`

The `ControlAffineSystem` type represents a nonlinear control affine system and can be constructed as

    `Σ = ControlAffineSystem(n, m, f, g)`

where `n` is an integer representing the state dimension, `m` is an integer representing the control dimension, `f(x)` is a function representing the drift dynamics, and `g(x)` is a function representing the control directions.

## Questions and Contributions
If you have any questions about the toolbox, have suggestions for improvements, or would like to make your own contribution to the toolbox feel free to reach out to the repo's owner at maxcohen@bu.edu.
 
