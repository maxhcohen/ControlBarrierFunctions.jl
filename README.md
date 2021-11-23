# CBFToolbox.jl
Control barrier functions in Julia

## Purpose
The purpose of this repository is to provide a unified framework for designing CBF-based controllers in Julia

## Laundry list of things to do
 - Discuss relationship with python CBF toolbox (what worked well and what needed improvement)
 - Decide on best method to construct and solve QPs. Currently using Convex.jl as a modeling framework and ECOS solver. Another popular modeling framework is JuMP.jl, should we consider that?
 - How to organize data structures? Currently the main ones are ControlAffineSystem, ControlBarrierFunction, and ControlLyapunovFunction.
 
