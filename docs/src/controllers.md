# Controllers

## Background

## Quadratic Programming Safety Filters
```@docs
ExplicitSafetyFilter
ExplicitSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)
ExplicitSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem)
QPSafetyFilter
QPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)
QPSafetyFilter(cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function)
TunableQPSafetyFilter
```

## Smooth Safety Filters
```@docs
SmoothSafetyFilter
ISSfSmoothSafetyFilter
```
