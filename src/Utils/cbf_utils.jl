"""
	Obstacle

Abstract type representing obstacle.
"""
abstract type Obstacle end

"""
	CircularObstacle <: Obstacle

Circular obstacle with center c and radius r
"""
struct CircularObstacle <: Obstacle
	c::Vector{Float64}
	r::Float64
end

"""
	ControlBarrierFunction(O::CircularObstacle, α)

Construct a CBF for a ciruclar obstacle.
"""
function ControlBarrierFunction(O::CircularObstacle, α)
	h(x) = (x[1] - O.c[1])^2 + (x[2] - O.c[2])^2 - O.r^2
	return ControlBarrierFunction(h, α)
end