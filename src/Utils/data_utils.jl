"""
	vec2mat(x::Vector{Vector{Float64}})
	vec2mat(x::Vector{Float64})

Convert vector of vectors into matrix. If input is a single vector do nothing.
"""
vec2mat(x::Vector{Vector{Float64}}) = reduce(hcat, x)
vec2mat(x::Vector{Float64}) = x

"""
    vec2mat!(Σ::ControlAffineSystem)

In-place modification of system trajectory
"""
function vec2mat!(Σ::ControlAffineSystem)
    Σ.xs = vec2mat(Σ.xs)

    return Σ
end
