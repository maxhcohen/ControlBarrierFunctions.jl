"""
	FeedbackPolicy

Struct representing a feedback control policy.

# Fields
- `control`: function u = control(x) describing the control law.
"""
struct FeedbackPolicy <: Policy
    control
end

"""
	(k::FeedbackPolicy)(x)

Evaluate FeedbackPolicy at state x.
"""
function (k::FeedbackPolicy)(x)
    return k.control(x)
end

"""
	TimeVaryingFeedbackPolicy

Struct representing a non-stationary feedback control policy.

# Fields
- `control`: function u = control(x, t) describing the control law.
"""
struct TimeVaryingFeedbackPolicy <: Policy
    control
end

"""
	(k::TimeVaryingFeedbackPolicy)(x)

Evaluate TimeVaryingFeedbackPolicy at state x and time t.
"""
function (k::TimeVaryingFeedbackPolicy)(x, t)
    return k.control(x, t)
end
