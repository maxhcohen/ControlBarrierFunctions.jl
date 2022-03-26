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