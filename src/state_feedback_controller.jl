"""
    StateFeedbackController <: FeedbackController

Custom state feedback control policy of the form u = k(x), where k is a function of x
representing the custom control law.
"""
struct StateFeedbackController <: FeedbackController
    control_law::Function
end

(k::StateFeedbackController)(x) = k.control_law(x)

