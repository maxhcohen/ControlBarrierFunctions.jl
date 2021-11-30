"""
    cbf_obstacle()
    
Construct a CBF for a circular obstacle centered at xo with radius r.

Optional argument allows one to specify the extended class K function used in the CBF
definition. Defaults to the identity, i.e., α(h(x))=h(x).
"""
function cbf_obstacle(xo, r; α=r->r)
    h(x) = (x[1] - xo[1])^2 + (x[2] - xo[2])^2 - r^2
    return ControlBarrierFunction(h, α=α)
end