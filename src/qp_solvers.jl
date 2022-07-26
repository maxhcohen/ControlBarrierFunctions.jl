"""
    qp_model(solver::String)

Instantiate a JuMP QP model using a specific solver.

Solver options:
 - OSQP
"""
function qp_model(solver::String)
    if solver == "OSQP"
        model = Model(OSQP.Optimizer)
    elseif solver == "Ipopt"
        model = Model(Ipopt.Optimizer)
    elseif solver == "Gurobi"
        model = Model(Gurobi.Optimizer)
    end

    set_silent(model)
    return model
end