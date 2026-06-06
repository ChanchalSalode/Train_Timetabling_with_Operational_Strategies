# objective.py

from gurobipy import quicksum, GRB


def set_objective(mdl, vars, sets):

    y = vars['y']

    K_u = sets['K_u']
    K_d = sets['K_d']

    mdl.setObjective(

        quicksum(
            y[l, k, 6]
            for l in K_u
            for k in K_d
        )

        +

        quicksum(
            y[l, k, 8]
            for l in K_u
            for k in K_d
        )

        +

        quicksum(
            y[l, k, 14]
            for l in K_d
            for k in K_u
        )

        +

        quicksum(
            y[l, k, 16]
            for l in K_d
            for k in K_u
        ),

        GRB.MAXIMIZE
    )