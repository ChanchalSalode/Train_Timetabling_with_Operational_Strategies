# objective.py

from gurobipy import (
    GRB,
    quicksum
)


# =========================================
# SET OBJECTIVE FUNCTION
# =========================================

def set_objective(
    mdl,
    vars
):

    # =====================================
    # Variables
    # =====================================

    y = vars['y']

    K_u = vars['K_u']
    K_d = vars['K_d']

    # =====================================
    # Objective:
    # Maximize Successful Turnbacks
    # =====================================

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
