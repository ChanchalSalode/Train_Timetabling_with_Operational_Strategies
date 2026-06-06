# objective.py

from gurobipy import GRB, quicksum


def set_objective(
    mdl,
    vars,
    sets
):

    # =========================================
    # Variables
    # =========================================

    zz = vars['zz']

    d = vars['d']

    # =========================================
    # Sets
    # =========================================

    K_u = sets['K_u']
    K_d = sets['K_d']

    S_u = sets['S_u']
    S_d = sets['S_d']

    # =========================================
    # Objective Function
    # =========================================

    objective = (

        quicksum(
            zz[k,1,6]
            for k in K_u
        )

        +

        quicksum(
            zz[k,1,8]
            for k in K_u
        )

        +

        quicksum(
            zz[k,3,6]
            for k in K_u
        )

        +

        quicksum(
            zz[k,3,8]
            for k in K_u
        )

        +

        quicksum(
            zz[k,9,14]
            for k in K_d
        )

        +

        quicksum(
            zz[k,9,16]
            for k in K_d
        )

        +

        quicksum(
            zz[k,11,14]
            for k in K_d
        )

        +

        quicksum(
            zz[k,11,16]
            for k in K_d
        )

        +

        quicksum(
            d[k,i] - d[k-1,i]
            for k in K_u
            for i in S_u
            if k != 1
        )

        +

        quicksum(
            d[k,i] - d[k-1,i]
            for k in K_d
            for i in S_d
            if k != 7
        )

    )

    # =========================================
    # Set Objective
    # =========================================

    mdl.setObjective(
        objective,
        GRB.MINIMIZE
    )