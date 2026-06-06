# pwl_constraints.py

from gurobipy import GRB


def add_pwl_constraints(
    mdl,
    vars,
    sets,
    params,
    pwl_up,
    pwl_down,
    backlog
):

    # =========================================
    # Variables
    # =========================================

    d = vars['d']

    w = vars['w']
    v = vars['v']

    t_up = vars['t_up']
    t_down = vars['t_down']

    A_up = vars['A_up']
    A_dn = vars['A_dn']

    # =========================================
    # Sets
    # =========================================

    K_u = sets['K_u']
    K_d = sets['K_d']

    S_u = sets['S_u']
    S_d = sets['S_d']

    # =========================================
    # Parameters
    # =========================================

    H = params['H']

    T_WINDOW = params['T_WINDOW']

    # =========================================
    # Time Transformation Variables
    # =========================================

    mdl.addConstrs(t_up[k, i] == d[k, i] - H for k in K_u for i in S_u)
    mdl.addConstrs(t_down[k, i] == d[k, i] - H for k in K_d for i in S_d)

    for k in K_u:
        for i in S_u:
            for j in S_u:
                if i < j and (i, j) in pwl_up:
                    x_breaks, y_breaks = pwl_up[(i, j)]
                    A_ki = mdl.addVar(lb=0.0, name=f"A_up_{k}_{i}_{j}")
                    mdl.addGenConstrPWL(t_up[k, i], A_ki, x_breaks, y_breaks, name=f"PWL_up_{k}_{i}_{j}")
                    A_up[(k, i, j)] = A_ki
                    if k == 1:
                        base = backlog.get((i, j), 0.0)
                        mdl.addConstr(w[k, i, j] == base + A_ki)
                    else:
                        A_prev = mdl.addVar(lb=0.0, name=f"A_up_prev_{k}_{i}_{j}")
                        mdl.addGenConstrPWL(t_up[k - 1, i], A_prev, x_breaks, y_breaks, name=f"PWL_up_prev_{k}_{i}_{j}")
                        mdl.addConstr(w[k, i, j] == v[k - 1, i, j] + (A_ki - A_prev))
    for k in K_d:
        for i in S_d:
            for j in S_d:
                if i < j and (i, j) in pwl_down:
                    x_breaks, y_breaks = pwl_down[(i, j)]
                    A_ki = mdl.addVar(lb=0.0, name=f"A_dn_{k}_{i}_{j}")
                    mdl.addGenConstrPWL(t_down[k, i], A_ki, x_breaks, y_breaks, name=f"PWL_dn_{k}_{i}_{j}")
                    A_dn[(k, i, j)] = A_ki
                    if k == 7:
                        base = backlog.get((i, j), 0.0)
                        mdl.addConstr(w[k, i, j] == base + A_ki)
                    else:
                        A_prev = mdl.addVar(lb=0.0, name=f"A_dn_prev_{k}_{i}_{j}")
                        mdl.addGenConstrPWL(t_down[k - 1, i], A_prev, x_breaks, y_breaks, name=f"PWL_dn_prev_{k}_{i}_{j}")
                        mdl.addConstr(w[k, i, j] == v[k - 1, i, j] + (A_ki - A_prev))