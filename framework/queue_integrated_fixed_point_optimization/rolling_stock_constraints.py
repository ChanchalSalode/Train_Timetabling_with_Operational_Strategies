# rolling_stock_constraints.py

from gurobipy import quicksum


def add_rolling_stock_constraints(
    mdl,
    vars,
    sets,
    params
):

    # =========================================
    # Variables
    # =========================================

    tau = vars['tau']

    z = vars['z']

    y = vars['y']

    a1 = vars['a1']

    d = vars['d']

    alpha_up1 = vars['alpha_up1']
    alpha_dn2 = vars['alpha_dn2']

    alpha_up3 = vars['alpha_up3']
    alpha_dn4 = vars['alpha_dn4']

    beta_up2 = vars['beta_up2']
    beta_dn1 = vars['beta_dn1']

    beta_up4 = vars['beta_up4']
    beta_dn3 = vars['beta_dn3']

    RS1 = vars['RS1']
    RS2 = vars['RS2']
    RS3 = vars['RS3']
    RS4 = vars['RS4']

    # =========================================
    # Sets
    # =========================================

    K_u = sets['K_u']
    K_d = sets['K_d']

    # =========================================
    # Parameters
    # =========================================

    delta_min = params['delta_min']

    M_time = params['M_time']

    RS = params['RS']

    # =========================================
    # Rolling Stock Compatibility Constraints
    # =========================================

    mdl.addConstrs(

        2 * y[k, l, 6]

        <=

        z[k, 1, 6]
        + z[l, 11, 16]
        + z[k, 3, 6]
        + z[l, 11, 14]

        for k in K_u
        for l in K_d

    )

    mdl.addConstrs(

        2 * y[k, l, 8]

        <=

        z[k, 1, 8]
        + z[l, 9, 16]
        + z[k, 3, 8]
        + z[l, 9, 14]

        for k in K_u
        for l in K_d

    )

    mdl.addConstrs(

        2 * y[k, l, 14]

        <=

        z[k, 9, 14]
        + z[l, 3, 8]
        + z[k, 11, 14]
        + z[l, 3, 6]

        for k in K_d
        for l in K_u

    )

    mdl.addConstrs(

        2 * y[k, l, 16]

        <=

        z[k, 9, 16]
        + z[l, 1, 8]
        + z[k, 11, 16]
        + z[l, 1, 6]

        for k in K_d
        for l in K_u

    )

    # =========================================
    # Turnaround Time Constraints
    # =========================================

    mdl.addConstrs(

        a1[l, 9]
        - d[k, 8]

        >=

        delta_min
        + M_time * (y[k, l, 8] - 1)

        for k in K_u
        for l in K_d

    )

    mdl.addConstrs(

        a1[l, 11]
        - d[k, 6]

        >=

        delta_min
        + M_time * (y[k, l, 6] - 1)

        for k in K_u
        for l in K_d

    )

    mdl.addConstrs(

        a1[l, 1]
        - d[k, 16]

        >=

        delta_min
        + M_time * (y[k, l, 16] - 1)

        for k in K_d
        for l in K_u

    )

    mdl.addConstrs(

        a1[l, 3]
        - d[k, 14]

        >=

        delta_min
        + M_time * (y[k, l, 14] - 1)

        for k in K_d
        for l in K_u

    )

    # =========================================
    # Incoming Rolling Stock Flow
    # =========================================

    mdl.addConstrs(

        quicksum(

            y[l, k, 16]

            for l in K_d

        )

        +

        quicksum(

            y[l, k, 14]

            for l in K_d

        )

        +

        alpha_up1[k]
        + alpha_up3[k]

        ==

        tau[k]

        for k in K_u

    )

    mdl.addConstrs(

        quicksum(

            y[l, k, 6]

            for l in K_u

        )

        +

        quicksum(

            y[l, k, 8]

            for l in K_u

        )

        +

        alpha_dn2[k]
        + alpha_dn4[k]

        ==

        tau[k]

        for k in K_d

    )

    # =========================================
    # Outgoing Rolling Stock Flow
    # =========================================

    mdl.addConstrs(

        quicksum(

            y[k, l, 6]

            for l in K_d

        )

        +

        quicksum(

            y[k, l, 8]

            for l in K_d

        )

        +

        beta_up2[k]
        + beta_up4[k]

        ==

        tau[k]

        for k in K_u

    )

    mdl.addConstrs(

        quicksum(

            y[k, l, 16]

            for l in K_u

        )

        +

        quicksum(

            y[k, l, 14]

            for l in K_u

        )

        +

        beta_dn1[k]
        + beta_dn3[k]

        ==

        tau[k]

        for k in K_d

    )

    # =========================================
    # Total Rolling Stock Limit
    # =========================================

    mdl.addConstr(

        quicksum(

            alpha_up1[k]

            for k in K_u

        )

        +

        quicksum(

            alpha_dn2[k]

            for k in K_d

        )

        +

        quicksum(

            alpha_up3[k]

            for k in K_u

        )

        +

        quicksum(

            alpha_dn4[k]

            for k in K_d

        )

        <= RS

    )

    # =========================================
    # Depot-wise Rolling Stock Limits
    # =========================================

    mdl.addConstr(

        quicksum(

            alpha_up1[k]

            for k in K_u

        )

        <= RS1

    )

    mdl.addConstr(

        quicksum(

            alpha_dn2[k]

            for k in K_d

        )

        <= RS2

    )

    mdl.addConstr(

        quicksum(

            alpha_up3[k]

            for k in K_u

        )

        <= RS3

    )

    mdl.addConstr(

        quicksum(

            alpha_dn4[k]

            for k in K_d

        )

        <= RS4

    )

    # =========================================
    # Rolling Stock Activation Constraints
    # =========================================

    mdl.addConstrs(

        alpha_up1[k]

        <=

        z[k, 1, 8]
        + z[k, 1, 6]

        for k in K_u

    )

    mdl.addConstrs(

        alpha_dn2[k]

        <=

        z[k, 9, 16]
        + z[k, 9, 14]

        for k in K_d

    )

    mdl.addConstrs(

        beta_dn1[k]

        <=

        z[k, 9, 16]
        + z[k, 11, 16]

        for k in K_d

    )

    mdl.addConstrs(

        beta_up2[k]

        <=

        z[k, 1, 8]
        + z[k, 3, 8]

        for k in K_u

    )

    mdl.addConstrs(

        alpha_up3[k]

        <=

        z[k, 3, 8]
        + z[k, 3, 6]

        for k in K_u

    )

    mdl.addConstrs(

        alpha_dn4[k]

        <=

        z[k, 11, 16]
        + z[k, 11, 14]

        for k in K_d

    )

    mdl.addConstrs(

        beta_dn3[k]

        <=

        z[k, 9, 14]
        + z[k, 11, 14]

        for k in K_d

    )

    mdl.addConstrs(

        beta_up4[k]

        <=

        z[k, 1, 6]
        + z[k, 3, 6]

        for k in K_u

    )