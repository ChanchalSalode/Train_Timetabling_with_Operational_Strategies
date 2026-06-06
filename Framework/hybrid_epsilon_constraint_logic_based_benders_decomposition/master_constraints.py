# master_constraints.py

from gurobipy import quicksum


def add_master_constraints(
    mdl,
    vars,
    sets,
    data,
    params
):

    # =========================================
    # Variables
    # =========================================

    tau = vars['tau']

    z = vars['z']
    y = vars['y']

    x = vars['x']

    a1 = vars['a1']
    d = vars['d']

    h = vars['h']

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

    S_u = sets['S_u']
    S_d = sets['S_d']

    # =========================================
    # Parameters
    # =========================================

    M = params['M']

    M_time = params['M_time']

    h_min = params['h_min']
    h_max = params['h_max']

    delta_min = params['delta_min']

    RS = params['RS']

    # =========================================
    # Data
    # =========================================

    pr = data['pr']

    ac = data['ac']
    dc = data['dc']

    e = data['e']

    # =========================================
    # Operation Planning Constraints
    # =========================================

    mdl.addConstrs(

        z[k, 1, 6]
        + z[k, 1, 8]
        + z[k, 3, 6]
        + z[k, 3, 8]

        == tau[k]

        for k in K_u

    )

    mdl.addConstrs(

        z[k, 9, 14]
        + z[k, 9, 16]
        + z[k, 11, 14]
        + z[k, 11, 16]

        == tau[k]

        for k in K_d

    )

    # =========================================
    # x and tau relationship
    # =========================================

    mdl.addConstrs(

        x[k, i] <= tau[k]

        for i in S_u
        for k in K_u

    )

    mdl.addConstrs(

        x[k, i] <= tau[k]

        for i in S_d
        for k in K_d

    )

    # =========================================
    # Skip-stop logic
    # =========================================

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_u

            if i > 6
        )

        <= M * (1 - z[k, 1, 6])

        for k in K_u

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_u

            if i < 3
        )

        <= M * (1 - z[k, 3, 6])

        for k in K_u

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_u

            if i > 6
        )

        <= M * (1 - z[k, 3, 6])

        for k in K_u

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_u

            if i < 3
        )

        <= M * (1 - z[k, 3, 8])

        for k in K_u

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_d

            if i > 14
        )

        <= M * (1 - z[k, 9, 14])

        for k in K_d

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_d

            if i < 11
        )

        <= M * (1 - z[k, 11, 14])

        for k in K_d

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_d

            if i > 14
        )

        <= M * (1 - z[k, 11, 14])

        for k in K_d

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_d

            if i < 11
        )

        <= M * (1 - z[k, 11, 16])

        for k in K_d

    )

    # =========================================
    # Consecutive service constraints
    # =========================================

    mdl.addConstrs(

        x[k - 1, i]
        + x[k, i]

        >= tau[k]

        for k in K_u
        for i in S_u

        if k != 1

    )

    mdl.addConstrs(

        x[k - 1, i]
        + x[k, i]

        >= tau[k]

        for k in K_d
        for i in S_d

        if k != 10

    )

    # =========================================
    # Initial departure time
    # =========================================

    mdl.addConstr(
        d[1, 1] == 27000
    )

    mdl.addConstr(
        d[10, 9] == 27000
    )

    # =========================================
    # Time bounds
    # =========================================

    mdl.addConstrs(

        d[k, 1] <= 28800

        for k in K_u

    )

    mdl.addConstrs(

        d[k, 3] <= 28800

        for k in K_u

    )

    mdl.addConstrs(

        d[k, 9] <= 28800

        for k in K_d

    )

    mdl.addConstrs(

        d[k, 11] <= 28800

        for k in K_d

    )

    # =========================================
    # Running time constraints
    # =========================================

    mdl.addConstrs(

        a1[k, i]
        - d[k, i - 1]

        ==

        pr[i - 1, i]

        + x[k, i - 1] * ac[i - 1, i]

        + x[k, i] * dc[i - 1, i]

        for k in K_u
        for i in S_u

        if i != 1

    )

    mdl.addConstrs(

        a1[k, i]
        - d[k, i - 1]

        ==

        pr[i - 1, i]

        + x[k, i - 1] * ac[i - 1, i]

        + x[k, i] * dc[i - 1, i]

        for k in K_d
        for i in S_d

        if i != 9

    )

    # =========================================
    # Dwelling time constraints
    # =========================================

    mdl.addConstrs(

        d[k, i]
        - a1[k, i]

        >= e[i] * x[k, i]

        for k in K_u
        for i in S_u

    )

    mdl.addConstrs(

        d[k, i]
        - a1[k, i]

        >= e[i] * x[k, i]

        for k in K_d
        for i in S_d

    )

    # =========================================
    # Minimum served stations
    # =========================================

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_u
        )

        >= 4

        for k in K_u

    )

    mdl.addConstrs(

        quicksum(
            x[k, i]

            for i in S_d
        )

        >= 4

        for k in K_d

    )

    # =========================================
    # Headway constraints
    # =========================================

    mdl.addConstrs(

        h_min * tau[k]

        <= h[k - 1, k]

        for k in K_u

        if k != 1

    )

    mdl.addConstrs(

        h[k - 1, k]

        <= h_max * tau[k]

        for k in K_u

        if k != 1

    )

    mdl.addConstrs(

        h_min * tau[k]

        <= h[k - 1, k]

        for k in K_d

        if k != 10

    )

    mdl.addConstrs(

        h[k - 1, k]

        <= h_max * tau[k]

        for k in K_d

        if k != 10

    )

    # =========================================
    # Timetable propagation
    # =========================================

    mdl.addConstrs(

        d[k, i]

        ==

        d[k - 1, i]
        + h[k - 1, k]

        for k in K_u
        for i in S_u

        if k != 1

    )

    mdl.addConstrs(

        d[k, i]

        ==

        d[k - 1, i]
        + h[k - 1, k]

        for k in K_d
        for i in S_d

        if k != 10

    )

    # =========================================
    # Turnback feasibility constraints
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
    # Transfer time constraints
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
    # Rolling stock flow balance
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

        == tau[k]

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

        == tau[k]

        for k in K_d

    )

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

        == tau[k]

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

        == tau[k]

        for k in K_d

    )

    # =========================================
    # Rolling stock availability
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
    # Alpha/Beta linking constraints
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