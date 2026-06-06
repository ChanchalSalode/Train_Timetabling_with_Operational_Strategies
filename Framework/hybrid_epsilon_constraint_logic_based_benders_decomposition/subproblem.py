# subproblem.py

from gurobipy import (
    Model,
    GRB,
    quicksum
)


def solve_passenger_subproblem(
    sets,
    params,
    data,
    x_val,
    d_val,
    z_val
):
    

    # =========================================
    # Sets
    # =========================================

    B = sets['B']
    C1 = sets['C1']

    K_u = sets['K_u']
    K_d = sets['K_d']

    S_u = sets['S_u']
    S_d = sets['S_d']

    # =========================================
    # Parameters
    # =========================================

    C = params['C']

    M = params['M']

    # =========================================
    # Data
    # =========================================

    p = data['p']

    # =========================================
    # Subproblem Model
    # =========================================

    sub = Model(
        "Passenger_Subproblem"
    )

    sub.setParam(
        'OutputFlag',
        1
    )

    # =========================================
    # Variables
    # =========================================

    w = sub.addVars(
        B,
        lb=0,
        name='w'
    )

    w_b = sub.addVars(
        C1,
        lb=0,
        name='w_b'
    )

    w_b1 = sub.addVars(
        B,
        lb=0,
        name='w_b1'
    )

    n_b = sub.addVars(
        C1,
        lb=0,
        name='n_b'
    )

    n_b1 = sub.addVars(
        B,
        lb=0,
        name='n_b1'
    )

    n_a = sub.addVars(
        C1,
        lb=0,
        name='n_a'
    )

    n1 = sub.addVars(
        C1,
        lb=0,
        ub=C,
        name='n1'
    )

    v = sub.addVars(
        B,
        lb=0,
        name='v'
    )

    sai = sub.addVars(
        C1,
        vtype=GRB.BINARY,
        name='sai'
    )

    # =========================================
    # Passenger Demand Constraints
    # =========================================

    sub.addConstrs(

        w[1, i, j]

        ==

        (p[i, j] / 1800) * 120

        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        w[k, i, j]

        ==

        v[k - 1, i, j]

        +

        (p[i, j] / 1800)

        *

        (
            d_val[(k, i)]
            - d_val[(k - 1, i)]
        )

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

        if k != 1

    )

    sub.addConstrs(

        w[10, i, j]

        ==

        (p[i, j] / 1800) * 120

        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        w[k, i, j]

        ==

        v[k - 1, i, j]

        +

        (p[i, j] / 1800)

        *

        (
            d_val[(k, i)]
            - d_val[(k - 1, i)]
        )

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

        if k != 10

    )

    # =========================================
    # Boarding / Alighting
    # =========================================

    sub.addConstrs(

        w_b[k, i]

        ==

        quicksum(

            w_b1[k, i, j]

            for j in S_u

            if i < j

        )

        for k in K_u
        for i in S_u

    )

    sub.addConstrs(

        w_b[k, i]

        ==

        quicksum(

            w_b1[k, i, j]

            for j in S_d

            if i < j

        )

        for k in K_d
        for i in S_d

    )

    sub.addConstrs(

        n_a[k, i]

        ==

        quicksum(

            n_b1[k, j, i]

            for j in S_u

            if j < i

        )

        for k in K_u
        for i in S_u

    )

    sub.addConstrs(

        n_a[k, i]

        ==

        quicksum(

            n_b1[k, j, i]

            for j in S_d

            if j < i

        )

        for k in K_d
        for i in S_d

    )

    # =========================================
    # Passenger Flow Constraints
    # =========================================

    sub.addConstrs(

        n_b[k, i]
        - n_b1[k, i, j]

        ==

        w_b[k, i]
        - w_b1[k, i, j]

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        n_b[k, i]
        - n_b1[k, i, j]

        ==

        w_b[k, i]
        - w_b1[k, i, j]

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    # =========================================
    # Linearization Constraints
    # =========================================

    sub.addConstrs(

        w_b1[k, i, j] >= 0

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        <=

        w[k, i, j]

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        <=

        M * x_val[(k, i)]

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        <=

        M * x_val[(k, j)]

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        >=

        w[k, i, j]

        -

        M * (

            2
            - x_val[(k, i)]
            - x_val[(k, j)]

        )

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    # =========================================
    # Down Direction
    # =========================================

    sub.addConstrs(

        w_b1[k, i, j] >= 0

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        <=

        w[k, i, j]

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        <=

        M * x_val[(k, i)]

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        <=

        M * x_val[(k, j)]

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        w_b1[k, i, j]

        >=

        w[k, i, j]

        -

        M * (

            2
            - x_val[(k, i)]
            - x_val[(k, j)]

        )

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    # =========================================
    # Capacity Constraints
    # =========================================

    sub.addConstrs(

        n_b[k, i]

        <=

        w_b[k, i]

        for k in K_u
        for i in S_u

    )

    sub.addConstrs(

        n_b[k, i]

        <=

        C
        - n1[k, i - 1]
        + n_a[k, i]

        for k in K_u
        for i in S_u

        if i != 1

    )

    sub.addConstrs(

        n_b[k, i]

        >=

        w_b[k, i]
        - M * (1 - sai[k, i])

        for k in K_u
        for i in S_u

    )

    sub.addConstrs(

        n_b[k, i]

        >=

        C
        - n1[k, i - 1]
        + n_a[k, i]
        - M * sai[k, i]

        for k in K_u
        for i in S_u

        if i != 1

    )

    # =========================================
    # Down Direction Capacity
    # =========================================

    sub.addConstrs(

        n_b[k, i]

        <=

        w_b[k, i]

        for k in K_d
        for i in S_d

    )

    sub.addConstrs(

        n_b[k, i]

        <=

        C
        - n1[k, i - 1]
        + n_a[k, i]

        for k in K_d
        for i in S_d

        if i != 9

    )

    sub.addConstrs(

        n_b[k, i]

        >=

        w_b[k, i]
        - M * (1 - sai[k, i])

        for k in K_d
        for i in S_d

    )

    sub.addConstrs(

        n_b[k, i]

        >=

        C
        - n1[k, i - 1]
        + n_a[k, i]
        - M * sai[k, i]

        for k in K_d
        for i in S_d

        if i != 9

    )

    # =========================================
    # Terminal Constraints
    # =========================================

    sub.addConstrs(

        n1[k, 6]

        <=

        C * (1 - z_val[k, 1, 6])

        for k in K_u

    )

    sub.addConstrs(

        n1[k, 6]

        <=

        C * (1 - z_val[k, 3, 6])

        for k in K_u

    )

    sub.addConstrs(

        n1[k, 8]

        <=

        C * (1 - z_val[k, 1, 8])

        for k in K_u

    )

    sub.addConstrs(

        n1[k, 8]

        <=

        C * (1 - z_val[k, 3, 8])

        for k in K_u

    )

    sub.addConstrs(

        n1[k, 14]

        <=

        C * (1 - z_val[k, 9, 14])

        for k in K_d

    )

    sub.addConstrs(

        n1[k, 14]

        <=

        C * (1 - z_val[k, 11, 14])

        for k in K_d

    )

    sub.addConstrs(

        n1[k, 16]

        <=

        C * (1 - z_val[k, 9, 16])

        for k in K_d

    )

    sub.addConstrs(

        n1[k, 16]

        <=

        C * (1 - z_val[k, 11, 16])

        for k in K_d

    )

    # =========================================
    # sai Constraints
    # =========================================

    sub.addConstrs(

        sai[k, i]

        >=

        (
            (
                C
                - n1[k, i - 1]
                + n_a[k, i]
            )

            -

            w_b[k, i]

        ) / M

        for k in K_u
        for i in S_u

        if i != 1

    )

    sub.addConstrs(

        sai[k, i]

        <=

        1

        +

        (
            (
                C
                - n1[k, i - 1]
                + n_a[k, i]
            )

            -

            w_b[k, i]

        ) / M

        for k in K_u
        for i in S_u

        if i != 1

    )

    sub.addConstrs(

        sai[k, i]

        >=

        (
            (
                C
                - n1[k, i - 1]
                + n_a[k, i]
            )

            -

            w_b[k, i]

        ) / M

        for k in K_d
        for i in S_d

        if i != 9

    )

    sub.addConstrs(

        sai[k, i]

        <=

        1

        +

        (
            (
                C
                - n1[k, i - 1]
                + n_a[k, i]
            )

            -

            w_b[k, i]

        ) / M

        for k in K_d
        for i in S_d

        if i != 9

    )

    # =========================================
    # Passenger Balance
    # =========================================

    sub.addConstrs(

        n1[k, 1]

        ==

        n_b[k, 1]

        for k in K_u

    )

    sub.addConstrs(

        n1[k, 9]

        ==

        n_b[k, 9]

        for k in K_d

    )

    sub.addConstrs(

        n1[k, i]

        ==

        n1[k, i - 1]
        - n_a[k, i]
        + n_b[k, i]

        for k in K_u
        for i in S_u

        if i != 1

    )

    sub.addConstrs(

        n1[k, i]

        ==

        n1[k, i - 1]
        - n_a[k, i]
        + n_b[k, i]

        for k in K_d
        for i in S_d

        if i != 9

    )

    # =========================================
    # Left Behind Passengers
    # =========================================

    sub.addConstrs(

        v[k, i, j]

        <=

        w[k, i, j]
        - n_b1[k, i, j]

        +

        M * (

            2
            - x_val[(k, i)]
            - x_val[(k, j)]

        )

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        v[k, i, j]

        >=

        w[k, i, j]
        - n_b1[k, i, j]

        -

        M * (

            2
            - x_val[(k, i)]
            - x_val[(k, j)]

        )

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        v[k, i, j]

        <=

        M * x_val[(k, i)]

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    sub.addConstrs(

        v[k, i, j]

        <=

        M * x_val[(k, j)]

        for k in K_u
        for i in S_u
        for j in S_u

        if i < j

    )

    # =========================================
    # Down Direction Left Behind
    # =========================================

    sub.addConstrs(

        v[k, i, j]

        <=

        w[k, i, j]
        - n_b1[k, i, j]

        +

        M * (

            2
            - x_val[(k, i)]
            - x_val[(k, j)]

        )

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        v[k, i, j]

        >=

        w[k, i, j]
        - n_b1[k, i, j]

        -

        M * (

            2
            - x_val[(k, i)]
            - x_val[(k, j)]

        )

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        v[k, i, j]

        <=

        M * x_val[(k, i)]

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    sub.addConstrs(

        v[k, i, j]

        <=

        M * x_val[(k, j)]

        for k in K_d
        for i in S_d
        for j in S_d

        if i < j

    )

    # =========================================
    # Dummy Objective
    # =========================================

    sub.setObjective(
        0,
        GRB.MINIMIZE
    )

    # =========================================
    # Optimize
    # =========================================

    sub.optimize()

    # =========================================
    # Infeasible
    # =========================================

    if sub.Status != GRB.OPTIMAL:

        return None

    # =========================================
    # Return Results
    # =========================================

    return {

        "w": {

            (k, i, j): w[k, i, j].X

            for (k, i, j) in B

        },

        "w_b": {

            (k, i): w_b[k, i].X

            for (k, i) in C1

        },

        "n1": {

            (k, i): n1[k, i].X

            for (k, i) in C1

        },

        "v": {

            (k, i, j): v[k, i, j].X

            for (k, i, j) in B

        },

        "n_b1": {

            (k, i, j): n_b1[k, i, j].X

            for (k, i, j) in B

        }

    }