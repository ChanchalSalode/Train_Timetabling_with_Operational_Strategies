# linearization_constraints.py


def add_linearization_constraints(
    mdl,
    vars,
    sets,
    params
):

    # =========================================
    # Variables
    # =========================================

    zz = vars['zz']

    z = vars['z']

    d = vars['d']

    # =========================================
    # Sets
    # =========================================

    K_u = sets['K_u']
    K_d = sets['K_d']

    # =========================================
    # Parameters
    # =========================================

    M = params['M']

    # =========================================
    # UP DIRECTION
    # =========================================

    mdl.addConstrs(

        zz[k, 1, 6]

        >=

        (
            d[k, 6]
            - d[k, 1]
        )

        -

        M * (

            1
            - z[k, 1, 6]

        )

        for k in K_u

    )

    mdl.addConstrs(

        zz[k, 1, 8]

        >=

        (
            d[k, 8]
            - d[k, 1]
        )

        -

        M * (

            1
            - z[k, 1, 8]

        )

        for k in K_u

    )

    mdl.addConstrs(

        zz[k, 3, 6]

        >=

        (
            d[k, 6]
            - d[k, 3]
        )

        -

        M * (

            1
            - z[k, 3, 6]

        )

        for k in K_u

    )

    mdl.addConstrs(

        zz[k, 3, 8]

        >=

        (
            d[k, 8]
            - d[k, 3]
        )

        -

        M * (

            1
            - z[k, 3, 8]

        )

        for k in K_u

    )

    # =========================================
    # DOWN DIRECTION
    # =========================================

    mdl.addConstrs(

        zz[k, 9, 14]

        >=

        (
            d[k, 14]
            - d[k, 9]
        )

        -

        M * (

            1
            - z[k, 9, 14]

        )

        for k in K_d

    )

    mdl.addConstrs(

        zz[k, 9, 16]

        >=

        (
            d[k, 16]
            - d[k, 9]
        )

        -

        M * (

            1
            - z[k, 9, 16]

        )

        for k in K_d

    )

    mdl.addConstrs(

        zz[k, 11, 14]

        >=

        (
            d[k, 14]
            - d[k, 11]
        )

        -

        M * (

            1
            - z[k, 11, 14]

        )

        for k in K_d

    )

    mdl.addConstrs(

        zz[k, 11, 16]

        >=

        (
            d[k, 16]
            - d[k, 11]
        )

        -

        M * (

            1
            - z[k, 11, 16]

        )

        for k in K_d

    )

    # =========================================
    # UPPER BOUNDS
    # =========================================

    mdl.addConstrs(

        zz[k, 1, 6]

        <=

        M * z[k, 1, 6]

        for k in K_u

    )

    mdl.addConstrs(

        zz[k, 1, 8]

        <=

        M * z[k, 1, 8]

        for k in K_u

    )

    mdl.addConstrs(

        zz[k, 3, 6]

        <=

        M * z[k, 3, 6]

        for k in K_u

    )

    mdl.addConstrs(

        zz[k, 3, 8]

        <=

        M * z[k, 3, 8]

        for k in K_u

    )

    mdl.addConstrs(

        zz[k, 9, 14]

        <=

        M * z[k, 9, 14]

        for k in K_d

    )

    mdl.addConstrs(

        zz[k, 9, 16]

        <=

        M * z[k, 9, 16]

        for k in K_d

    )

    mdl.addConstrs(

        zz[k, 11, 14]

        <=

        M * z[k, 11, 14]

        for k in K_d

    )

    mdl.addConstrs(

        zz[k, 11, 16]

        <=

        M * z[k, 11, 16]

        for k in K_d

    )