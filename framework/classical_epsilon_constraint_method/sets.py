# sets.py


def create_sets(params):

    # =========================================
    # Parameters
    # =========================================

    S_U = params['S_U']
    S_D = params['S_D']

    K_U = params['K_U']
    K_D = params['K_D']

    # =========================================
    # Station Sets
    # =========================================

    S_u = [
        i
        for i in range(1, S_U + 1)
    ]

    S_d = [
        i
        for i in range(9, S_D + 1)
    ]

    S_dd = [
        i
        for i in range(1, S_D + 1)
    ]

    # =========================================
    # Train Sets
    # =========================================

    K_u = [
        i
        for i in range(1, K_U + 1)
    ]

    K_d = [
        i
        for i in range(7, K_D + 1)
    ]

    K_dd = [
        i
        for i in range(1, K_D + 1)
    ]

    # =========================================
    # Multi-dimensional Sets
    # =========================================

    B = [

        (k, m, n)

        for k in K_dd
        for m in S_dd
        for n in S_dd

    ]

    A = [

        (k, l, m)

        for k in K_dd
        for l in K_dd
        for m in S_dd

    ]

    C1 = [

        (k, i)

        for k in K_dd
        for i in S_dd

    ]

    D = [

        (k, l)

        for k in K_dd
        for l in K_dd

    ]

    # =========================================
    # Return All Sets
    # =========================================

    return {

        'S_u': S_u,
        'S_d': S_d,
        'S_dd': S_dd,

        'K_u': K_u,
        'K_d': K_d,
        'K_dd': K_dd,

        'B': B,
        'A': A,
        'C1': C1,
        'D': D

    }