# sets.py

from config import *


def create_sets():

    S_d = [i for i in range(9, S_dn + 1)]
    S_dd = [i for i in range(1, S_dn + 1)]

    K_d = [i for i in range(7, K_dn + 1)]
    K_dd = [i for i in range(1, K_dn + 1)]

    K_u = [i for i in range(1, K_up + 1)]
    S_u = [i for i in range(1, S_up + 1)]

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