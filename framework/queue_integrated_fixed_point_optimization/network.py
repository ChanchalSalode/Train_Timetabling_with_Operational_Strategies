# network.py

from config import (
    S_U,
    S_D,
    K_U,
    K_D
)


def create_network_sets():
    """Create station, train, and index sets for the QFPO model."""

    S_u = [i for i in range(1, S_U + 1)]
    S_d = [i for i in range(S_U + 1, S_D + 1)]
    S_dd = [i for i in range(1, S_D + 1)]

    K_u = [i for i in range(1, K_U + 1)]
    K_d = [i for i in range(K_U + 1, K_D + 1)]
    K_dd = [i for i in range(1, K_D + 1)]

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


def create_parameter_dict():
    """Expose scalar model parameters as an explicit dependency bundle."""

    from config import (
        h_min,
        h_max,
        delta_min,
        RS,
        C,
        M_time,
        M_flow,
        M,
        H,
        T_WINDOW,
        BIN,
        _times,
        _n_bins,
        eps
    )

    return {
        'h_min': h_min,
        'h_max': h_max,
        'delta_min': delta_min,
        'RS': RS,
        'C': C,
        'M_time': M_time,
        'M_flow': M_flow,
        'M': M,
        'H': H,
        'T_WINDOW': T_WINDOW,
        'BIN': BIN,
        'times': _times,
        'n_bins': _n_bins,
        'eps': eps
    }
