# performance_metrics.py


def compute_performance_metrics(
    vars,
    sets,
    params
):

    # =========================================
    # Variables
    # =========================================

    tau = vars['tau']
    h = vars['h']

    x = vars['x']

    w = vars['w']
    w_b = vars['w_b']

    v = vars['v']

    n1 = vars['n1']

    # =========================================
    # Sets
    # =========================================

    B = sets['B']

    # =========================================
    # Parameters
    # =========================================

    C = params['C']
    h_min = params['h_min']

    # =========================================
    # Constants
    # =========================================

    eps = 1e-6

    cap = C

    # =========================================
    # 1. Peak Load Factor (PLF)
    # =========================================

    plf_by_ki = {}

    for (k, i) in n1.keys():

        if k in tau.keys():

            if tau[k].X > 0.5:

                val = n1[k, i].X

                plf = (
                    val / float(cap)
                    if cap > eps
                    else float('inf')
                )

                plf_by_ki[(k, i)] = plf

    max_plf = (
        max(plf_by_ki.values())
        if plf_by_ki
        else 0.0
    )

    # =========================================
    # Overloaded Services
    # =========================================

    overloaded = [

        (k, i, p)

        for (k, i), p
        in plf_by_ki.items()

        if p > 0.85

    ]

    # =========================================
    # 2. Congestion Duration
    # =========================================

    cong_seconds = 0.0

    services_with_cong = set(

        [k for (k, i, p) in overloaded]

    )

    for k in services_with_cong:

        idx = (k - 1, k)

        if idx in h.keys():

            cong_seconds += h[idx].X

        else:

            cong_seconds += 120.0

    cong_minutes = cong_seconds / 60.0

    # =========================================
    # 3. Left-Behind Passenger Rate
    # =========================================

    sum_v = 0.0

    sum_w = 0.0

    valid_B = [
        (k, i, j)
        for (k, i, j) in B
        if i < j
    ]

    for (k, i, j) in valid_B:

        if ((k, i) in x.keys()) and ((k, j) in x.keys()):

            if x[k, i].X > 0.5 and x[k, j].X > 0.5:

                sum_v += float(v[k, i, j].X)

                sum_w += float(w[k, i, j].X)

    left_behind_rate = (

        sum_v / sum_w

        if sum_w > eps

        else 0.0

    )

    # =========================================
    # 4. Demand-to-Capacity Ratio (DCR)
    # =========================================

    dcr_by_ki = {}

    hw_by_k = {}

    for (k, i) in w_b.keys():

        if k in tau.keys():

            if tau[k].X > 0.5:

                arrivals = float(w_b[k, i].X)

                idx = (k - 1, k)

                if idx in h.keys():

                    hw = max(
                        h[idx].X,
                        eps
                    )

                else:

                    hw = h_min

                hw_by_k[k] = hw

                cap_per_time = (
                    float(cap) / hw
                    if hw > eps
                    else float('inf')
                )

                dcr_val = (
                    arrivals / cap_per_time
                    if cap_per_time > eps
                    else float('inf')
                )

                dcr_by_ki[(k, i)] = dcr_val

    peak_dcr = (

        max(dcr_by_ki.values())

        if dcr_by_ki

        else 0.0

    )
    print('\nPerformance metrics:')

    print(
        'Peak Load Factor (PLF) = max_k,i n1[k,i]/C = {:.4f}'.format(
            max_plf
        )
    )

    if overloaded:

        print(
            'Overloaded (PLF>0.85) (k,i,PLF):'
        )

        for (k, i, p) in overloaded:

            print(
                '  {}: {:.4f}'.format(
                    (k, i),
                    p
                )
            )

    else:

        print(
            'No overloaded (PLF>0.85) stations/services.'
        )

    print(
        'Congestion duration (seconds): {:.1f}  | minutes: {:.2f}'.format(
            cong_seconds,
            cong_minutes
        )
    )

    print(
        'Left-Behind Passenger Rate (sum v / sum w): {:.4f}'.format(
            left_behind_rate
        )
    )

    print(
        'Peak Demand-to-Capacity Ratio (DCR_peak): {:.4f}'.format(
            peak_dcr
        )
    )

    # =========================================
    # Return Metrics
    # =========================================

    return {

        'max_plf': max_plf,

        'overloaded': overloaded,

        'cong_seconds': cong_seconds,

        'cong_minutes': cong_minutes,

        'left_behind_rate': left_behind_rate,

        'peak_dcr': peak_dcr

    }