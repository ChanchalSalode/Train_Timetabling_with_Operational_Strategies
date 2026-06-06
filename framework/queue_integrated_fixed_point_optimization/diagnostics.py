# diagnostics.py

from collections import defaultdict


# =========================================
# SAFE VARIABLE VALUE
# =========================================

def var_value_safe(var, default=0.0):

    try:

        if var is None:
            return float(default)

        v = getattr(var, 'X', None)

        if v is None:
            v = getattr(var, 'x', None)

        return float(v) if v is not None else float(default)

    except Exception:

        return float(default)


# =========================================
# PERFORMANCE METRICS
# =========================================

def collect_performance_metrics_from_model(
    model,
    var_dicts,
    C,
    h_min,
    eps=1e-9
):

    try:

        tau = var_dicts.get('tau', {})
        n1 = var_dicts.get('n1', {})
        w = var_dicts.get('w', {})
        v = var_dicts.get('v', {})
        w_b = var_dicts.get('w_b', {})
        h = var_dicts.get('h', {})
        x = var_dicts.get('x', {})
        B = var_dicts.get('B', [])

    except Exception:

        return {
            'max_plf': 0.0,
            'overloaded': [],
            'cong_seconds': 0.0,
            'cong_minutes': 0.0,
            'left_behind_rate': 0.0,
            'dcr_by_ki': {},
            'dcr_details': {},
            'peak_dcr': 0.0
        }

    cap = float(C)

    # =====================================
    # Passenger Load Factor
    # =====================================

    plf_by_ki = {}

    for (k, i) in list(n1.keys()):

        tau_val = var_value_safe(
            tau.get(k)
            if isinstance(tau, dict)
            else None,
            1.0
        )

        if tau_val > 0.5:

            plf_by_ki[(k, i)] = (

                var_value_safe(
                    n1.get((k, i)),
                    0.0
                )

                /

                (cap if cap > eps else 1.0)

            )

    max_plf = (

        max(plf_by_ki.values())

        if plf_by_ki
        else 0.0

    )

    overloaded = [

        (k, i, round(p, 4))

        for (k, i), p
        in plf_by_ki.items()

        if p > 0.85

    ]

    # =====================================
    # Congestion
    # =====================================

    cong_seconds = 0.0

    services_with_cong = set(

        [
            k
            for (k, i, p)
            in overloaded
        ]

    )

    for k in services_with_cong:

        idx = (k - 1, k)

        if idx in h:

            cong_seconds += var_value_safe(
                h.get(idx)
            )

        else:

            cong_seconds += float(h_min)

    cong_minutes = cong_seconds / 60.0

    # =====================================
    # Left Behind Rate
    # =====================================

    sum_v = 0.0
    sum_w = 0.0

    for (k, i, j) in B:

        xi = var_value_safe(
            x.get((k, i))
        )

        if xi > 0.5:

            sum_v += var_value_safe(
                v.get((k, i, j))
            )

            sum_w += var_value_safe(
                w.get((k, i, j))
            )

    left_behind_rate = (

        (sum_v / sum_w)

        if sum_w > eps

        else 0.0

    )

    # =====================================
    # Demand Capacity Ratio
    # =====================================

    dcr_by_ki = {}

    dcr_details = {}

    peak_dcr = 0.0

    for (k, i) in list(w_b.keys()):

        tau_val = var_value_safe(
            tau.get(k)
            if isinstance(tau, dict)
            else None,
            1.0
        )

        arrivals = var_value_safe(
            w_b.get((k, i)),
            0.0
        )

        idx = (k - 1, k)

        if idx in h:

            headway = var_value_safe(
                h.get(idx),
                float(h_min)
            )

        else:

            headway = (

                float(h_min)

                if tau_val > 0.5

                else 0.0

            )

        if headway > eps:

            cap_per_time = cap / headway

            dcr_val = (

                arrivals / cap_per_time

                if cap_per_time > eps

                else float('inf')

            )

        else:

            dcr_val = (

                0.0

                if arrivals <= eps

                else float('inf')

            )

        dcr_by_ki[(k, i)] = float(dcr_val)

        dcr_details[(k, i)] = (

            float(dcr_val),
            float(arrivals),
            float(headway),
            float(cap)

        )

        if (

            dcr_val != float('inf')

            and

            dcr_val > peak_dcr

        ):

            peak_dcr = float(dcr_val)

    # =====================================
    # Return Metrics
    # =====================================

    return {

        'max_plf': float(max_plf),

        'overloaded': overloaded,

        'cong_seconds': float(cong_seconds),

        'cong_minutes': float(cong_minutes),

        'left_behind_rate': float(left_behind_rate),

        'dcr_by_ki': dcr_by_ki,

        'dcr_details': dcr_details,

        'peak_dcr': float(peak_dcr)

    }


# =========================================
# WRITE METRICS COMPARISON
# =========================================

def write_metrics_comparison(
    metrics_before,
    metrics_after,
    output_file
):

    with open(output_file, 'w') as cf:

        cf.write(
            'Comparison of system performance metrics\n\n'
        )

        cf.write(
            'METRIC | BEFORE_ITERATIVE | AFTER_ITERATIVE\n'
        )

        cf.write(
            '----------------------------------------------------------\n'
        )

        cf.write(

            'Peak Load Factor (max PLF) | '

            '{:.6f} | {:.6f}\n'.format(

                metrics_before.get(
                    'max_plf',
                    0.0
                ),

                metrics_after.get(
                    'max_plf',
                    0.0
                )

            )

        )

        cf.write(

            'Congestion minutes         | '

            '{:.3f} | {:.3f}\n'.format(

                metrics_before.get(
                    'cong_minutes',
                    0.0
                ),

                metrics_after.get(
                    'cong_minutes',
                    0.0
                )

            )

        )

        cf.write(

            'Left-Behind Rate           | '

            '{:.6f} | {:.6f}\n'.format(

                metrics_before.get(
                    'left_behind_rate',
                    0.0
                ),

                metrics_after.get(
                    'left_behind_rate',
                    0.0
                )

            )

        )

        cf.write(

            'Number overloaded (PLF>0.85) | '

            '{} | {}\n'.format(

                len(
                    metrics_before.get(
                        'overloaded',
                        []
                    )
                ),

                len(
                    metrics_after.get(
                        'overloaded',
                        []
                    )
                )

            )

        )

    print(
        f"[Saved] Metrics comparison -> {output_file}"
    )
# =========================================
# PREPROCESSOR DIAGNOSTICS
# =========================================

def print_preprocessing_diagnostics(
    iteration,
    station_id,
    queue_results,
    total_boarded,
    total_left
):

    print("\n----------------------------------------")

    print(
        f"ITERATION {iteration} | STATION {station_id}"
    )

    print("----------------------------------------")

    try:

        queue_len = len(queue_results)

    except Exception:

        queue_len = 0

    print(
        f"Queue records          : {queue_len}"
    )

    print(
        f"Total boarded pax      : "
        f"{float(total_boarded):.2f}"
    )

    print(
        f"Total left-behind pax  : "
        f"{float(total_left):.2f}"
    )

    print("----------------------------------------")

# =========================================
# WRITE MODEL VARIABLE VALUES
# =========================================

def write_model_vars_dump(
    model,
    var_dicts,
    basename="QFPO",
    stage="final"
):

    filename = f"{basename}_{stage}_variables.txt"

    with open(filename, "w") as f:

        f.write(
            "=====================================\n"
        )

        f.write(
            "MODEL VARIABLE VALUES\n"
        )

        f.write(
            "=====================================\n\n"
        )

        try:

            for v in model.getVars():

                try:

                    if abs(v.X) > 1e-6:

                        f.write(
                            f"{v.VarName} = {v.X}\n"
                        )

                except Exception:

                    continue

        except Exception as ex:

            f.write(
                f"\nError writing variables: {ex}\n"
            )

    print(
        f"[Saved] Variable dump -> {filename}"
    )