# metrics.py

from helpers import (
    var_value_safe
)

from config import (
    C,
    h_min,
    eps
)


# =========================================
# COLLECT PERFORMANCE METRICS
# =========================================

def collect_performance_metrics_from_model(

    model,

    var_dicts

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

    # =====================================
    # Passenger Load Factor (PLF)
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

                (C if C > eps else 1.0)

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
    # Congestion Time
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
    # Left-Behind Passenger Rate
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

        sum_v / sum_w

        if sum_w > eps

        else 0.0

    )

    # =====================================
    # Demand-Capacity Ratio (DCR)
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

        # ================================
        # DCR Calculation
        # ================================

        if headway > eps:

            cap_per_time = C / headway

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

            float(C)

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
# PRINT METRICS SUMMARY
# =========================================

def print_metrics_summary(

    metrics

):

    print("\n=================================")

    print("PERFORMANCE METRICS SUMMARY")

    print("=================================")

    print(

        f"Peak Load Factor : "

        f"{metrics.get('max_plf', 0.0):.4f}"

    )

    print(

        f"Congestion Minutes : "

        f"{metrics.get('cong_minutes', 0.0):.2f}"

    )

    print(

        f"Left-Behind Rate : "

        f"{metrics.get('left_behind_rate', 0.0):.4f}"

    )

    print(

        f"Peak DCR : "

        f"{metrics.get('peak_dcr', 0.0):.4f}"

    )

    print(

        f"Overloaded Sections : "

        f"{len(metrics.get('overloaded', []))}"

    )


# =========================================
# EXPORT METRICS COMPARISON
# =========================================

def export_metrics_comparison(

    metrics_before,

    metrics_after,

    output_file

):

    with open(

        output_file,

        'w'

    ) as cf:

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

                metrics_before.get('max_plf', 0.0),

                metrics_after.get('max_plf', 0.0)

            )

        )

        cf.write(

            'Congestion minutes | '

            '{:.3f} | {:.3f}\n'.format(

                metrics_before.get('cong_minutes', 0.0),

                metrics_after.get('cong_minutes', 0.0)

            )

        )

        cf.write(

            'Left-Behind Rate | '

            '{:.6f} | {:.6f}\n'.format(

                metrics_before.get('left_behind_rate', 0.0),

                metrics_after.get('left_behind_rate', 0.0)

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

        f"[Saved] Metrics comparison -> "

        f"{output_file}"

    )