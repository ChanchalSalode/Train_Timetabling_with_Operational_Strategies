# performance_metrics.py


from config import (
    C,
    h_min
)


# =========================================
# ANALYZE PASSENGER SOLUTION
# =========================================

def analyze_passenger_solution(
    solution_id,
    tau_val,
    x_val,
    d_val,
    z_val,
    h_val,
    sub_vars,
    output_prefix="FINAL"
):

    # =====================================
    # Extract Subproblem Variables
    # =====================================

    w = sub_vars["w"]

    w_b = sub_vars["w_b"]

    n1 = sub_vars["n1"]

    v = sub_vars["v"]

    n_b1 = sub_vars["n_b1"]

    # =====================================
    # Constants
    # =====================================

    eps = 1e-6

    cap = C

    # =====================================
    # 1. Peak Load Factor (PLF)
    # =====================================

    plf = max(

        (

            n1[(k, i)] / cap

            for (k, i)
            in n1

            if tau_val.get(k, 0) > 0.5

        ),

        default=0.0

    )

    # =====================================
    # Overloaded Services
    # =====================================

    overloaded = [

        (
            k,
            i,
            n1[(k, i)] / cap
        )

        for (k, i)
        in n1

        if (
            tau_val.get(k, 0) > 0.5
            and
            n1[(k, i)] / cap > 0.85
        )

    ]

    # =====================================
    # 2. Congestion Duration
    # =====================================

    cong_seconds = sum(

        h_val.get(
            (k - 1, k),
            h_min
        )

        for (k, _, _)
        in overloaded

    )

    cong_minutes = cong_seconds / 60.0

    # =====================================
    # 3. Left-Behind Passenger Rate
    # =====================================

    sum_v = sum(

        v[key]

        for key in v

        if (
            x_val.get(
                (key[0], key[1]),
                0
            ) > 0.5

            and

            x_val.get(
                (key[0], key[2]),
                0
            ) > 0.5
        )

    )

    sum_w = sum(

        w[key]

        for key in w

        if (
            x_val.get(
                (key[0], key[1]),
                0
            ) > 0.5

            and

            x_val.get(
                (key[0], key[2]),
                0
            ) > 0.5
        )

    )

    left_behind = (

        sum_v / sum_w

        if sum_w > eps

        else 0.0

    )

    # =====================================
    # 4. Demand-to-Capacity Ratio
    # =====================================

    dcr = {}

    for (k, i), arr in w_b.items():

        if tau_val.get(k, 0) < 0.5:

            continue

        hw = h_val.get(
            (k - 1, k),
            h_min
        )

        dcr[(k, i)] = (

            arr * hw / cap

        )

    peak_dcr = max(

        dcr.values(),

        default=0.0

    )

    # =====================================
    # Print Summary
    # =====================================

    print("\n====================================")
    print(f"{output_prefix} SOLUTION {solution_id}")
    print("====================================")

    print("\nPerformance Metrics:")

    print(
        f"Peak Load Factor (PLF): "
        f"{plf:.4f}"
    )

    print(
        f"Congestion Duration: "
        f"{cong_seconds:.2f} seconds "
        f"({cong_minutes:.2f} minutes)"
    )

    print(
        f"Left-Behind Passenger Rate: "
        f"{left_behind:.4f}"
    )

    print(
        f"Peak Demand-to-Capacity Ratio: "
        f"{peak_dcr:.4f}"
    )

    # =====================================
    # Overloaded Stations
    # =====================================

    if overloaded:

        print("\nOverloaded Stations:")

        for (k, i, val) in overloaded:

            print(

                f"(Train {k}, Station {i}) "
                f"PLF = {val:.4f}"

            )

    else:

        print(
            "\nNo overloaded stations."
        )

    # =====================================
    # Return Metrics
    # =====================================

    return {

        "plf": plf,

        "overloaded": overloaded,

        "cong_seconds": cong_seconds,

        "cong_minutes": cong_minutes,

        "left_behind": left_behind,

        "peak_dcr": peak_dcr

    }