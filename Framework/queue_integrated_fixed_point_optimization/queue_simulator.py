# queue_simulator.py

import math

import pandas as pd

from collections import defaultdict


# =========================================
# INTEGRATE ARRIVALS OVER TIME INTERVAL
# =========================================

def integrate_bin(
    arrivals,
    t1_min,
    t2_min
):

    if t2_min <= t1_min:

        return 0.0

    total = 0.0

    start_min = int(
        math.floor(t1_min)
    )

    end_min = int(
        math.floor(t2_min)
    )

    # =====================================
    # Same bin
    # =====================================

    if start_min == end_min:

        frac = (
            t2_min - t1_min
        )

        if 0 <= start_min < len(arrivals):

            return (

                arrivals[start_min]
                * frac

            )

        return 0.0

    # =====================================
    # First partial minute
    # =====================================

    first_frac = (

        start_min + 1
        - t1_min

    )

    if 0 <= start_min < len(arrivals):

        total += (

            arrivals[start_min]
            * first_frac

        )

    # =====================================
    # Full minutes
    # =====================================

    for m in range(

        start_min + 1,
        end_min

    ):

        if 0 <= m < len(arrivals):

            total += arrivals[m]

    # =====================================
    # Last partial minute
    # =====================================

    last_frac = (

        t2_min - end_min

    )

    if 0 <= end_min < len(arrivals):

        total += (

            arrivals[end_min]
            * last_frac

        )

    return total


# =========================================
# FLUID QUEUE SIMULATOR
# =========================================

def simulate_queues_station(

    arrivals_by_od,
    services_at_station,

    mu_board=0.5,

    dwell_base_sec=20.0,

    dwell_per_pax_sec=0.25,

    initial_backlog=None

):

    # =====================================
    # Boarding rate
    # =====================================

    mu_board_per_min = (
        mu_board * 60.0
    )

    # =====================================
    # Queue state
    # =====================================

    Q = defaultdict(float)

    # =====================================
    # Initial backlog
    # =====================================

    if initial_backlog:

        for od, val in initial_backlog.items():

            try:

                Q[od] = float(val)

            except Exception:

                continue

    # =====================================
    # Sort services
    # =====================================

    services_sorted = sorted(

        services_at_station,

        key=lambda s:
        s['depart_time_min']

    )

    last_time = 0.0

    records = []

    # =====================================
    # OD pairs
    # =====================================

    od_keys = set(
        arrivals_by_od.keys()
    )

    if initial_backlog:

        od_keys.update(
            initial_backlog.keys()
        )

    # =====================================
    # Process each train
    # =====================================

    for s in services_sorted:

        k = s['k']

        depart_time = float(
            s['depart_time_min']
        )

        x_stop = int(
            s.get('x_stop', 1)
        )

        capacity = float(
            s.get('capacity', 1e9)
        )

        onboard_arrival = float(

            s.get(
                'onboard_arrival',
                0.0
            )

        )

        # =================================
        # New arrivals
        # =================================

        A_by_od = {}

        for od in od_keys:

            series = arrivals_by_od.get(
                od,
                None
            )

            if series is None:

                A = 0.0

            else:

                try:

                    A = integrate_bin(

                        series,

                        last_time,

                        depart_time

                    )

                except Exception:

                    A = 0.0

            A_by_od[od] = A

            Q[od] = (

                Q.get(od, 0.0)
                + A

            )

        # =================================
        # Queue before boarding
        # =================================

        Q_before = dict(Q)

        nb_by_od = {

            od: 0.0

            for od in od_keys

        }

        v_by_od = {

            od: Q_before.get(
                od,
                0.0
            )

            for od in od_keys

        }

        nb_total = 0.0

        # =================================
        # Skip-stop train
        # =================================

        if x_stop == 0:

            pass

        else:

            total_queue = sum(

                Q_before.get(
                    od,
                    0.0
                )

                for od in od_keys

            )

            if total_queue > 0.0:

                # =========================
                # Initial dwell estimate
                # =========================

                dwell_min = (
                    dwell_base_sec / 60.0
                )

                approx_capacity_by_time = (

                    mu_board_per_min
                    * dwell_min

                )

                approx_capacity = min(

                    max(
                        0.0,
                        capacity
                        - onboard_arrival
                    ),

                    approx_capacity_by_time

                )

                pred_board = min(

                    total_queue,
                    approx_capacity

                )

                # =========================
                # Updated dwell time
                # =========================

                dwell_sec = (

                    dwell_base_sec

                    +

                    dwell_per_pax_sec
                    * pred_board

                )

                dwell_min = (
                    dwell_sec / 60.0
                )

                capacity_by_time = (

                    mu_board_per_min
                    * dwell_min

                )

                C = max(

                    0.0,

                    min(

                        max(
                            0.0,
                            capacity
                            - onboard_arrival
                        ),

                        capacity_by_time

                    )

                )

                # =========================
                # Positive queues
                # =========================

                positive_queue_od = [

                    od

                    for od in od_keys

                    if Q_before.get(
                        od,
                        0.0
                    ) > 0

                ]

                total_pos_queue = sum(

                    Q_before.get(
                        od,
                        0.0
                    )

                    for od in positive_queue_od

                )

                # =========================
                # Allocate boarding
                # =========================

                if total_pos_queue > 0.0:

                    for od in od_keys:

                        qval = Q_before.get(
                            od,
                            0.0
                        )

                        if qval <= 0.0:

                            nb_by_od[od] = 0.0

                        else:

                            alloc = (

                                C

                                *

                                (

                                    qval
                                    / total_pos_queue

                                )

                            )

                            nb_by_od[od] = min(
                                qval,
                                alloc
                            )

                    nb_total = sum(
                        nb_by_od.values()
                    )

                    # =====================
                    # Remaining capacity
                    # =====================

                    leftover = (
                        C - nb_total
                    )

                    if leftover > 0.5:

                        rem_list = sorted(

                            positive_queue_od,

                            key=lambda od:

                            (

                                Q_before.get(
                                    od,
                                    0.0
                                )

                                -

                                nb_by_od.get(
                                    od,
                                    0.0
                                )

                            ),

                            reverse=True

                        )

                        for od in rem_list:

                            if leftover <= 0.5:

                                break

                            can_take = (

                                Q_before.get(
                                    od,
                                    0.0
                                )

                                -

                                nb_by_od.get(
                                    od,
                                    0.0
                                )

                            )

                            give = min(
                                can_take,
                                leftover
                            )

                            nb_by_od[od] += give

                            leftover -= give

                        nb_total = sum(
                            nb_by_od.values()
                        )

                    # =====================
                    # Left behind
                    # =====================

                    v_by_od = {

                        od:

                        Q_before.get(
                            od,
                            0.0
                        )

                        -

                        nb_by_od.get(
                            od,
                            0.0
                        )

                        for od in od_keys

                    }

        # =================================
        # Store results
        # =================================

        for od in sorted(Q_before.keys()):

            i, j = od

            records.append({

                'k': k,

                'i': i,

                'j': j,

                'depart_time_min': depart_time,

                'w_before':

                float(
                    Q_before.get(
                        od,
                        0.0
                    )
                ),

                'A_arrivals':

                float(
                    A_by_od.get(
                        od,
                        0.0
                    )
                ),

                'nb_boarded':

                float(
                    nb_by_od.get(
                        od,
                        0.0
                    )
                ),

                'v_left':

                float(
                    v_by_od.get(
                        od,
                        0.0
                    )
                ),

                'nb_total':

                float(nb_total)

            })

        # =================================
        # Update queue
        # =================================

        for od in Q_before.keys():

            Q[od] = v_by_od.get(

                od,

                Q_before.get(
                    od,
                    0.0
                )

            )

        last_time = depart_time

    # =====================================
    # Return DataFrame
    # =====================================

    df = pd.DataFrame.from_records(
        records
    )

    return df