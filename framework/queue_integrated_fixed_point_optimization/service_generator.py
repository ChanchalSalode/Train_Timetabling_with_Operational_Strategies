# service_generator.py

import pandas as pd

from helpers import var_value_safe


def generate_initial_services(
    var_dicts,
    H,
    C
):

    # =========================================
    # Variables
    # =========================================

    d = var_dicts.get('d', {})

    x = var_dicts.get('x', {})

    # =========================================
    # Sets
    # =========================================

    K_u = var_dicts.get('K_u', [])
    K_d = var_dicts.get('K_d', [])

    S_u = var_dicts.get('S_u', [])
    S_d = var_dicts.get('S_d', [])

    # =========================================
    # Service Storage
    # =========================================

    services_rows = []

    # =========================================
    # Up Direction Services
    # =========================================

    for k in K_u:

        for i in S_u:

            # =================================
            # Stop Decision
            # =================================

            xstop = 1

            try:

                if (
                    (k, i) in x.keys()
                    and
                    hasattr(x[(k, i)], 'X')
                ):

                    xstop = int(
                        round(
                            x[(k, i)].X
                        )
                    )

                elif (k, i) in x.keys():

                    xstop = int(
                        x[(k, i)].LB
                    )

            except Exception:

                xstop = 1

            # =================================
            # Departure Time
            # =================================

            depart_min = 0.0

            try:

                if (
                    (k, i) in d.keys()
                    and
                    hasattr(d[(k, i)], 'X')
                ):

                    depart_min = max(

                        0.0,

                        (
                            float(d[(k, i)].X)
                            - float(H)
                        ) / 60.0

                    )

            except Exception:

                depart_min = 0.0

            # =================================
            # Store Service
            # =================================

            services_rows.append({

                'station_i': i,

                'k': k,

                'depart_time_min': depart_min,

                'x_stop': xstop,

                'capacity': C,

                'onboard_arrival': 0.0

            })

    # =========================================
    # Down Direction Services
    # =========================================

    for k in K_d:

        for i in S_d:

            # =================================
            # Stop Decision
            # =================================

            xstop = 1

            try:

                if (
                    (k, i) in x.keys()
                    and
                    hasattr(x[(k, i)], 'X')
                ):

                    xstop = int(
                        round(
                            x[(k, i)].X
                        )
                    )

                elif (k, i) in x.keys():

                    xstop = int(
                        x[(k, i)].LB
                    )

            except Exception:

                xstop = 1

            # =================================
            # Departure Time
            # =================================

            depart_min = 0.0

            try:

                if (
                    (k, i) in d.keys()
                    and
                    hasattr(d[(k, i)], 'X')
                ):

                    depart_min = max(

                        0.0,

                        (
                            float(d[(k, i)].X)
                            - float(H)
                        ) / 60.0

                    )

            except Exception:

                depart_min = 0.0

            # =================================
            # Store Service
            # =================================

            services_rows.append({

                'station_i': i,

                'k': k,

                'depart_time_min': depart_min,

                'x_stop': xstop,

                'capacity': C,

                'onboard_arrival': 0.0

            })

    # =========================================
    # Convert to DataFrame
    # =========================================

    services_df = pd.DataFrame(
        services_rows
    )

    return services_df


def update_services_from_solution(
    services_src,
    var_dicts,
    H,
    T_WINDOW
):
    """Update service departures and stopping decisions from model variables."""

    d = var_dicts.get('d', {})
    x = var_dicts.get('x', {})

    if isinstance(services_src, str):

        try:

            services_df = pd.read_csv(services_src)

        except Exception:

            return

    else:

        services_df = services_src

    for idx, row in services_df.iterrows():

        try:

            k = row['k']
            station = row['station_i']

            if (k, station) in d.keys():

                depart_min = (
                    var_value_safe(d[(k, station)])
                    - float(H)
                ) / 60.0

                depart_min = max(0.0, depart_min)

                if depart_min > T_WINDOW / 60.0:

                    depart_min = T_WINDOW / 60.0

                services_df.at[idx, 'depart_time_min'] = depart_min

            if (k, station) in x.keys():

                services_df.at[idx, 'x_stop'] = int(
                    round(
                        var_value_safe(
                            x[(k, station)],
                            row.get('x_stop', 1)
                        )
                    )
                )

        except Exception:

            continue

    if isinstance(services_src, str):

        services_df.to_csv(
            services_src,
            index=False
        )
