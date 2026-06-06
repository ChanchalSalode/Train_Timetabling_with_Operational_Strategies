# demand_processing.py

import numpy as np
import pandas as pd

from config import TIME_DEPENDENT_FILE


# =========================================
# BUILD PHI FROM EXCEL
# =========================================

def build_phi_from_excel(
    path=TIME_DEPENDENT_FILE,
    sheet_name='Arrival_profile_per_min',
    n_stations=16,
    n_bins=30
):
    """
    Load station-wise temporal arrival profile.

    Returns
    -------
    dict
        phi[i] = list of arrival shares by bin
    """

    try:

        df_phi = pd.read_excel(
            path,
            sheet_name
        )

        col_station = (

            'station'

            if 'station' in df_phi.columns

            else (

                'i'

                if 'i' in df_phi.columns

                else None

            )

        )

        col_bin = (

            'bin_idx'

            if 'bin_idx' in df_phi.columns

            else (

                'minute_idx'

                if 'minute_idx' in df_phi.columns

                else None

            )

        )

        col_phi = 'phi'

        if (
            col_station is None
            or
            col_bin is None
            or
            col_phi not in df_phi.columns
        ):

            raise ValueError(
                "Arrival profile columns not recognized"
            )

        phi = {}

        for i in range(1, n_stations + 1):

            vec = [0.0] * n_bins

            sub = df_phi[
                df_phi[col_station] == i
            ]

            for _, row in sub.iterrows():

                b = int(row[col_bin])

                if 0 <= b < n_bins:

                    vec[b] = float(
                        row[col_phi]
                    )

            s = sum(vec)

            phi[i] = [

                v / s

                if s > 0

                else 1.0 / n_bins

                for v in vec

            ]

        return phi

    except Exception as e:

        print(
            f"[WARNING] Failed to load phi: {e}"
        )

        return {

            i: [1.0 / n_bins] * n_bins

            for i in range(1, n_stations + 1)

        }


# =========================================
# BUILD PWL BREAKPOINTS
# =========================================

def build_pwl_breaks_for_dir(
    demand,
    phi,
    stations,
    times,
    n_bins,
    od_share_per_bin=None
):
    """
    Construct PWL cumulative demand functions.

    Parameters
    ----------
    demand : dict
        OD demand matrix.
    phi : dict
        Station-wise temporal profile.
    stations : list
        Stations in the direction.
    times : list
        Time breakpoints.
    n_bins : int
        Number of bins.
    od_share_per_bin : dict | None
        Optional OD-specific temporal shares.

    Returns
    -------
    dict
        (i,j) -> (x_breaks, y_breaks)
    """

    pwl = {}

    for i in stations:

        J = [

            j

            for j in stations

            if j > i
            and (i, j) in demand

        ]

        if not J:

            continue

        row_sum = sum(

            demand.get((i, j), 0.0)

            for j in J

        )

        if row_sum <= 0:

            for j in J:

                pwl[(i, j)] = (

                    times,

                    [0.0] * (n_bins + 1)

                )

            continue

        bin_mass_i = [

            row_sum * phi[i][b]

            for b in range(n_bins)

        ]

        for j in J:

            pij = demand.get((i, j), 0.0)

            if pij <= 0:

                pwl[(i, j)] = (

                    times,

                    [0.0] * (n_bins + 1)

                )

                continue

            if (
                od_share_per_bin is not None
                and (i, j) in od_share_per_bin
            ):

                s_vec = od_share_per_bin[(i, j)]

                raw = [

                    bin_mass_i[b]
                    *
                    max(0.0, float(s_vec[b]))

                    for b in range(n_bins)

                ]

                raw_sum = sum(raw)

                if raw_sum > 0.0:

                    mass_ij_per_bin = [

                        raw[b] * (pij / raw_sum)

                        for b in range(n_bins)

                    ]

                else:

                    mass_ij_per_bin = [

                        pij / float(n_bins)

                    ] * n_bins

            else:

                share = (

                    pij / row_sum

                    if row_sum > 0

                    else 0.0

                )

                mass_ij_per_bin = [

                    bin_mass_i[b] * share

                    for b in range(n_bins)

                ]

            cum = np.cumsum(
                mass_ij_per_bin
            )

            y_breaks = [0.0] + list(cum)

            y_breaks[-1] = pij

            pwl[(i, j)] = (

                times,
                y_breaks

            )

    return pwl


# =========================================
# CONVERT PWL TO MINUTE SERIES
# =========================================

def pwl_to_minute_series(
    pwl_entry
):
    """
    Convert cumulative PWL to arrivals/bin.
    """

    x_breaks, y_breaks = pwl_entry

    masses = []

    for b in range(1, len(y_breaks)):

        masses.append(

            float(

                y_breaks[b]
                - y_breaks[b - 1]

            )

        )

    return masses


# =========================================
# LOAD OD TEMPORAL SHARES
# =========================================

def load_od_share_per_bin(
    path=TIME_DEPENDENT_FILE,
    sheet_name='OD_share_per_bin',
    demand=None,
    n_stations=16,
    n_bins=30
):
    """
    Load OD-specific temporal demand shares.

    Returns
    -------
    dict | None
        (i,j) -> list of shares by bin
    """

    try:

        df_share = pd.read_excel(
            path,
            sheet_name
        )

        temp = {}

        if demand is not None:

            for (i, j), _ in demand.items():

                temp[(int(i), int(j))] = [0.0] * n_bins

        for _, row in df_share.iterrows():

            i = int(row['i'])
            j = int(row['j'])
            b = int(row['bin_idx'])
            val = float(row['share'])

            if (
                (i, j) in temp
                and 0 <= b < n_bins
            ):

                temp[(i, j)][b] = max(0.0, val)

        for i in range(1, n_stations + 1):

            Js = [
                jj
                for (ii, jj) in temp.keys()
                if ii == i and jj > i
            ]

            for b in range(n_bins):

                ssum = sum(
                    temp[(i, jj)][b]
                    for jj in Js
                ) if Js else 0.0

                if ssum > 0:

                    for jj in Js:

                        temp[(i, jj)][b] /= ssum

                else:

                    row_sum = sum(
                        demand.get((i, jj), 0.0)
                        for jj in Js
                    ) if demand is not None else 0.0

                    for jj in Js:

                        temp[(i, jj)][b] = (
                            demand.get((i, jj), 0.0) / row_sum
                            if row_sum > 0
                            else 0.0
                        )

        return temp

    except Exception as e:

        print(
            f"[WARNING] Failed to load OD shares: {e}"
        )

        return None


# =========================================
# BUILD ARRIVALS BY STATION
# =========================================

def build_arrivals_by_station(
    pwl_up,
    pwl_down
):
    """
    Create arrivals grouped by origin station.
    """

    arrivals_by_station_od = {}

    for (i, j), pwl in pwl_up.items():

        arrivals_by_station_od.setdefault(
            i,
            {}
        )[(i, j)] = pwl_to_minute_series(
            pwl
        )

    for (i, j), pwl in pwl_down.items():

        arrivals_by_station_od.setdefault(
            i,
            {}
        )[(i, j)] = pwl_to_minute_series(
            pwl
        )

    return arrivals_by_station_od
