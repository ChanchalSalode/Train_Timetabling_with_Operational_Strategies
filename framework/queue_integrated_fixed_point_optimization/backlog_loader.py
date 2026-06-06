# backlog_loader.py

import pandas as pd

from config import MAIN_DATA_FILE


# =========================================
# LOAD BACKLOG DATA
# =========================================

def load_backlog_from_excel(
    path=MAIN_DATA_FILE,
    sheet_name='Backlog_at_H'
):
    """
    Load initial passenger backlog matrix.

    Returns
    -------
    dict or None
        {(i,j): backlog passengers}
    """

    try:

        df = pd.read_excel(
            path,
            sheet_name=sheet_name
        )

        backlog = {}

        for _, row in df.iterrows():

            backlog[
                (
                    int(row['i']),
                    int(row['j'])
                )
            ] = float(row['P_ij'])

        print(
            "Backlog data loaded successfully."
        )

        return backlog

    except Exception as e:

        print(
            f"[WARNING] Failed to load backlog data: {e}"
        )

        return None


# =========================================
# BUILD DEFAULT BACKLOG
# =========================================

def build_default_backlog(
    demand,
    pre_back_minutes=2.0,
    horizon_minutes=30.0
):
    """
    Build fallback backlog approximation
    when Excel sheet is unavailable.

    Parameters
    ----------
    demand : dict
        OD demand dictionary

    pre_back_minutes : float
        Assumed pre-existing backlog duration

    horizon_minutes : float
        Total demand horizon duration

    Returns
    -------
    dict
        {(i,j): backlog passengers}
    """

    factor = (
        float(pre_back_minutes)
        /
        float(horizon_minutes)
    )

    backlog = {}

    for (i, j), val in demand.items():

        backlog[
            (
                int(i),
                int(j)
            )
        ] = factor * float(val)

    print(
        "Default backlog approximation created."
    )

    return backlog


# =========================================
# GET BACKLOG MATRIX
# =========================================

def get_backlog_matrix(
    demand,
    excel_path=MAIN_DATA_FILE,
    sheet_name='Backlog_at_H',
    pre_back_minutes=2.0,
    horizon_minutes=30.0
):
    """
    Unified interface for backlog loading.

    First tries Excel.
    If unavailable, builds default backlog.
    """

    backlog = load_backlog_from_excel(
        path=excel_path,
        sheet_name=sheet_name
    )

    if backlog is None:

        backlog = build_default_backlog(
            demand=demand,
            pre_back_minutes=pre_back_minutes,
            horizon_minutes=horizon_minutes
        )

    return backlog
