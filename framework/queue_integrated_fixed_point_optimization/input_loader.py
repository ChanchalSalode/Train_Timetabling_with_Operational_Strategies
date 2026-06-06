# input_loader.py

import pandas as pd

from config import MAIN_DATA_FILE


# =========================================
# LOAD RUNNING TIME
# =========================================

def load_running_time(
    path=MAIN_DATA_FILE
):

    try:

        df = pd.read_excel(
            path,
            'Running_time'
        )

        r = dict(

            (
                tuple((a, b)),
                c
            )

            for a, b, c
            in df.values

        )

        print(
            "Running time data loaded."
        )

        return r

    except Exception as e:

        print(
            f"Failed to load running times: {e}"
        )

        return {}


# =========================================
# LOAD DWELL TIME
# =========================================

def load_dwell_time(
    path=MAIN_DATA_FILE
):

    try:

        df = pd.read_excel(
            path,
            'Dwelling_time'
        )

        my_tuple = [

            tuple(x)

            for x
            in df.values

        ]

        e = dict(

            (a, c)

            for a, c
            in my_tuple

        )

        print(
            "Dwelling time data loaded."
        )

        return e

    except Exception as ex:

        print(
            f"Failed to load dwell times: {ex}"
        )

        return {}


# =========================================
# LOAD DEMAND MATRIX
# =========================================

def load_demand_matrix(
    path=MAIN_DATA_FILE
):

    try:

        df = pd.read_excel(
            path,
            'Demand_30(7.30-8.00am)mint'
        )

        p = dict(

            (
                tuple((a, b)),
                c
            )

            for a, b, c
            in df.values

        )

        print(
            "Demand matrix loaded."
        )

        return p

    except Exception as ex:

        print(
            f"Failed to load demand matrix: {ex}"
        )

        return {}


# =========================================
# LOAD ALL INPUTS
# =========================================

def load_all_inputs(
    path=MAIN_DATA_FILE
):

    print(
        "\n================================="
    )

    print(
        "LOADING INPUT DATA"
    )

    print(
        "================================="
    )

    r = load_running_time(path)

    e = load_dwell_time(path)

    p = load_demand_matrix(path)

    print(
        "\nAll input data loaded successfully."
    )

    return {

        'r': r,
        'e': e,
        'p': p

    }


# =========================================
# LOAD SERVICES CSV
# =========================================

def load_services_data(
    csv_path
):

    try:

        services_df = pd.read_csv(
            csv_path
        )

        print(
            f"Services data loaded from {csv_path}"
        )

        return services_df

    except Exception as ex:

        print(
            f"Failed to load services data: {ex}"
        )

        return pd.DataFrame()


# =========================================
# LOAD EXCEL SHEET
# =========================================

def load_excel_sheet(
    file_path,
    sheet_name
):

    try:

        df = pd.read_excel(
            file_path,
            sheet_name
        )

        print(
            f"Loaded sheet: {sheet_name}"
        )

        return df

    except Exception as ex:

        print(
            f"Failed to load sheet {sheet_name}: {ex}"
        )

        return pd.DataFrame()
