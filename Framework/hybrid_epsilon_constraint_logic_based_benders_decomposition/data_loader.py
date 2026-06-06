# data_loader.py

import pandas as pd

from config import (
    DATA_FILE
)


# =========================================
# PURE RUNNING TIME
# =========================================

def load_pure_running_time():

    df = pd.read_excel(
        DATA_FILE,
        sheet_name='Pure_running_time'
    )

    return {

        (a, b): c

        for a, b, c
        in df.values

    }


# =========================================
# ACCELERATION TIME
# =========================================

def load_acceleration_time():

    df = pd.read_excel(
        DATA_FILE,
        sheet_name='Acceleration_time'
    )

    return {

        (a, b): c

        for a, b, c
        in df.values

    }


# =========================================
# DECELERATION TIME
# =========================================

def load_deceleration_time():

    df = pd.read_excel(
        DATA_FILE,
        sheet_name='Deceleration_time'
    )

    return {

        (a, b): c

        for a, b, c
        in df.values

    }


# =========================================
# DEMAND DATA
# =========================================

def load_demand_data():

    df = pd.read_excel(
        DATA_FILE,
        sheet_name='Demand_30(7.30-8.00am)mint_peak'
    )

    return {

        (a, b): c

        for a, b, c
        in df.values

    }


# =========================================
# DWELLING TIME
# =========================================

def load_dwelling_time():

    df = pd.read_excel(
        DATA_FILE,
        sheet_name='Dwelling_time'
    )

    return {

        a: c

        for a, c
        in df.values

    }


# =========================================
# LOAD ALL DATA
# =========================================

def load_all_data():

    pr = load_pure_running_time()

    ac = load_acceleration_time()

    dc = load_deceleration_time()

    p = load_demand_data()

    e = load_dwelling_time()

    return {

        'pr': pr,

        'ac': ac,

        'dc': dc,

        'p': p,

        'e': e

    }