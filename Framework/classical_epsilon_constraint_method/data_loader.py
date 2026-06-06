# data_loader.py

import pandas as pd


def load_running_times(file_path):
    """
    Load running time data from Excel.
    Returns dictionary:
    r[(i,j)] = running time
    """

    df = pd.read_excel(
        file_path,
        sheet_name='Running_time'
    )

    r = {
        (a, b): c
        for a, b, c in df.values
    }

    return r


def load_dwelling_times(file_path):
    """
    Load dwelling time data from Excel.
    Returns dictionary:
    e[i] = dwelling time
    """

    df = pd.read_excel(
        file_path,
        sheet_name='Dwelling_time'
    )

    e = {
        a: c
        for a, c in df.values
    }

    return e


def load_demand(file_path):
    """
    Load passenger demand data from Excel.
    Returns dictionary:
    p[(i,j)] = passenger demand
    """

    df = pd.read_excel(
        file_path,
        sheet_name='Demand_30(7.30-8.00am)mint'
    )

    p = {
        (a, b): c
        for a, b, c in df.values
    }

    return p


def load_all_data(file_path):
    """
    Load all input data together.
    """

    r = load_running_times(file_path)

    e = load_dwelling_times(file_path)

    p = load_demand(file_path)

    return {
        'r': r,
        'e': e,
        'p': p
    }