# data_loader.py

import pandas as pd


def load_running_times(path):

    df = pd.read_excel(path, 'Running_time')

    return {
        (a, b): c
        for a, b, c in df.values
    }


def load_dwelling_times(path):

    df = pd.read_excel(path, 'Dwelling_time')

    return {
        a: c
        for a, c in df.values
    }


def load_demand(path):

    df = pd.read_excel(
        path,
        'Demand_30(7.30-8.00am)mint'
    )

    return {
        (a, b): c
        for a, b, c in df.values
    }