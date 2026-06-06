from config import DATA_FILE, S_D
from data_loader import (
    load_acceleration_time,
    load_deceleration_time,
    load_demand_data,
    load_dwelling_time,
    load_pure_running_time
)

from .schema import ModelInputData, default_backlog, dense_od_matrix, uniform_phi


def load_excel_inputs(path=DATA_FILE):
    # Existing data_loader functions use DATA_FILE from config; path is retained
    # for interface consistency across adapters.
    pr = load_pure_running_time()
    ac = load_acceleration_time()
    dc = load_deceleration_time()
    p = dense_od_matrix(
        load_demand_data(),
        S_D
    )
    e = load_dwelling_time()

    return ModelInputData(
        pr=pr,
        ac=ac,
        dc=dc,
        e=e,
        p=p,
        phi=uniform_phi(S_D, 30),
        P_backlog=default_backlog(p)
    )

