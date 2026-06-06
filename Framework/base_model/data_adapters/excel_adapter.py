from config import DATA_FILE, S_dn
from data_loader import load_demand, load_dwelling_times, load_running_times

from .schema import ModelInputData, default_backlog, dense_od_matrix, uniform_phi


def load_excel_inputs(path=DATA_FILE):
    r = load_running_times(path)
    e = load_dwelling_times(path)
    p = dense_od_matrix(
        load_demand(path),
        S_dn
    )

    return ModelInputData(
        r=r,
        e=e,
        p=p,
        phi=uniform_phi(S_dn, 30),
        P_backlog=default_backlog(p)
    )
