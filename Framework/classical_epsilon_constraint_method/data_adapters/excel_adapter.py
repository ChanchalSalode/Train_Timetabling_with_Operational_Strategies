from config import INPUT_FILE, S_D
from data_loader import load_demand, load_dwelling_times, load_running_times

from .schema import ModelInputData, default_backlog, dense_od_matrix, uniform_phi


def load_excel_inputs(path=INPUT_FILE):
    r = load_running_times(path)
    e = load_dwelling_times(path)
    p = dense_od_matrix(
        load_demand(path),
        S_D
    )

    return ModelInputData(
        r=r,
        e=e,
        p=p,
        phi=uniform_phi(S_D, 30),
        P_backlog=default_backlog(p)
    )

