from backlog_loader import get_backlog_matrix
from config import MAIN_DATA_FILE, S_D, _n_bins
from demand_processing import build_phi_from_excel, load_od_share_per_bin
from input_loader import load_demand_matrix, load_dwell_time, load_running_time

from .schema import ModelInputData


def load_excel_inputs(path=MAIN_DATA_FILE):
    r = load_running_time(path)
    e = load_dwell_time(path)
    p = load_demand_matrix(path)

    phi = build_phi_from_excel(
        n_stations=S_D,
        n_bins=_n_bins
    )

    od_share_per_bin = load_od_share_per_bin(
        demand=p,
        n_stations=S_D,
        n_bins=_n_bins
    )

    P_backlog = get_backlog_matrix(
        p,
        excel_path=path
    )

    return ModelInputData(
        r=r,
        e=e,
        p=p,
        phi=phi,
        P_backlog=P_backlog,
        od_share_per_bin=od_share_per_bin
    )

