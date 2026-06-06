import os

from config import DATA_FILE, S_D

from .afc_adapter import load_afc_demand
from .excel_adapter import load_excel_inputs
from .gtfs_adapter import load_gtfs_network
from .schema import ModelInputData, dense_od_matrix, zero_link_times


def load_model_inputs(input_type="excel", input_path=None, afc_path=None):
    input_path = input_path or DATA_FILE

    if input_type == "excel":
        return load_excel_inputs(input_path)

    if input_type != "gtfs_afc":
        raise ValueError("input_type must be one of: excel, gtfs_afc")

    gtfs = load_gtfs_network(input_path)

    if afc_path is None:
        candidate = os.path.join(input_path, "afc_taps.csv")
        afc_path = candidate if os.path.exists(candidate) else None

    afc = load_afc_demand(
        afc_path=afc_path,
        station_map=gtfs["station_map"]
    )

    pr = gtfs["r"]

    return ModelInputData(
        pr=pr,
        ac=zero_link_times(pr),
        dc=zero_link_times(pr),
        e=gtfs["e"],
        p=dense_od_matrix(afc["p"], S_D),
        phi=afc["phi"],
        P_backlog=dense_od_matrix(afc["P_backlog"], S_D),
        station_map=gtfs["station_map"],
        line_direction=gtfs["line_direction"],
        service_candidates=gtfs["service_candidates"]
    )
