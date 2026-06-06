import os

from config import MAIN_DATA_FILE

from .afc_adapter import load_afc_demand
from .excel_adapter import load_excel_inputs
from .gtfs_adapter import load_gtfs_network
from .schema import ModelInputData


def load_model_inputs(input_type="excel", input_path=MAIN_DATA_FILE, afc_path=None):
    input_path = input_path or MAIN_DATA_FILE

    if input_type == "excel":
        return load_excel_inputs(input_path)

    if input_type != "gtfs_afc":
        raise ValueError(
            "input_type must be one of: excel, gtfs_afc"
        )

    gtfs_data = load_gtfs_network(input_path)

    if afc_path is None:
        candidate = os.path.join(input_path, "afc_taps.csv")
        afc_path = candidate if os.path.exists(candidate) else None

    afc_data = load_afc_demand(
        afc_path=afc_path,
        station_map=gtfs_data["station_map"]
    )

    return ModelInputData(
        r=gtfs_data["r"],
        e=gtfs_data["e"],
        p=afc_data["p"],
        phi=afc_data["phi"],
        P_backlog=afc_data["P_backlog"],
        od_share_per_bin=afc_data["od_share_per_bin"],
        station_map=gtfs_data["station_map"],
        line_direction=gtfs_data["line_direction"],
        service_candidates=gtfs_data["service_candidates"]
    )
