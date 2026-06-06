import os
from collections import defaultdict

import pandas as pd

from config import BIN, DWELL_BASE_SEC, H, S_D, _n_bins


def _seconds(value):
    if pd.isna(value):
        return None

    parts = str(value).split(":")
    if len(parts) != 3:
        return None

    h, m, s = (int(float(p)) for p in parts)
    return h * 3600 + m * 60 + s


def _read_required_csv(folder, filename):
    path = os.path.join(folder, filename)

    if not os.path.exists(path):
        raise FileNotFoundError(f"Required GTFS file not found: {path}")

    return pd.read_csv(path)


def load_gtfs_network(gtfs_folder, n_stations=S_D):
    stops = _read_required_csv(gtfs_folder, "stops.txt")
    trips = _read_required_csv(gtfs_folder, "trips.txt")
    stop_times = _read_required_csv(gtfs_folder, "stop_times.txt")

    if "stop_id" not in stops.columns:
        raise ValueError("stops.txt must contain stop_id")

    ordered_stop_ids = list(stops["stop_id"].astype(str))
    ordered_stop_ids = ordered_stop_ids[:n_stations]

    station_map = {
        stop_id: idx + 1
        for idx, stop_id in enumerate(ordered_stop_ids)
    }

    stop_name_map = {}
    if "stop_name" in stops.columns:
        stop_name_map = {
            str(row.stop_id): row.stop_name
            for row in stops.itertuples(index=False)
        }

    stop_times = stop_times.copy()
    stop_times["stop_id"] = stop_times["stop_id"].astype(str)
    stop_times = stop_times[stop_times["stop_id"].isin(station_map)]
    stop_times["arrival_sec"] = stop_times["arrival_time"].map(_seconds)
    stop_times["departure_sec"] = stop_times["departure_time"].map(_seconds)
    stop_times = stop_times.dropna(
        subset=["arrival_sec", "departure_sec", "stop_sequence"]
    )

    r_samples = defaultdict(list)
    e_samples = defaultdict(list)
    service_candidates = []

    for trip_id, grp in stop_times.groupby("trip_id"):
        grp = grp.sort_values("stop_sequence")
        seq = []

        for row in grp.itertuples(index=False):
            station = station_map[str(row.stop_id)]
            arr = int(row.arrival_sec)
            dep = int(row.departure_sec)
            seq.append((station, arr, dep))
            e_samples[station].append(max(0.0, dep - arr))

        for prev, curr in zip(seq, seq[1:]):
            i, _, dep_prev = prev
            j, arr_curr, _ = curr

            if j == i + 1:
                r_samples[(i, j)].append(max(0.0, arr_curr - dep_prev))

        if seq:
            service_candidates.append(
                {
                    "trip_id": trip_id,
                    "stations": [station for station, _, _ in seq],
                    "start_time": seq[0][2],
                    "end_time": seq[-1][1],
                }
            )

    r = {}
    for i in range(1, n_stations):
        samples = r_samples.get((i, i + 1), [])
        r[(i, i + 1)] = float(pd.Series(samples).median()) if samples else 120.0

    e = {}
    for i in range(1, n_stations + 1):
        samples = e_samples.get(i, [])
        dwell = float(pd.Series(samples).median()) if samples else DWELL_BASE_SEC
        e[i] = dwell if dwell > 0 else DWELL_BASE_SEC

    line_direction = {}
    if "route_id" in trips.columns:
        for row in trips.itertuples(index=False):
            direction = getattr(row, "direction_id", None)
            line_direction[str(row.trip_id)] = direction

    return {
        "r": r,
        "e": e,
        "station_map": station_map,
        "station_names": stop_name_map,
        "line_direction": line_direction,
        "service_candidates": service_candidates,
    }

