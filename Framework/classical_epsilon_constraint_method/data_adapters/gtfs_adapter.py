import os
from collections import defaultdict

import pandas as pd

from config import S_D


def _seconds(value):
    if pd.isna(value):
        return None

    parts = str(value).split(":")
    if len(parts) != 3:
        return None

    h, m, s = (int(float(p)) for p in parts)
    return h * 3600 + m * 60 + s


def _read_required(folder, filename):
    path = os.path.join(folder, filename)

    if not os.path.exists(path):
        raise FileNotFoundError(f"Required GTFS file not found: {path}")

    return pd.read_csv(path)


def load_gtfs_network(gtfs_folder, n_stations=S_D):
    stops = _read_required(gtfs_folder, "stops.txt")
    trips = _read_required(gtfs_folder, "trips.txt")
    stop_times = _read_required(gtfs_folder, "stop_times.txt")

    stop_ids = list(stops["stop_id"].astype(str))[:n_stations]
    station_map = {
        stop_id: idx + 1
        for idx, stop_id in enumerate(stop_ids)
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

    for trip_id, group in stop_times.groupby("trip_id"):
        sequence = []

        for row in group.sort_values("stop_sequence").itertuples(index=False):
            station = station_map[str(row.stop_id)]
            arrival = int(row.arrival_sec)
            departure = int(row.departure_sec)
            sequence.append((station, arrival, departure))
            e_samples[station].append(max(0.0, departure - arrival))

        for previous, current in zip(sequence, sequence[1:]):
            i, _, dep_i = previous
            j, arr_j, _ = current

            if j == i + 1:
                r_samples[(i, j)].append(max(0.0, arr_j - dep_i))

        if sequence:
            service_candidates.append(
                {
                    "trip_id": trip_id,
                    "stations": [station for station, _, _ in sequence],
                    "start_time": sequence[0][2],
                    "end_time": sequence[-1][1],
                }
            )

    r = {}
    for i in range(1, n_stations):
        samples = r_samples.get((i, i + 1), [])
        r[(i, i + 1)] = float(pd.Series(samples).median()) if samples else 120.0

    e = {}
    for i in range(1, n_stations + 1):
        samples = e_samples.get(i, [])
        dwell = float(pd.Series(samples).median()) if samples else 20.0
        e[i] = dwell if dwell > 0 else 20.0

    line_direction = {}
    for row in trips.itertuples(index=False):
        direction = getattr(row, "direction_id", None)
        line_direction[str(row.trip_id)] = direction

    return {
        "r": r,
        "e": e,
        "station_map": station_map,
        "line_direction": line_direction,
        "service_candidates": service_candidates,
    }

