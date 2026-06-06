import os

import pandas as pd

from config import S_D

from .schema import default_backlog, normalize_phi, uniform_phi


def _seconds(value):
    if pd.isna(value):
        return None

    if isinstance(value, (int, float)):
        return int(value)

    parts = str(value).split(":")
    if len(parts) != 3:
        return None

    h, m, s = (int(float(p)) for p in parts)
    return h * 3600 + m * 60 + s


def _station(value, station_map):
    text = str(value)

    if text in station_map:
        return int(station_map[text])

    return int(float(value))


def load_afc_demand(
    afc_path,
    station_map=None,
    n_stations=S_D,
    horizon_start=27000,
    horizon_seconds=1800,
    n_bins=30
):
    station_map = station_map or {}

    if afc_path is None or not os.path.exists(afc_path):
        print("[WARNING] AFC file not found; using zero demand.")
        return {
            "p": {},
            "phi": uniform_phi(n_stations, n_bins),
            "P_backlog": {},
        }

    df = pd.read_csv(afc_path)

    origin_col = next(
        (c for c in ["origin", "origin_station", "origin_stop_id", "tap_in_stop_id"] if c in df.columns),
        None
    )
    destination_col = next(
        (c for c in ["destination", "destination_station", "destination_stop_id", "tap_out_stop_id"] if c in df.columns),
        None
    )
    time_col = next(
        (c for c in ["tap_in_time", "time", "timestamp", "entry_time"] if c in df.columns),
        None
    )
    count_col = next(
        (c for c in ["count", "passengers", "weight"] if c in df.columns),
        None
    )

    if origin_col is None or destination_col is None:
        raise ValueError("AFC input must include origin and destination columns.")

    p = {}
    backlog = {}
    phi_mass = {
        i: [0.0] * n_bins
        for i in range(1, n_stations + 1)
    }
    bin_seconds = horizon_seconds / float(n_bins)

    for row in df.itertuples(index=False):
        values = row._asdict()
        i = _station(values[origin_col], station_map)
        j = _station(values[destination_col], station_map)

        if i == j:
            continue

        if i > j:
            i, j = j, i

        count = float(values[count_col]) if count_col else 1.0
        seconds = _seconds(values[time_col]) if time_col else None

        if seconds is not None and seconds < horizon_start:
            backlog[(i, j)] = backlog.get((i, j), 0.0) + count
            continue

        p[(i, j)] = p.get((i, j), 0.0) + count

        if seconds is not None and horizon_start <= seconds < horizon_start + horizon_seconds:
            b = int((seconds - horizon_start) // bin_seconds)
            if 0 <= b < n_bins:
                phi_mass[i][b] += count

    return {
        "p": p,
        "phi": normalize_phi(phi_mass, n_stations, n_bins),
        "P_backlog": {
            (i, j): backlog.get((i, j), 0.0)
            for (i, j) in p
        } if backlog else default_backlog(p),
    }

