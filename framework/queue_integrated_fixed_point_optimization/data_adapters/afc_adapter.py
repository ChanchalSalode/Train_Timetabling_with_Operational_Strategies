import os

import pandas as pd

from backlog_loader import build_default_backlog
from config import H, S_D, T_WINDOW, _n_bins

from .schema import normalize_backlog, normalize_phi, uniform_phi


def _seconds(value):
    if pd.isna(value):
        return None

    if isinstance(value, (int, float)):
        return int(value)

    text = str(value)
    if ":" not in text:
        return None

    parts = text.split(":")
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
    horizon_start=H,
    horizon_seconds=T_WINDOW,
    n_bins=_n_bins
):
    station_map = station_map or {}

    if afc_path is None or not os.path.exists(afc_path):
        print("[WARNING] AFC file not found; using zero demand.")
        p = {}
        return {
            "p": p,
            "phi": uniform_phi(n_stations, n_bins),
            "P_backlog": {},
            "od_share_per_bin": None,
        }

    df = pd.read_csv(afc_path)

    origin_col = next(
        (c for c in ["origin", "origin_station", "origin_stop_id", "tap_in_stop_id"] if c in df.columns),
        None
    )
    dest_col = next(
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

    if origin_col is None or dest_col is None:
        raise ValueError("AFC input must include origin and destination columns.")

    p = {}
    phi_mass = {
        i: [0.0] * n_bins
        for i in range(1, n_stations + 1)
    }
    backlog = {}

    bin_seconds = horizon_seconds / float(n_bins)

    for row in df.itertuples(index=False):
        row_data = row._asdict()
        i = _station(row_data[origin_col], station_map)
        j = _station(row_data[dest_col], station_map)

        if i == j:
            continue

        count = float(row_data[count_col]) if count_col else 1.0

        if i > j:
            i, j = j, i

        t = _seconds(row_data[time_col]) if time_col else None

        if t is not None and t < horizon_start:
            backlog[(i, j)] = backlog.get((i, j), 0.0) + count
            continue

        p[(i, j)] = p.get((i, j), 0.0) + count

        if t is None:
            continue

        if horizon_start <= t < horizon_start + horizon_seconds:
            b = int((t - horizon_start) // bin_seconds)
            if 0 <= b < n_bins:
                phi_mass[i][b] += count

    phi = normalize_phi(phi_mass, n_stations, n_bins)

    if backlog:
        P_backlog = normalize_backlog(backlog, p)
    else:
        P_backlog = build_default_backlog(p)

    return {
        "p": p,
        "phi": phi,
        "P_backlog": P_backlog,
        "od_share_per_bin": None,
    }

