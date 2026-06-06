from dataclasses import dataclass, field


@dataclass
class ModelInputData:
    r: dict
    e: dict
    p: dict
    phi: dict
    P_backlog: dict
    od_share_per_bin: dict | None = None
    station_map: dict = field(default_factory=dict)
    line_direction: dict = field(default_factory=dict)
    service_candidates: list = field(default_factory=list)


def uniform_phi(n_stations, n_bins):
    return {
        i: [1.0 / n_bins] * n_bins
        for i in range(1, n_stations + 1)
    }


def normalize_phi(phi, n_stations, n_bins):
    out = {}

    for i in range(1, n_stations + 1):
        vec = list(phi.get(i, []))[:n_bins]
        vec += [0.0] * (n_bins - len(vec))
        total = sum(max(0.0, float(v)) for v in vec)

        if total > 0:
            out[i] = [max(0.0, float(v)) / total for v in vec]
        else:
            out[i] = [1.0 / n_bins] * n_bins

    return out


def normalize_backlog(backlog, demand):
    return {
        (int(i), int(j)): float(backlog.get((i, j), 0.0))
        for (i, j) in demand
    }

