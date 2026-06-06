from dataclasses import dataclass, field


@dataclass
class ModelInputData:
    r: dict
    e: dict
    p: dict
    phi: dict = field(default_factory=dict)
    P_backlog: dict = field(default_factory=dict)
    station_map: dict = field(default_factory=dict)
    line_direction: dict = field(default_factory=dict)
    service_candidates: list = field(default_factory=list)


def uniform_phi(n_stations, n_bins):
    return {
        i: [1.0 / n_bins] * n_bins
        for i in range(1, n_stations + 1)
    }


def normalize_phi(phi_mass, n_stations, n_bins):
    phi = {}

    for i in range(1, n_stations + 1):
        vec = list(phi_mass.get(i, []))[:n_bins]
        vec += [0.0] * (n_bins - len(vec))
        total = sum(max(0.0, float(v)) for v in vec)

        if total > 0.0:
            phi[i] = [max(0.0, float(v)) / total for v in vec]
        else:
            phi[i] = [1.0 / n_bins] * n_bins

    return phi


def default_backlog(demand, pre_back_minutes=2.0, horizon_minutes=30.0):
    factor = pre_back_minutes / horizon_minutes

    return {
        (int(i), int(j)): factor * float(value)
        for (i, j), value in demand.items()
    }


def dense_od_matrix(demand, n_stations):
    return {
        (i, j): float(demand.get((i, j), 0.0))
        for i in range(1, n_stations + 1)
        for j in range(1, n_stations + 1)
    }
