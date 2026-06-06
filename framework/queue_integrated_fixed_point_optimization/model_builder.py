# model_builder.py

from gurobipy import Model

from backlog_loader import get_backlog_matrix
from config import S_U, S_D
from demand_processing import (
    build_phi_from_excel,
    build_pwl_breaks_for_dir,
    load_od_share_per_bin
)
from network import create_network_sets, create_parameter_dict
from objective import set_objective
from passenger_constraints import add_passenger_constraints
from pwl_constraints import add_pwl_constraints
from rolling_stock_constraints import add_rolling_stock_constraints
from service_patterns import add_service_pattern_constraints
from variables import create_variables


def build_demand_profiles(p, params, stations_up=None, stations_down=None):
    """Build time-dependent PWL arrivals for both operating directions."""

    stations_up = stations_up or list(range(1, S_U + 1))
    stations_down = stations_down or list(range(S_U + 1, S_D + 1))

    phi = build_phi_from_excel(
        n_stations=S_D,
        n_bins=params['n_bins']
    )

    od_share_per_bin = load_od_share_per_bin(
        demand=p,
        n_stations=S_D,
        n_bins=params['n_bins']
    )

    pwl_up = build_pwl_breaks_for_dir(
        p,
        phi,
        stations_up,
        params['times'],
        params['n_bins'],
        od_share_per_bin=od_share_per_bin
    )

    pwl_down = build_pwl_breaks_for_dir(
        p,
        phi,
        stations_down,
        params['times'],
        params['n_bins'],
        od_share_per_bin=od_share_per_bin
    )

    return pwl_up, pwl_down


def _public_var_dict(vars, sets, params):
    """Return the variable dictionary expected by reporting and coupling code."""

    var_dicts = dict(vars)
    var_dicts.update({
        'K_u': sets['K_u'],
        'K_d': sets['K_d'],
        'S_u': sets['S_u'],
        'S_d': sets['S_d'],
        'B': sets['B'],
        'C': params['C'],
        'Aset': sets['A'],
        'D': sets['D']
    })
    return var_dicts


def build_trb_model(p, r, e, backlog=None, pwl_up=None, pwl_down=None):
    """Assemble the QFPO optimization model from focused modules."""

    mdl = Model('TRB_Santiago_16station')
    mdl.setParam('OutputFlag', 0)

    sets = create_network_sets()
    params = create_parameter_dict()
    vars = create_variables(mdl, sets, params)

    if backlog is None:
        backlog = get_backlog_matrix(p)

    if pwl_up is None or pwl_down is None:
        pwl_up, pwl_down = build_demand_profiles(
            p,
            params,
            stations_up=sets['S_u'],
            stations_down=sets['S_d']
        )

    data = {
        'r': r,
        'e': e
    }

    add_service_pattern_constraints(
        mdl,
        vars,
        sets,
        params,
        data
    )

    add_rolling_stock_constraints(
        mdl,
        vars,
        sets,
        params
    )

    add_pwl_constraints(
        mdl,
        vars,
        sets,
        params,
        pwl_up,
        pwl_down,
        backlog
    )

    add_passenger_constraints(
        mdl,
        vars,
        sets,
        params
    )

    var_dicts = _public_var_dict(vars, sets, params)

    set_objective(
        mdl,
        var_dicts
    )

    return mdl, var_dicts
