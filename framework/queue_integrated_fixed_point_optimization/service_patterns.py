# service_patterns.py

from gurobipy import quicksum


def add_service_pattern_constraints(
    mdl,
    vars,
    sets,
    params,
    data
):
    """Add stopping pattern, route timing, and headway constraints."""

    tau = vars['tau']
    z = vars['z']
    x = vars['x']
    a1 = vars['a1']
    d = vars['d']
    h = vars['h']

    K_u = sets['K_u']
    K_d = sets['K_d']
    S_u = sets['S_u']
    S_d = sets['S_d']

    r = data['r']
    e = data['e']

    h_min = params['h_min']
    h_max = params['h_max']
    M = params['M']
    H = params['H']

    mdl.addConstrs(z[k, 1, 6] + z[k, 1, 8] + z[k, 3, 6] + z[k, 3, 8] == tau[k] for k in K_u)
    mdl.addConstrs(z[k, 9, 14] + z[k, 9, 16] + z[k, 11, 14] + z[k, 11, 16] == tau[k] for k in K_d)

    mdl.addConstrs(x[k, i] <= tau[k] for i in S_u for k in K_u)
    mdl.addConstrs(x[k, i] <= tau[k] for i in S_d for k in K_d)

    mdl.addConstrs(quicksum(x[k, i] for i in S_u if i > 6) <= M * (1 - z[k, 1, 6]) for k in K_u)
    mdl.addConstrs(x[k, i] >= z[k, 1, 6] for k in K_u for i in S_u if 1 <= i <= 6)
    mdl.addConstrs(x[k, i] >= z[k, 1, 8] for k in K_u for i in S_u if 1 <= i <= 8)
    mdl.addConstrs(quicksum(x[k, i] for i in S_u if i < 3) <= M * (1 - z[k, 3, 6]) for k in K_u)
    mdl.addConstrs(quicksum(x[k, i] for i in S_u if i > 6) <= M * (1 - z[k, 3, 6]) for k in K_u)
    mdl.addConstrs(x[k, i] >= z[k, 3, 6] for k in K_u for i in S_u if 3 <= i <= 6)
    mdl.addConstrs(quicksum(x[k, i] for i in S_u if i < 3) <= M * (1 - z[k, 3, 8]) for k in K_u)
    mdl.addConstrs(x[k, i] >= z[k, 3, 8] for k in K_u for i in S_u if 3 <= i <= 8)

    mdl.addConstrs(quicksum(x[k, i] for i in S_d if i > 14) <= M * (1 - z[k, 9, 14]) for k in K_d)
    mdl.addConstrs(x[k, i] >= z[k, 9, 14] for k in K_d for i in S_d if 9 <= i <= 14)
    mdl.addConstrs(x[k, i] >= z[k, 9, 16] for k in K_d for i in S_d if 9 <= i <= 16)
    mdl.addConstrs(quicksum(x[k, i] for i in S_d if i < 11) <= M * (1 - z[k, 11, 14]) for k in K_d)
    mdl.addConstrs(quicksum(x[k, i] for i in S_d if i > 14) <= M * (1 - z[k, 11, 14]) for k in K_d)
    mdl.addConstrs(x[k, i] >= z[k, 11, 14] for k in K_d for i in S_d if 11 <= i <= 14)
    mdl.addConstrs(quicksum(x[k, i] for i in S_d if i < 11) <= M * (1 - z[k, 11, 16]) for k in K_d)
    mdl.addConstrs(x[k, i] >= z[k, 11, 16] for k in K_d for i in S_d if 11 <= i <= 16)

    mdl.addConstrs(x[k - 1, i] + x[k, i] >= 1 for k in K_u for i in S_u if k != K_u[0])
    mdl.addConstrs(x[k - 1, i] + x[k, i] >= 1 for k in K_d for i in S_d if k != K_d[0])

    mdl.addConstr(d[K_u[0], S_u[0]] == H)
    mdl.addConstr(d[K_d[0], S_d[0]] == H)

    latest_departure = H + params['T_WINDOW']
    mdl.addConstrs(d[k, 1] <= latest_departure for k in K_u)
    mdl.addConstrs(d[k, 3] <= latest_departure for k in K_u)
    mdl.addConstrs(d[k, 9] <= latest_departure for k in K_d)
    mdl.addConstrs(d[k, 11] <= latest_departure for k in K_d)

    mdl.addConstrs(a1[k, i] - d[k, i - 1] == r[i - 1, i] for k in K_u for i in S_u if i != S_u[0])
    mdl.addConstrs(a1[k, i] - d[k, i - 1] == r[i - 1, i] for k in K_d for i in S_d if i != S_d[0])
    mdl.addConstrs(d[k, i] - a1[k, i] == e[i] for k in K_u for i in S_u)
    mdl.addConstrs(d[k, i] - a1[k, i] == e[i] for k in K_d for i in S_d)

    mdl.addConstrs(h_min * tau[k] <= h[k - 1, k] for k in K_u if k != K_u[0])
    mdl.addConstrs(h[k - 1, k] <= h_max * tau[k] for k in K_u if k != K_u[0])
    mdl.addConstrs(h_min * tau[k] <= h[k - 1, k] for k in K_d if k != K_d[0])
    mdl.addConstrs(h[k - 1, k] <= h_max * tau[k] for k in K_d if k != K_d[0])
    mdl.addConstrs(d[k, i] == d[k - 1, i] + h[k - 1, k] for k in K_u for i in S_u if k != K_u[0])
    mdl.addConstrs(d[k, i] == d[k - 1, i] + h[k - 1, k] for k in K_d for i in S_d if k != K_d[0])
