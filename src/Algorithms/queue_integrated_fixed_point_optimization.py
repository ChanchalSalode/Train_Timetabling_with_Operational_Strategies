import os
import math
import traceback
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#Model Parameters
S_U = 8
S_D = 16
K_U = 6
K_D = 12
h_max = 360
h_min = 90
delta_min = 135
RS = 14
C = 250
M_time = 100000.0
M_flow = 5000.0
M = M_flow
H = 27000                 # 07:30:00 seconds 
T_WINDOW = 1800           # 30 minutes in seconds
BIN = 60                  # 1-minute bins
_times = list(range(0, T_WINDOW + 1, BIN))
_n_bins = len(_times) - 1
eps = 1e-9
PREPROC_HISTORY = []
#PWL-building and demand utilities
def build_phi_from_excel(path='time_dependent_16stations_30mint.xlsx'):
    try:
        df_phi = pd.read_excel(path, 'Arrival_profile_per_min')
        col_station = 'station' if 'station' in df_phi.columns else ('i' if 'i' in df_phi.columns else None)
        col_bin = 'bin_idx' if 'bin_idx' in df_phi.columns else ('minute_idx' if 'minute_idx' in df_phi.columns else None)
        col_phi = 'phi'
        if col_station is None or col_bin is None or col_phi not in df_phi.columns:
            raise ValueError("Arrival_profile_per_min columns not recognized")
        phi = {}
        for i in range(1, S_D + 1):
            vec = [0.0] * _n_bins
            sub = df_phi[df_phi[col_station] == i]
            for _, r0 in sub.iterrows():
                b = int(r0[col_bin])
                if 0 <= b < _n_bins:
                    vec[b] = float(r0[col_phi])
            s = sum(vec)
            phi[i] = [v / s if s > 0 else 1.0 / _n_bins for v in vec]
        return phi
    except Exception:
        return {i: [1.0 / _n_bins] * _n_bins for i in range(1, S_D + 1)}
def load_backlog_from_excel(path='TRB_16_stations.xlsx'):
    try:
        dfP = pd.read_excel(path, 'Backlog_at_H')
        P_backlog = {}
        for _, r0 in dfP.iterrows():
            P_backlog[(int(r0['i']), int(r0['j']))] = float(r0['P_ij'])
        return P_backlog
    except Exception:
        return None
def build_pwl_breaks_for_dir(p, phi, s_dir):
    pwl = {}
    for i in s_dir:
        J = [j for j in s_dir if j > i and (i, j) in p]
        if not J:
            continue
        row_sum = sum(p.get((i, j), 0.0) for j in J)
        if row_sum <= 0:
            for j in J:
                pwl[(i, j)] = (_times, [0.0] * (_n_bins + 1))
            continue
        bin_mass_i = [row_sum * phi[i][b] for b in range(_n_bins)]
        for j in J:
            pij = p.get((i, j), 0.0)
            if pij <= 0:
                pwl[(i, j)] = (_times, [0.0] * (_n_bins + 1))
                continue
            share = pij / row_sum if row_sum > 0 else 0.0
            mass_ij_per_bin = [bin_mass_i[b] * share for b in range(_n_bins)]
            cum = np.cumsum(mass_ij_per_bin)
            y_breaks = [0.0] + list(cum)
            y_breaks[-1] = pij  # snap final
            pwl[(i, j)] = (_times, y_breaks)
    return pwl
def pwl_to_minute_series(pwl_entry):
    x_breaks, y_breaks = pwl_entry
    masses = []
    for b in range(1, len(y_breaks)):
        masses.append(float(y_breaks[b] - y_breaks[b-1]))
    return masses
#Precomputed params orchestration (compute_precomputed_params + add constraints)
def compute_precomputed_params(all_stations, arrivals_by_station_od, services_df_or_csvpath, mu_board=0.5, dwell_base_sec=20.0, dwell_per_pax_sec=0.25,initial_backlogs=None):
    pre_w = {}
    pre_nb = {}
    pre_v = {}
    pre_nb_total = {}
    if isinstance(services_df_or_csvpath, str):
        services_all = pd.read_csv(services_df_or_csvpath)
    else:
        services_all = services_df_or_csvpath.copy()

    for station in all_stations:
        arrivals_by_od = arrivals_by_station_od.get(station, {})
        df_s = services_all[services_all['station_i'] == station]
        services_list = []
        for _, row in df_s.iterrows():
            services_list.append({
                'k': row['k'],
                'depart_time_min': float(row['depart_time_min']),
                'x_stop': int(row.get('x_stop', 1)),
                'capacity': float(row.get('capacity', 1e9)),
                'onboard_arrival': float(row.get('onboard_arrival', 0.0))
            })
        init_back = initial_backlogs.get(station) if initial_backlogs else None
        df_out = simulate_queues_station(arrivals_by_od, services_list,
                                         mu_board=mu_board, dwell_base_sec=dwell_base_sec,
                                         dwell_per_pax_sec=dwell_per_pax_sec, initial_backlog=init_back)
        for _, r in df_out.iterrows():
            k = r['k']; i = r['i']; j = r['j']
            pre_w[(k,i,j)] = float(r['w_before'])
            pre_nb[(k,i,j)] = float(r['nb_boarded'])
            pre_v[(k,i,j)] = float(r['v_left'])
            pre_nb_total[(k,i)] = float(r['nb_total'])
    return pre_w, pre_nb, pre_v, pre_nb_total
#Compute adaptive pre_penalty for iterative preprocessor
def compute_pre_penalty_from_model(model, preproc_history, alpha=0.1, min_penalty=1e-3, beta=1.0):
    try:
        Obj_base = abs(model.ObjVal) if getattr(model, 'ObjVal', None) is not None else None
    except Exception:
        Obj_base = None
    if Obj_base is None:
        Obj_base = 0.0
        if preproc_history:
            last = preproc_history[-1]
            pnb = last.get('pre_nb', {}) if isinstance(last, dict) else {}
            pnt = last.get('pre_nb_total', {}) if isinstance(last, dict) else {}
            Obj_base = float(max(1.0, sum(pnb.values()) + sum(pnt.values())))
        else:
            Obj_base = 1.0
    if preproc_history:
        last = preproc_history[-1]
        N_slack_est = max(1, len(last.get('pre_nb_total', {}) if isinstance(last, dict) else {}))
    else:
        N_slack_est = 1
    per_slack = alpha * float(Obj_base) / float(N_slack_est)
    max_total = beta * float(max(1.0, Obj_base))
    per_slack = min(per_slack, max_total / float(N_slack_est))
    per_slack = max(per_slack, float(min_penalty))
    return float(per_slack)

from gurobipy import quicksum
def add_precomputed_constraints_to_model(model, pre_w, pre_nb, pre_nb_total,var_maps, mode='ub', pre_penalty=10000.0):
    try:
        to_remove = [c for c in model.getConstrs() if c.ConstrName and c.ConstrName.startswith('PRE_')]
        if to_remove:
            model.remove(to_remove)
            model.update()
    except Exception:
        pass
    try:
        prev_slacks = [v for v in model.getVars() if v.VarName and v.VarName.startswith('PRE_s_')]
        if prev_slacks:
            model.remove(prev_slacks)
            model.update()
    except Exception:
        pass
    created_slacks = []
    nb_dict = var_maps.get('nb_var', {})
    nbtot_dict = var_maps.get('nb_total_var', {})
    wdict = var_maps.get('w_var', {})
    for (k,i,j), val in pre_nb.items():
        if (k,i,j) in nb_dict:
            var = nb_dict[(k,i,j)]
            if mode == 'ub':
                cname = f"PRE_nb_ub_{k}_{i}_{j}"
                model.addConstr(var <= val + 1e-6, name=cname)
            elif mode == 'fix':
                cname = f"PRE_nb_eq_{k}_{i}_{j}"
                model.addConstr(var == val, name=cname)
            elif mode == 'soft':
                s = model.addVar(lb=0.0, name=f"PRE_s_nb_{k}_{i}_{j}")
                created_slacks.append(s)
                cname = f"PRE_nb_soft_{k}_{i}_{j}"
                model.addConstr(var <= val + s, name=cname)
    for (k,i), val in pre_nb_total.items():
        if (k,i) in nbtot_dict:
            var = nbtot_dict[(k,i)]
            if mode == 'ub':
                cname = f"PRE_nbtotal_ub_{k}_{i}"
                model.addConstr(var <= val + 1e-6, name=cname)
            elif mode == 'fix':
                cname = f"PRE_nbtotal_eq_{k}_{i}"
                model.addConstr(var == val, name=cname)
            elif mode == 'soft':
                s = model.addVar(lb=0.0, name=f"PRE_s_nbtotal_{k}_{i}")
                created_slacks.append(s)
                cname = f"PRE_nbtotal_soft_{k}_{i}"
                model.addConstr(var <= val + s, name=cname)
    model.update()
    if created_slacks:
        try:
            cur_obj = model.getObjective()
        except Exception:
            cur_obj = None
        penalty_expr = pre_penalty * quicksum(created_slacks)
        if cur_obj is not None:
            model.setObjective(cur_obj - penalty_expr)
        else:
            model.setObjective(-penalty_expr)
        model.update()
    return created_slacks
#Iterative coupling loop
from gurobipy import GRB
def run_iterative_with_preprocessor(model, all_stations, arrivals_by_station_od, services_csv_or_df,var_maps, max_iters=3, tol=1e-3, mu_board=0.5, dwell_base_sec=20.0,dwell_per_pax_sec=0.25, initial_backlogs=None, mode='soft', pre_penalty=10000.0):
    prev_pre_nb = None
    services_src = services_csv_or_df
    for itr in range(1, max_iters + 1):
        print(f"[Preproc Iter {itr}] Computing precomputed params...")
        pre_w, pre_nb, pre_v, pre_nb_total = compute_precomputed_params(all_stations, arrivals_by_station_od,services_src,mu_board=mu_board, dwell_base_sec=dwell_base_sec,dwell_per_pax_sec=dwell_per_pax_sec,initial_backlogs=initial_backlogs)
        try:
            PREPROC_HISTORY.append({'iter': itr,'pre_w': pre_w.copy(),'pre_nb': pre_nb.copy(),'pre_v': pre_v.copy(),'pre_nb_total': pre_nb_total.copy()})
        except Exception:
            PREPROC_HISTORY.append({'iter': itr})
        print(f"[Preproc Iter {itr}] Adding constraints to model...")
        add_precomputed_constraints_to_model(model, pre_w, pre_nb, pre_nb_total, var_maps, mode=mode)
        print(f"[Preproc Iter {itr}] Solving MIP...") 
        model.optimize()
        try:
            n_pre = sum(1 for c in model.getConstrs() if c.ConstrName and c.ConstrName.startswith('PRE_'))
            n_slack_vars = sum(1 for v in model.getVars() if v.VarName and v.VarName.startswith('PRE_s_'))
            penalty_contrib = 0.0
            for v in model.getVars():
                if v.VarName and v.VarName.startswith('PRE_s_'):
                    xv = getattr(v, 'X', None)
                    if xv is None:
                        xv = getattr(v, 'x', 0.0)
                    try:
                        penalty_contrib += float(xv)
                    except Exception:
                        pass
            penalty_contrib *= float(pre_penalty)
            print(f"[Preproc Diagnostics] PRE constraints={n_pre}, PRE slacks={n_slack_vars}, pre_penalty={pre_penalty}, penalty_contrib={penalty_contrib}")
        except Exception as _diag_e:
            print('Failed to compute preproc diagnostics:', _diag_e)
        status = model.Status
        if status == GRB.INFEASIBLE:
            print(f"[Preproc Iter {itr}] Model became INFEASIBLE after PRE constraints. Computing IIS and removing PRE constraints...")
            try:
                iis_name = f"mdl_preproc_iter{itr}.ilp"
                model.computeIIS()
                model.write(iis_name)
                print(f"[Preproc Iter {itr}] IIS written to {iis_name}.")
            except Exception as e:
                print(f"[Preproc Iter {itr}] Failed to compute/write IIS: {e}")
            try:
                to_remove = [c for c in model.getConstrs() if c.ConstrName and c.ConstrName.startswith('PRE_')]
                if to_remove:
                    model.remove(to_remove)
                    model.update()
            except Exception:
                pass
            print(f"[Preproc Iter {itr}] PRE_ constraints removed. Stopping preprocessor to recover baseline model.")
            break

        if status not in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL, GRB.INTERRUPTED):
            print(f"[Preproc Iter {itr}] Solver returned status {status}. Stopping iterative preprocessor.")
            break
        try:
            update_services_from_solution(model, services_src)
        except Exception as e:
            print("NOTE: update_services_from_solution not implemented or failed. Continue with unchanged services.")
            traceback.print_exc()
        # convergence check on pre_nb averages
        if prev_pre_nb is not None:
            keys = set(pre_nb.keys()).intersection(prev_pre_nb.keys())
            if len(keys) > 0:
                diff = sum(abs(pre_nb[k] - prev_pre_nb[k]) for k in keys) / len(keys)
                print(f"[Preproc Iter {itr}] avg change in pre_nb = {diff:.6f}")
                if diff < tol:
                    print(f"[Preproc Iter {itr}] Converged (tol={tol}). Stopping iter.")
                    break
        prev_pre_nb = pre_nb
    print("Iterative preprocessor+solve finished.")
    return
#Gurobi model builder
from gurobipy import Model
def build_trb_model(p, r, e):
    mdl = Model('TRB_Santiago_16station')
    mdl.setParam('OutputFlag', 0)
    #Index Sets
    S_d = [i for i in range(9, S_D + 1)]
    S_dd = [i for i in range(1, S_D + 1)]
    K_d = [i for i in range(7, K_D + 1)]
    K_dd = [i for i in range(1, K_D + 1)]
    K_u = [i for i in range(1, K_U + 1)]
    S_u = [i for i in range(1, S_U + 1)]
    B = [(k, m, n) for k in K_dd for m in S_dd for n in S_dd]
    A = [(k, l, m) for k in K_dd for l in K_dd for m in S_dd]
    C1 = [(k, i) for k in K_dd for i in S_dd]
    D = [(k, l) for k in K_dd for l in K_dd]

    #Decision Variables
    tau = mdl.addVars(K_dd, vtype=GRB.BINARY, name='tau')
    z = mdl.addVars(B, vtype=GRB.BINARY, name='z')
    y = mdl.addVars(A, vtype=GRB.BINARY, name='y')
    x = mdl.addVars(C1, vtype=GRB.BINARY, name='x')
    a1 = mdl.addVars(C1, vtype=GRB.CONTINUOUS, name='a1')
    d = mdl.addVars(C1, vtype=GRB.CONTINUOUS, name='d')
    h = mdl.addVars(D, vtype=GRB.CONTINUOUS, name='h')
    alpha_up1 = mdl.addVars(K_u, vtype=GRB.BINARY, name='alpha_u1')
    alpha_dn2 = mdl.addVars(K_d, vtype=GRB.BINARY, name='alpha_d2')
    alpha_up3 = mdl.addVars(K_u, vtype=GRB.BINARY, name='alpha_u3')
    alpha_dn4 = mdl.addVars(K_d, vtype=GRB.BINARY, name='alpha_d4')
    beta_up2 = mdl.addVars(K_u, vtype=GRB.BINARY, name='beta_u2')
    beta_dn1 = mdl.addVars(K_d, vtype=GRB.BINARY, name='beta_d1')
    beta_up4 = mdl.addVars(K_u, vtype=GRB.BINARY, name='beta_u4')
    beta_dn3 = mdl.addVars(K_d, vtype=GRB.BINARY, name='beta_d3')
    RS1 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS1')
    RS2 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS2')
    RS3 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS3')
    RS4 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS4')
    w = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='w')
    w_b = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='w_b')
    w_b1 = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='w_b1')
    n_b = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='n_b')
    n_b1 = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='n_b1')
    n_a = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='n_a')
    n1 = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, ub=C, name='n1')
    v = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='v')
    sai = mdl.addVars(C1, vtype=GRB.BINARY, name='sai')
    phi = {}
    try:
        df_phi = pd.read_excel('time_dependent_16stations_30mint.xlsx', 'Arrival_profile_per_min')
        col_station = 'station' if 'station' in df_phi.columns else ('i' if 'i' in df_phi.columns else None)
        col_bin = 'bin_idx' if 'bin_idx' in df_phi.columns else ('minute_idx' if 'minute_idx' in df_phi.columns else None)
        col_phi = 'phi'
        if col_station is None or col_bin is None or col_phi not in df_phi.columns:
            raise ValueError("Arrival_profile_per_min columns not recognized")
        for i in range(1, S_D + 1):
            vec = [0.0] * _n_bins
            sub = df_phi[df_phi[col_station] == i]
            for _, r0 in sub.iterrows():
                b = int(r0[col_bin])
                if 0 <= b < _n_bins:
                    vec[b] = float(r0[col_phi])
            s = sum(vec)
            phi[i] = [v / s if s > 0 else 1.0 / _n_bins for v in vec]
    except Exception:
        phi = {i: [1.0 / _n_bins] * _n_bins for i in range(1, S_D + 1)}
    P_backlog = {}
    try:
        dfP = pd.read_excel('TRB_16_stations.xlsx', 'Backlog_at_H')
        for _, r0 in dfP.iterrows():
            P_backlog[(int(r0['i']), int(r0['j']))] = float(r0['P_ij'])
    except Exception:
        PRE_BACK_MIN = 2.0
        factor = PRE_BACK_MIN / 30.0
        P_backlog = {(int(i), int(j)): factor * float(val) for (i, j), val
    in p.items()}
    s_share_per_bin = None
    try:
        df_sh = pd.read_excel('time_dependent_16stations_30mint.xlsx', 'OD_share_per_bin')
        tmp = {}
        for (ii, jj), _val in p.items():
            tmp[(int(ii), int(jj))] = [0.0] * _n_bins
        for _, r0 in df_sh.iterrows():
            ii = int(r0['i']); jj = int(r0['j']); bb = int(r0['bin_idx']); val = float(r0['share'])
            if (ii, jj) in tmp and 0 <= bb < _n_bins:
                tmp[(ii, jj)][bb] = max(0.0, val)
        for ii in range(1, S_D + 1):
            Js_all = [jj for (i2, jj) in p.keys() if i2 == ii and jj > ii]
            for bb in range(_n_bins):
                ssum = sum(tmp[(ii, jj)][bb] for jj in Js_all) if Js_all else 0.0
                if ssum > 0:
                    for jj in Js_all:
                        tmp[(ii, jj)][bb] /= ssum
                else:
                    row_sum_all = sum(p.get((ii, jj), 0.0) for jj in Js_all)
                    for jj in Js_all:
                        tmp[(ii, jj)][bb] = (p.get((ii, jj), 0.0) /
    row_sum_all) if row_sum_all > 0 else 0.0
        s_share_per_bin = tmp
    except Exception:
        s_share_per_bin = None
    def _build_pwl_breaks_for_dir(S_dir):
        pwl = {}
        for i in S_dir:
            J = [j for j in S_dir if j > i and (i, j) in p]
            if not J:
                continue
            row_sum = sum(p.get((i, j), 0.0) for j in J)
            if row_sum <= 0:
                for j in J:
                    pwl[(i, j)] = (_times, [0.0] * (_n_bins + 1))
                continue
            bin_mass_i = [row_sum * phi[i][b] for b in range(_n_bins)]

            for j in J:
                pij = p.get((i, j), 0.0)
                if pij <= 0:
                    pwl[(i, j)] = (_times, [0.0] * (_n_bins + 1))
                    continue
                if s_share_per_bin is not None and (i, j) in s_share_per_bin:
                    s_vec = s_share_per_bin[(i, j)]
                    raw = [bin_mass_i[b] * max(0.0, float(s_vec[b])) for b in range(_n_bins)]
                    raw_sum = sum(raw)
                    if raw_sum > 0.0:
                        mass_ij_per_bin = [raw[b] * (pij / raw_sum) for b in range(_n_bins)]
                    else:
                        mass_ij_per_bin = [pij / float(_n_bins)] * _n_bins
                else:
                    share = pij / row_sum if row_sum > 0.0 else 0.0
                    mass_ij_per_bin = [bin_mass_i[b] * share for b in range(_n_bins)]
                cum = np.cumsum(mass_ij_per_bin)
                y_breaks = [0.0] + list(cum)
                y_breaks[-1] = pij
                pwl[(i, j)] = (_times, y_breaks)
        return pwl
    _pwl_up = _build_pwl_breaks_for_dir(S_u)
    _pwl_down = _build_pwl_breaks_for_dir(S_d)
    t_up = mdl.addVars(K_u, S_u, lb=0.0, ub=T_WINDOW, name="t_up")
    t_down = mdl.addVars(K_d, S_d, lb=0.0, ub=T_WINDOW, name="t_down")
    A_up = {}
    A_dn = {}
    #Operation Planning Constraints
    mdl.addConstrs(t_up[k, i] == d[k, i] - H for k in K_u for i in S_u)
    mdl.addConstrs(t_down[k, i] == d[k, i] - H for k in K_d for i in S_d)
    mdl.addConstrs(z[k,1,6] + z[k,1,8] + z[k,3,6]+ z[k,3,8]  == tau[k] for k in K_u)
    mdl.addConstrs(z[k,9,14] + z[k,9,16] + z[k,11,14] + z[k,11,16]  == tau[k] for k in K_d)
    mdl.addConstrs(x[k,i] <= tau[k] for i in S_u for k in K_u)
    mdl.addConstrs(x[k,i] <= tau[k] for i in S_d for k in K_d)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M*(1-z[k,1,6])for k in K_u)
    mdl.addConstrs(x[k,i] >= z[k,1,6] for k in K_u for i in S_u if i >= 1 if i <= 6)
    mdl.addConstrs(x[k,i] >= z[k,1,8] for k in K_u for i in S_u if i >= 1 if i <= 8)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M*(1-z[k,3,6])for k in K_u)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M*(1-z[k,3,6])for k in K_u)
    mdl.addConstrs(x[k,i] >= z[k,3,6] for k in K_u for i in S_u if i >= 3 if i <= 6)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M*(1-z[k,3,8])for k in K_u)
    mdl.addConstrs(x[k,i] >= z[k,3,8] for k in K_u for i in S_u if i >= 3 if i <= 8)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M*(1-z[k,9,14])for k in K_d)
    mdl.addConstrs(x[k,i] >= z[k,9,14] for k in K_d for i in S_d if i >= 9 if i <= 14)
    mdl.addConstrs(x[k,i] >= z[k,9,16] for k in K_d for i in S_d if i >= 9 if i <= 16)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M*(1-z[k,11,14])for k in K_d)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M*(1-z[k,11,14])for k in K_d)
    mdl.addConstrs(x[k,i] >= z[k,11,14] for k in K_d for i in S_d if i >= 11 if i <= 14)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M*(1-z[k,11,16])for k in K_d)
    mdl.addConstrs(x[k,i] >= z[k,11,16] for k in K_d for i in S_d if i >= 11 if i <= 16)
    mdl.addConstrs(x[k-1,i] + x[k,i] >= 1 for k in K_u for i in S_u if k != 1)
    mdl.addConstrs(x[k-1,i] + x[k,i] >= 1 for k in K_d for i in S_d if k != 7)
    mdl.addConstr(d[1,1] == 27000)      
    mdl.addConstr(d[7,9] == 27000)     
    mdl.addConstrs(d[k,1] <= 28800 for k in K_u)
    mdl.addConstrs(d[k,3] <= 28800 for k in K_u)
    mdl.addConstrs(d[k,9] <= 28800 for k in K_d)
    mdl.addConstrs(d[k,11] <= 28800 for k in K_d)
    mdl.addConstrs(a1[k,i] - d[k,i-1] == r[i-1,i] for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(a1[k,i] - d[k,i-1] == r[i-1,i] for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(d[k,i] - a1[k,i] == e[i] for k in K_u for i in S_u)
    mdl.addConstrs(d[k,i] - a1[k,i] == e[i] for k in K_d for i in S_d)
    mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_u if k!= 1)
    mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_u if k!= 1)
    mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_d if k!= 7)
    mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_d if k!= 7)
    mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_u for i in S_u if k != 1)
    mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_d for i in S_d if k != 7)
    mdl.addConstrs(2*y[k,l,6] <= z[k,1,6] +z[l,11,16] + z[k,3,6] + z[l,11,14] for k in K_u for l in K_d)
    mdl.addConstrs(2*y[k,l,8] <= z[k,1,8] + z[l,9,16] +z[k,3,8] + z[l,9,14] for k in K_u for l in K_d)
    mdl.addConstrs(2*y[k,l,14] <= z[k,9,14] + z[l,3,8]+z[k,11,14] + z[l,3,6] for k in K_d for l in K_u)
    mdl.addConstrs(2*y[k,l,16] <= z[k,9,16] + z[l,1,8]+z[k,11,16] + z[l,1,6] for k in K_d for l in K_u)
    mdl.addConstrs(a1[l,9] - d[k,8] >= delta_min + M_time*(y[k,l,8]-1) for k in K_u for l in K_d)
    mdl.addConstrs(a1[l,11] - d[k,6] >= delta_min + M_time*(y[k,l,6]-1) for k in K_u for l in K_d)
    mdl.addConstrs(a1[l,1] - d[k,16] >= delta_min + M_time*(y[k,l,16]-1) for k in K_d for l in K_u)
    mdl.addConstrs(a1[l,3] - d[k,14] >= delta_min + M_time*(y[k,l,14]-1) for k in K_d for l in K_u)
    mdl.addConstrs(quicksum(y[l,k,16]for l in K_d)+quicksum(y[l,k,14]for l in K_d) + alpha_up1[k] + alpha_up3[k] == tau[k] for k in K_u)
    mdl.addConstrs(quicksum(y[l,k,6]for l in K_u)+quicksum(y[l,k,8]for l in K_u) + alpha_dn2[k] + alpha_dn4[k] == tau[k] for k in K_d)
    mdl.addConstrs(quicksum(y[k,l,6]for l in K_d)+quicksum(y[k,l,8]for l in K_d) + beta_up2[k]+ beta_up4[k] == tau[k] for k in K_u)
    mdl.addConstrs(quicksum(y[k,l,16]for l in K_u)+quicksum(y[k,l,14]for l in K_u) + beta_dn1[k]+ beta_dn3[k] == tau[k] for k in K_d)
    mdl.addConstr(quicksum(alpha_up1[k] for k in K_u) + quicksum(alpha_dn2[k] for k in K_d)+quicksum(alpha_up3[k] for k in K_u) + quicksum(alpha_dn4[k] for k in K_d)<= RS)
    mdl.addConstr(quicksum(alpha_up1[k] for k in K_u) <= RS1)
    mdl.addConstr(quicksum(alpha_dn2[k] for k in K_d) <= RS2)
    mdl.addConstr(quicksum(alpha_up3[k] for k in K_u) <= RS3)
    mdl.addConstr(quicksum(alpha_dn4[k] for k in K_d) <= RS4)
    mdl.addConstrs(alpha_up1[k] <= z[k,1,8] + z[k,1,6] for k in K_u)
    mdl.addConstrs(alpha_dn2[k] <= z[k,9,16] + z[k,9,14] for k in K_d)
    mdl.addConstrs(beta_dn1[k] <= z[k,9,16]+ z[k,11,16]for k in K_d)
    mdl.addConstrs(beta_up2[k] <= z[k,1,8] + z[k,3,8] for k in K_u)
    mdl.addConstrs(alpha_up3[k] <= z[k,3,8] + z[k,3,6] for k in K_u)
    mdl.addConstrs(alpha_dn4[k] <= z[k,11,16] + z[k,11,14] for k in K_d)
    mdl.addConstrs(beta_dn3[k] <= z[k,9,14]+z[k,11,14] for k in K_d)
    mdl.addConstrs(beta_up4[k] <= z[k,1,6] + z[k,3,6] for k in K_u)
    #Passenger Demand Constraints (Time Dependent)
    for k in K_u:
        for i in S_u:
            for j in S_u:
                if i < j and (i, j) in _pwl_up:
                    x_breaks, y_breaks = _pwl_up[(i, j)]
                    A_ki = mdl.addVar(lb=0.0, name=f"A_up_{k}_{i}_{j}")
                    mdl.addGenConstrPWL(t_up[k, i], A_ki, x_breaks, y_breaks, name=f"PWL_up_{k}_{i}_{j}")
                    A_up[(k, i, j)] = A_ki
                    if k == 1:
                        base = P_backlog.get((i, j), 0.0)
                        mdl.addConstr(w[k, i, j] == base + A_ki)
                    else:
                        A_prev = mdl.addVar(lb=0.0, name=f"A_up_prev_{k}_{i}_{j}")
                        mdl.addGenConstrPWL(t_up[k - 1, i], A_prev, x_breaks, y_breaks, name=f"PWL_up_prev_{k}_{i}_{j}")
                        mdl.addConstr(w[k, i, j] == v[k - 1, i, j] + (A_ki - A_prev))
    for k in K_d:
        for i in S_d:
            for j in S_d:
                if i < j and (i, j) in _pwl_down:
                    x_breaks, y_breaks = _pwl_down[(i, j)]
                    A_ki = mdl.addVar(lb=0.0, name=f"A_dn_{k}_{i}_{j}")
                    mdl.addGenConstrPWL(t_down[k, i], A_ki, x_breaks, y_breaks, name=f"PWL_dn_{k}_{i}_{j}")
                    A_dn[(k, i, j)] = A_ki
                    if k == 7:
                        base = P_backlog.get((i, j), 0.0)
                        mdl.addConstr(w[k, i, j] == base + A_ki)
                    else:
                        A_prev = mdl.addVar(lb=0.0, name=f"A_dn_prev_{k}_{i}_{j}")
                        mdl.addGenConstrPWL(t_down[k - 1, i], A_prev, x_breaks, y_breaks, name=f"PWL_dn_prev_{k}_{i}_{j}")
                        mdl.addConstr(w[k, i, j] == v[k - 1, i, j] + (A_ki - A_prev))
    mdl.addConstrs(w_b[k,i] == quicksum(w_b1[k,i,j] for j in S_u if i < j) for k in K_u for i in S_u)
    mdl.addConstrs(w_b[k,i] == quicksum(w_b1[k,i,j] for j in S_d if i < j) for k in K_d for i in S_d)
    mdl.addConstrs(n_a[k,i] == quicksum(n_b1[k,j,i] for j in S_u if j < i) for k in K_u for i in S_u)
    mdl.addConstrs(n_a[k,i] == quicksum(n_b1[k,j,i] for j in S_d if j < i) for k in K_d for i in S_d)
    mdl.addConstrs(n_b[k,i] - n_b1[k,i,j] == w_b[k,i] - w_b1[k,i,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(n_b[k,i] - n_b1[k,i,j] == w_b[k,i] - w_b1[k,i,j] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= 0 for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= w[k,i,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,i] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M*(2-x[k,i]-x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= 0 for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= w[k,i,j] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,i] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,j] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M*(2-x[k,i]-x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_u for i in S_u)
    mdl.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_u for i in S_u if i!= 1)
    mdl.addConstrs(n_b[k,i] >= w_b[k,i] - M*(1-sai[k,i]) for k in K_u for i in S_u)
    mdl.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M*sai[k,i] for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_d for i in S_d)
    mdl.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_d for i in S_d if i!= 9)
    mdl.addConstrs(n_b[k,i] >= w_b[k,i] - M*(1-sai[k,i]) for k in K_d for i in S_d)
    mdl.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M*sai[k,i] for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(n1[k,6] <= C * (1 - z[k,1,6]) for k in K_u)
    mdl.addConstrs(n1[k,6] <= C * (1 - z[k,3,6]) for k in K_u)
    mdl.addConstrs(n1[k,8] <= C * (1 - z[k,1,8]) for k in K_u)
    mdl.addConstrs(n1[k,8] <= C * (1 - z[k,3,8]) for k in K_u)
    mdl.addConstrs(n1[k,14] <= C * (1 - z[k,9,14]) for k in K_d)
    mdl.addConstrs(n1[k,14] <= C * (1 - z[k,11,14]) for k in K_d)
    mdl.addConstrs(n1[k,16] <= C * (1 - z[k,9,16]) for k in K_d)
    mdl.addConstrs(n1[k,16] <= C * (1 - z[k,11,16]) for k in K_d)
    mdl.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(n1[k,1] == n_b[k,1] for k in K_u)
    mdl.addConstrs(n1[k,9] == n_b[k,9] for k in K_d)
    mdl.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M*(2 - x[k,i] - x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M*(2 - x[k,i] - x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,i] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M*(2 - x[k,i] - x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M*(2 - x[k,i] - x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,i] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,j] for k in K_d for i in S_d for j in S_d if i < j)
    #Objective Function
    mdl.setObjective(quicksum(y[l,k,6] for l in K_u for k in K_d) + quicksum(y[l,k,8] for l in K_u for k in K_d) +
                     quicksum(y[l,k,14] for l in K_d for k in K_u) + quicksum(y[l,k,16] for l in K_d for k in K_u), GRB.MAXIMIZE)
    #Return model and variables
    var_dicts = {
        'tau': tau, 'z': z, 'y': y, 'x': x, 'a1': a1, 'd': d, 'h': h,
        'w': w, 'w_b': w_b, 'w_b1': w_b1, 'n_b': n_b, 'n_b1': n_b1, 'n_a': n_a,
        'n1': n1, 'v': v, 'sai': sai, 't_up': t_up, 't_down': t_down,
        'A_up': A_up, 'A_dn': A_dn, 'K_u': K_u, 'K_d': K_d, 'S_u': S_u, 'S_d': S_d,
        'B': B, 'C': C, 'Aset': A, 'D': D
    }
    return mdl, var_dicts
#Utility helpers
def integrate_bin(arrivals, t1_min, t2_min):
    if t2_min <= t1_min:
        return 0.0
    total = 0.0
    start_min = int(math.floor(t1_min))
    end_min = int(math.floor(t2_min))
    if start_min == end_min:
        frac = (t2_min - t1_min)
        if 0 <= start_min < len(arrivals):
            return arrivals[start_min] * frac
        return 0.0
    first_frac = (start_min + 1 - t1_min)
    if 0 <= start_min < len(arrivals):
        total += arrivals[start_min] * first_frac
    for m in range(start_min + 1, end_min):
        if 0 <= m < len(arrivals):
            total += arrivals[m]
    last_frac = (t2_min - end_min)
    if 0 <= end_min < len(arrivals):
        total += arrivals[end_min] * last_frac
    return total
def var_value_safe(var, default=0.0):
    try:
        if var is None:
            return float(default)
        v = getattr(var, 'X', None)
        if v is None:
            v = getattr(var, 'x', None)
        return float(v) if v is not None else float(default)
    except Exception:
        return float(default)

#Fluid-queue simulator per station
def simulate_queues_station(arrivals_by_od, services_at_station,mu_board=0.5,dwell_base_sec=20.0,dwell_per_pax_sec=0.25,initial_backlog=None):
    mu_board_per_min = mu_board * 60.0
    Q = defaultdict(float)
    if initial_backlog:
        for od, val in initial_backlog.items():
            try:
                Q[od] = float(val)
            except Exception:
                continue
    services_sorted = sorted(services_at_station, key=lambda s: s['depart_time_min'])
    last_time = 0.0
    records = []
    # od keys to track
    od_keys = set(arrivals_by_od.keys())
    if initial_backlog:
        od_keys.update(initial_backlog.keys())
    for s in services_sorted:
        k = s['k']; depart_time = float(s['depart_time_min']); x_stop = int(s.get('x_stop',1))
        capacity = float(s.get('capacity', 1e9)); onboard_arrival = float(s.get('onboard_arrival', 0.0))
        A_by_od = {}
        for od in od_keys:
            series = arrivals_by_od.get(od, None)
            if series is None:
                A = 0.0
            else:
                try:
                    A = integrate_bin(series, last_time, depart_time)
                except Exception:
                    A = 0.0
            A_by_od[od] = A
            Q[od] = Q.get(od, 0.0) + A
        Q_before = dict(Q)
        nb_by_od = {od: 0.0 for od in od_keys}
        v_by_od = {od: Q_before.get(od, 0.0) for od in od_keys}
        nb_total = 0.0
        if x_stop == 0:
            pass
        else:
            total_queue = sum(Q_before.get(od, 0.0) for od in od_keys)
            if total_queue <= 0.0:
                pass
            else:
                dwell_min = dwell_base_sec / 60.0
                approx_capacity_by_time = mu_board_per_min * dwell_min
                approx_capacity = min(max(0.0, capacity - onboard_arrival), approx_capacity_by_time)
                pred_board = min(total_queue, approx_capacity)
                dwell_sec = dwell_base_sec + dwell_per_pax_sec * pred_board
                dwell_min = dwell_sec / 60.0
                capacity_by_time = mu_board_per_min * dwell_min
                C = max(0.0, min(max(0.0, capacity - onboard_arrival), capacity_by_time))
                positive_queue_od = [od for od in od_keys if Q_before.get(od, 0.0) > 0]
                total_pos_queue = sum(Q_before.get(od, 0.0) for od in positive_queue_od)
                if total_pos_queue <= 0.0:
                    pass
                else:
                    for od in od_keys:
                        qval = Q_before.get(od, 0.0)
                        if qval <= 0.0:
                            nb_by_od[od] = 0.0
                        else:
                            alloc = C * (qval / total_pos_queue)
                            nb_by_od[od] = min(qval, alloc)
                    nb_total = sum(nb_by_od.values())
                    leftover = C - nb_total
                    if leftover > 0.5:
                        rem_list = sorted(positive_queue_od, key=lambda od: (Q_before.get(od,0.0) - nb_by_od.get(od,0.0)), reverse=True)
                        for od in rem_list:
                            if leftover <= 0.5:
                                break
                            can_take = Q_before.get(od,0.0) - nb_by_od.get(od,0.0)
                            give = min(can_take, leftover)
                            nb_by_od[od] += give
                            leftover -= give
                        nb_total = sum(nb_by_od.values())
                    v_by_od = {od: Q_before.get(od,0.0) - nb_by_od.get(od,0.0) for od in od_keys}
        for od in sorted(Q_before.keys()):
            i, j = od
            records.append({
                'k': k, 'i': i, 'j': j, 'depart_time_min': depart_time,
                'w_before': float(Q_before.get(od, 0.0)), 'A_arrivals': float(A_by_od.get(od, 0.0)),
                'nb_boarded': float(nb_by_od.get(od, 0.0)), 'v_left': float(v_by_od.get(od, 0.0)),
                'nb_total': float(nb_total)
            })
        for od in Q_before.keys():
            Q[od] = v_by_od.get(od, Q_before.get(od,0.0))
        last_time = depart_time

    df = pd.DataFrame.from_records(records)
    return df
#Metrics / diagnostics (single merged safe implementation)
def collect_performance_metrics_from_model(model, var_dicts):
    try:
        tau = var_dicts.get('tau', {})
        n1 = var_dicts.get('n1', {})
        w = var_dicts.get('w', {})
        v = var_dicts.get('v', {})
        w_b = var_dicts.get('w_b', {})
        h = var_dicts.get('h', {})
        x = var_dicts.get('x', {})
        B = var_dicts.get('B', [])
    except Exception:
        return {'max_plf':0.0,'overloaded':[],'cong_seconds':0.0,'cong_minutes':0.0,'left_behind_rate':0.0,'dcr_by_ki':{},'dcr_details':{},'peak_dcr':0.0}
    cap = float(C)
    plf_by_ki = {}
    for (k,i) in list(n1.keys()):
        tau_val = var_value_safe(tau.get(k) if isinstance(tau, dict) else None, 1.0)
        if tau_val > 0.5:
            plf_by_ki[(k,i)] = var_value_safe(n1.get((k,i)), 0.0) / (cap if cap>eps else 1.0)
    max_plf = max(plf_by_ki.values()) if plf_by_ki else 0.0
    overloaded = [(k,i,round(p,4)) for (k,i),p in plf_by_ki.items() if p > 0.85]
    cong_seconds = 0.0
    services_with_cong = set([k for (k,i,p) in overloaded])
    for k in services_with_cong:
        idx = (k-1, k)
        if idx in h:
            cong_seconds += var_value_safe(h.get(idx))
        else:
            cong_seconds += float(h_min)
    cong_minutes = cong_seconds / 60.0
    sum_v = 0.0; sum_w = 0.0
    for (k,i,j) in B:
        xi = var_value_safe(x.get((k,i)))
        if xi > 0.5:
            sum_v += var_value_safe(v.get((k,i,j)))
            sum_w += var_value_safe(w.get((k,i,j)))
    left_behind_rate = (sum_v / sum_w) if sum_w > eps else 0.0
    dcr_by_ki = {}
    dcr_details = {}
    peak_dcr = 0.0
    for (k,i) in list(w_b.keys()):
        tau_val = var_value_safe(tau.get(k) if isinstance(tau, dict) else None, 1.0)
        arrivals = var_value_safe(w_b.get((k,i)), 0.0)
        idx = (k-1, k)
        if idx in h:
            headway = var_value_safe(h.get(idx), float(h_min))
        else:
            headway = float(h_min) if tau_val > 0.5 else 0.0
        if headway > eps:
            cap_per_time = cap / headway
            dcr_val = arrivals / cap_per_time if cap_per_time > eps else float('inf')
        else:
            dcr_val = 0.0 if arrivals <= eps else float('inf')
        dcr_by_ki[(k,i)] = float(dcr_val)
        dcr_details[(k,i)] = (float(dcr_val), float(arrivals), float(headway), float(cap))
        if dcr_val != float('inf') and dcr_val > peak_dcr:
            peak_dcr = float(dcr_val)
    return {
        'max_plf': float(max_plf),
        'overloaded': overloaded,
        'cong_seconds': float(cong_seconds),
        'cong_minutes': float(cong_minutes),
        'left_behind_rate': float(left_behind_rate),
        'dcr_by_ki': dcr_by_ki,
        'dcr_details': dcr_details,
        'peak_dcr': float(peak_dcr)
    }
def update_services_from_solution(model, services_src):
    if isinstance(services_src, str):
        try:
            services_src_df = pd.read_csv(services_src)
        except Exception:
            return
    else:
        services_src_df = services_src
    global d, x, H, T_WINDOW
    for idx, row in services_src_df.iterrows():
        try:
            k = row['k']
            station = row['station_i']
            if (k, station) in d.keys():
                var = d[k, station]
                if var is not None and getattr(var, 'X', None) is not None:
                    depart_min = (float(var.X) - float(H)) / 60.0
                    depart_min = max(0.0, depart_min)
                    if depart_min > T_WINDOW / 60.0:
                        depart_min = T_WINDOW / 60.0
                    services_src_df.at[idx, 'depart_time_min'] = depart_min
                    if (k, station) in x.keys():
                        try:
                            services_src_df.at[idx, 'x_stop'] = int(round(x[(k, station)].X))
                        except Exception:
                            pass
        except Exception:
            continue
    if not isinstance(services_src, str):
        for col in services_src_df.columns:
            services_src[col] = services_src_df[col]
    else:
        services_src_df.to_csv(services_src, index=False)
#Main execution: assemble PWLs, build model, solve, run preprocessor, metrics, plotting
def main():
    try:
        df_rt = pd.read_excel('TRB_16_stations.xlsx', 'Running_time')
        r = dict((int(a), int(b), float(c)) if False else None)  
    except Exception:
        pass
    try:
        df1 = pd.read_excel('TRB_16_stations.xlsx', 'Running_time')
        r = dict((tuple((a, b)), c) for a,b,c in df1.values)
    except Exception:
        r = {}
    try:
        df1 = pd.read_excel('TRB_16_stations.xlsx', 'Dwelling_time')
        my_tuple = [tuple(x) for x in df1.values]
        e = dict((a, c) for a, c in my_tuple)
    except Exception:
        e = {}
    try:
        df1 = pd.read_excel('TRB_16_stations.xlsx', 'Demand_30(7.30-8.00am)mint')
        p = dict((tuple((a, b)), c) for a,b,c in df1.values)
    except Exception:
        p = {}
    phi = build_phi_from_excel()
    P_backlog = load_backlog_from_excel()
    if P_backlog is None:
        PRE_BACK_MIN = 2.0
        factor = PRE_BACK_MIN / 30.0
        P_backlog = {(int(i), int(j)): factor * float(val) for (i, j), val in p.items()}
    _pwl_up = build_pwl_breaks_for_dir(p, phi, list(range(1, S_U+1)))
    _pwl_down = build_pwl_breaks_for_dir(p, phi, list(range(9, S_D+1)))
    arrivals_by_station_od = {}
    for (i,j), pwl in _pwl_up.items():
        arrivals_by_station_od.setdefault(i, {})[(i,j)] = pwl_to_minute_series(pwl)
    for (i,j), pwl in _pwl_down.items():
        arrivals_by_station_od.setdefault(i, {})[(i,j)] = pwl_to_minute_series(pwl)
    mdl, var_dicts = build_trb_model(p, r, e)
    global d, x
    d = var_dicts.get('d', {})
    x = var_dicts.get('x', {})
    try:
        mdl.optimize()
    except Exception as _e:
        print("Initial solve failed:", _e)
    services_rows = []
    K_u = var_dicts.get('K_u', [])
    K_d = var_dicts.get('K_d', [])
    S_u = var_dicts.get('S_u', [])
    S_d = var_dicts.get('S_d', [])
    for k in K_u:
        for i in S_u:
            xstop = 1
            try:
                if (k, i) in x.keys() and hasattr(x[(k,i)], 'X'):
                    xstop = int(round(x[(k, i)].X))
                elif (k, i) in x.keys():
                    xstop = int(x[(k, i)].LB)
            except Exception:
                xstop = 1
            depart_min = 0.0
            try:
                if (k, i) in d.keys() and hasattr(d[(k, i)], 'X'):
                    depart_min = max(0.0, (float(d[(k, i)].X) - float(H)) / 60.0)
            except Exception:
                depart_min = 0.0
            services_rows.append({'station_i': i, 'k': k, 'depart_time_min': depart_min, 'x_stop': xstop, 'capacity': C, 'onboard_arrival': 0.0})
    for k in K_d:
        for i in S_d:
            xstop = 1
            try:
                if (k, i) in x.keys() and hasattr(x[(k, i)], 'X'):
                    xstop = int(round(x[(k, i)].X))
                elif (k, i) in x.keys():
                    xstop = int(x[(k, i)].LB)
            except Exception:
                xstop = 1
            depart_min = 0.0
            try:
                if (k, i) in d.keys() and hasattr(d[(k, i)], 'X'):
                    depart_min = max(0.0, (float(d[(k, i)].X) - float(H)) / 60.0)
            except Exception:
                depart_min = 0.0
            services_rows.append({'station_i': i, 'k': k, 'depart_time_min': depart_min, 'x_stop': xstop, 'capacity': C, 'onboard_arrival': 0.0})

    services_df = pd.DataFrame(services_rows)
    services_df_before = services_df.copy()
    initial_backlogs_per_station = {}
    for st in sorted(arrivals_by_station_od.keys()):
        sub = {od: val for od, val in P_backlog.items() if od[0] == st}
        initial_backlogs_per_station[st] = sub
    var_maps = {
        'w_var': var_dicts.get('w', {}),
        'nb_var': var_dicts.get('n_b1', {}),
        'nb_total_var': var_dicts.get('n_b', {})
    }
    metrics_before = collect_performance_metrics_from_model(mdl, var_dicts)
    def write_model_vars_dump(model, var_dicts, basename, stage):
        out_path = f"{basename}_{stage}_iterative.txt"
        try:
            K_u = var_dicts.get('K_u', [])
            K_d = var_dicts.get('K_d', [])
        except Exception:
            K_u = []
            K_d = []
        vars_by_prefix = {}
        for v in model.getVars():
            name = v.VarName
            if name is None:
                continue
            if '[' in name:
                prefix = name.split('[')[0]
            else:
                prefix = name
            vars_by_prefix.setdefault(prefix, []).append(v)
        def vv(var):
            try:
                val = getattr(var, 'X', None)
                if val is None:
                    val = getattr(var, 'x', None)
                return float(val) if val is not None else None
            except Exception:
                return None
        with open(out_path, 'w', encoding='utf-8') as f:
            # If available, print iteration-wise precomputed params for 'after' stage BEFORE runtime
            if stage == 'after' and PREPROC_HISTORY:
                f.write('Iteration-wise precomputed queue parameters (used in next iteration):\n')
                for entry in PREPROC_HISTORY:
                    it = entry.get('iter')
                    f.write(f'Iteration {it}:\n')
                    pw = entry.get('pre_w', {})
                    pnb = entry.get('pre_nb', {})
                    pv = entry.get('pre_v', {})
                    pnt = entry.get('pre_nb_total', {})
                    f.write('  pre_nb (nb_boarded):\n')
                    for k in sorted(pnb.keys()):
                        if pnb[k] > 0:
                            f.write(f'    {k}: {pnb[k]}\n')
                    f.write('  pre_w (queue before):\n')
                    for k in sorted(pw.keys()):
                        if pw[k] > 0:
                            f.write(f'    {k}: {pw[k]}\n')
                    f.write('  pre_v (leftover v):\n')
                    for k in sorted(pv.keys()):
                        if pv[k] > 0:
                            f.write(f'    {k}: {pv[k]}\n')
                    f.write('  pre_nb_total (nb_total):\n')
                    for k in sorted(pnt.keys()):
                        if pnt[k] > 0:
                            f.write(f'    {k}: {pnt[k]}\n')
                    f.write('\n')
            rt = getattr(model, 'Runtime', None)
            if rt is not None:
                f.write(f"The run time is {rt}\n\n")
            else:
                f.write('The run time is N/A\n\n')
            try:
                obj = getattr(model, 'ObjVal', None)
                f.write('\nObjective Function Value:\n')
                f.write(f'ObjVal: {obj}\n')
            except Exception:
                f.write('\nObjective Function Value:\nObjVal: N/A\n')
            if 'tau' in vars_by_prefix:
                f.write('\nValues of binary variable tau_u:\n')
                tau_vars = {v.VarName: v for v in vars_by_prefix['tau']}
                for k in sorted(K_u):
                    name = f"tau[{k}]"
                    val = vv(tau_vars.get(name)) if name in tau_vars else None
                    if val is not None and round(val) == 1:
                        f.write(f"tau[{k}]: {val}\n")
                f.write('\nValues of binary variable tau_d:\n')
                for k in sorted(K_d):
                    name = f"tau[{k}]"
                    val = vv(tau_vars.get(name)) if name in tau_vars else None
                    if val is not None and round(val) == 1:
                        f.write(f"tau[{k}]: {val}\n")
            for rs in ['RS1', 'RS2', 'RS3', 'RS4']:
                f.write(f"\nValue of integer variable {rs}:\n")
                vlist = vars_by_prefix.get(rs, [])
                if vlist:
                    val = vv(vlist[0])
                    if val is not None and val > 0:
                        f.write(f"{rs}: {val}\n")
                else:
                    f.write(f"{rs}: \n")
            pref_order = ['h', 'z', 'y', 'x', 'a1', 'd', 'alpha_u1', 'alpha_u3', 'alpha_d2', 'alpha_d4',
                          'beta_u2', 'beta_u4', 'beta_d1', 'beta_d3', 'w', 'w_b', 'w_b1', 'n_b', 'n_b1',
                          'n_a', 'n1', 'v', 'sai']
            for prefix in pref_order:
                if prefix in vars_by_prefix:
                    f.write(f"\nValues of {'continuous' if prefix in ['h','a1','d','w','w_b','w_b1','n_b','n_b1','n_a','n1','v'] else 'binary'} variable {prefix}:\n")
                    def keyfn(var):
                        try:
                            name = var.VarName
                            inside = name[name.find('[')+1:name.rfind(']')]
                            parts = [p.strip() for p in inside.split(',')]
                            return tuple(float(p) if '.' in p else int(p) for p in parts)
                        except Exception:
                            return (name,)
                    for v in sorted(vars_by_prefix[prefix], key=keyfn):
                        val = vv(v)
                        is_binary = prefix in ['tau', 'z', 'y', 'x', 'alpha_u1', 'alpha_u3', 'alpha_d2', 'alpha_d4','beta_u2', 'beta_u4', 'beta_d1', 'beta_d3', 'sai']
                        if is_binary:
                            if val is not None and round(val) == 1:
                                f.write(f"{v.VarName}: {val}\n")
                        else:
                            if val is not None and val > 0:
                                f.write(f"{v.VarName}: {val}\n")
            other = [p for p in vars_by_prefix.keys() if p not in pref_order and p not in ('tau', 'RS1','RS2','RS3','RS4')]
            for prefix in sorted(other):
                f.write(f"\nValues of variable {prefix}:\n")
                for v in sorted(vars_by_prefix[prefix], key=lambda v: v.VarName):
                    val = vv(v)
                    if val is not None and val > 0:
                        f.write(f"{v.VarName}: {val}\n")
        print(f"[Saved] variable dump -> {out_path}")
    #Helper: compute adaptive pre_penalty with logging and fallback ---
    def compute_and_apply_pre_penalty(model, preproc_history, alpha=0.1, min_penalty=0.1, beta=1.0):
        try:
            pre = compute_pre_penalty_from_model(model, preproc_history, alpha=alpha, min_penalty=min_penalty, beta=beta)
        except Exception as e:
            print('[Preproc] Failed to compute adaptive pre_penalty, using fallback 1.0:', e)
            pre = 1.0
        try:
            print(f"[Preproc] adaptive pre_penalty = {pre:.6f}")
        except Exception:
            print('[Preproc] adaptive pre_penalty computed.')
        return float(pre)
    print('\n[Preprocessor] Running iterative coupling (max_iters=3)...')
    pre_penalty = compute_and_apply_pre_penalty(mdl, PREPROC_HISTORY, alpha=0.1, min_penalty=0.1, beta=1.0)
    run_iterative_with_preprocessor(mdl, sorted(arrivals_by_station_od.keys()), arrivals_by_station_od, services_df,var_maps, max_iters=3, tol=1e-3, mu_board=0.5, dwell_base_sec=20.0,dwell_per_pax_sec=0.25, initial_backlogs=initial_backlogs_per_station, mode='soft', pre_penalty=pre_penalty)
    print('[Preprocessor] Iterative coupling complete.\n')
    try:
        metrics_after = collect_performance_metrics_from_model(mdl, var_dicts)
        with open('Qmetrics_comparison_QTRB_30_mint_Model_1_Off_peak_Morning_14_trains.txt', 'w') as cf:
            cf.write('Comparison of system performance metrics\n\n')
            cf.write('METRIC | BEFORE_ITERATIVE | AFTER_ITERATIVE\n')
            cf.write('----------------------------------------------------------\n')
            cf.write('Peak Load Factor (max PLF) | {:.6f} | {:.6f}\n'.format(metrics_before.get('max_plf',0.0), metrics_after.get('max_plf',0.0)))
            cf.write('Congestion minutes         | {:.3f} | {:.3f}\n'.format(metrics_before.get('cong_minutes',0.0), metrics_after.get('cong_minutes',0.0)))
            cf.write('Left-Behind Rate           | {:.6f} | {:.6f}\n'.format(metrics_before.get('left_behind_rate',0.0), metrics_after.get('left_behind_rate',0.0)))
            cf.write('Number overloaded (PLF>0.85)| {} | {}\n'.format(len(metrics_before.get('overloaded',[])), len(metrics_after.get('overloaded',[]))))
    except Exception as _merr:
        print('Failed to write metrics after iterative preprocessor:', _merr)
if __name__ == "__main__":
    main()