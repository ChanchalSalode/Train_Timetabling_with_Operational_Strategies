# precompute_params.py

from gurobipy import quicksum

from queue_simulator import (
    simulate_queues_station
)


# =========================================
# COMPUTE PRECOMPUTED PARAMETERS
# =========================================

def compute_precomputed_params(

    all_stations,

    arrivals_by_station_od,

    services_df_or_csvpath,

    mu_board=0.5,

    dwell_base_sec=20.0,

    dwell_per_pax_sec=0.25,

    initial_backlogs=None

):

    import pandas as pd

    # =====================================
    # Storage Dictionaries
    # =====================================

    pre_w = {}

    pre_nb = {}

    pre_v = {}

    pre_nb_total = {}

    # =====================================
    # Load Services
    # =====================================

    if isinstance(

        services_df_or_csvpath,
        str

    ):

        services_all = pd.read_csv(
            services_df_or_csvpath
        )

    else:

        services_all = (
            services_df_or_csvpath.copy()
        )

    # =====================================
    # Process Each Station
    # =====================================

    for station in all_stations:

        arrivals_by_od = (

            arrivals_by_station_od.get(
                station,
                {}
            )

        )

        # ================================
        # Services for Station
        # ================================

        df_s = services_all[

            services_all['station_i']
            == station

        ]

        services_list = []

        for _, row in df_s.iterrows():

            services_list.append({

                'k':
                row['k'],

                'depart_time_min':
                float(
                    row['depart_time_min']
                ),

                'x_stop':
                int(
                    row.get(
                        'x_stop',
                        1
                    )
                ),

                'capacity':
                float(
                    row.get(
                        'capacity',
                        1e9
                    )
                ),

                'onboard_arrival':
                float(
                    row.get(
                        'onboard_arrival',
                        0.0
                    )
                )

            })

        # ================================
        # Initial Backlogs
        # ================================

        init_back = (

            initial_backlogs.get(station)

            if initial_backlogs

            else None

        )

        # ================================
        # Run Queue Simulator
        # ================================

        df_out = simulate_queues_station(

            arrivals_by_od,

            services_list,

            mu_board=mu_board,

            dwell_base_sec=dwell_base_sec,

            dwell_per_pax_sec=dwell_per_pax_sec,

            initial_backlog=init_back

        )

        # ================================
        # Extract Parameters
        # ================================

        for _, r in df_out.iterrows():

            k = r['k']

            i = r['i']

            j = r['j']

            pre_w[(k, i, j)] = float(
                r['w_before']
            )

            pre_nb[(k, i, j)] = float(
                r['nb_boarded']
            )

            pre_v[(k, i, j)] = float(
                r['v_left']
            )

            pre_nb_total[(k, i)] = float(
                r['nb_total']
            )

    return (

        pre_w,
        pre_nb,
        pre_v,
        pre_nb_total

    )


# =========================================
# COMPUTE ADAPTIVE PRE-PENALTY
# =========================================

def compute_pre_penalty_from_model(

    model,

    preproc_history,

    alpha=0.1,

    min_penalty=1e-3,

    beta=1.0

):

    # =====================================
    # Base Objective
    # =====================================

    try:

        Obj_base = abs(model.ObjVal)

        if Obj_base is None:

            Obj_base = None

    except Exception:

        Obj_base = None

    # =====================================
    # Fallback Objective
    # =====================================

    if Obj_base is None:

        Obj_base = 0.0

        if preproc_history:

            last = preproc_history[-1]

            pnb = (

                last.get('pre_nb', {})

                if isinstance(last, dict)

                else {}

            )

            pnt = (

                last.get(
                    'pre_nb_total',
                    {}
                )

                if isinstance(last, dict)

                else {}

            )

            Obj_base = float(

                max(

                    1.0,

                    sum(pnb.values())
                    +
                    sum(pnt.values())

                )

            )

        else:

            Obj_base = 1.0

    # =====================================
    # Slack Estimate
    # =====================================

    if preproc_history:

        last = preproc_history[-1]

        N_slack_est = max(

            1,

            len(

                last.get(
                    'pre_nb_total',
                    {}
                )

                if isinstance(last, dict)

                else {}

            )

        )

    else:

        N_slack_est = 1

    # =====================================
    # Penalty Calculation
    # =====================================

    per_slack = (

        alpha
        *
        float(Obj_base)
        /
        float(N_slack_est)

    )

    max_total = (

        beta
        *
        float(
            max(
                1.0,
                Obj_base
            )
        )

    )

    per_slack = min(

        per_slack,

        max_total
        /
        float(N_slack_est)

    )

    per_slack = max(

        per_slack,

        float(min_penalty)

    )

    return float(per_slack)


# =========================================
# ADD PRECOMPUTED CONSTRAINTS
# =========================================

def add_precomputed_constraints_to_model(

    model,

    pre_w,

    pre_nb,

    pre_nb_total,

    var_maps,

    mode='ub',

    pre_penalty=10000.0

):

    # =====================================
    # Remove Old PRE Constraints
    # =====================================

    try:

        to_remove = [

            c

            for c in model.getConstrs()

            if (

                c.ConstrName

                and

                c.ConstrName.startswith(
                    'PRE_'
                )

            )

        ]

        if to_remove:

            model.remove(to_remove)

            model.update()

    except Exception:

        pass

    # =====================================
    # Remove Old Slack Variables
    # =====================================

    try:

        prev_slacks = [

            v

            for v in model.getVars()

            if (

                v.VarName

                and

                v.VarName.startswith(
                    'PRE_s_'
                )

            )

        ]

        if prev_slacks:

            model.remove(prev_slacks)

            model.update()

    except Exception:

        pass

    # =====================================
    # Variable Dictionaries
    # =====================================

    created_slacks = []

    nb_dict = var_maps.get(
        'nb_var',
        {}
    )

    nbtot_dict = var_maps.get(
        'nb_total_var',
        {}
    )

    wdict = var_maps.get(
        'w_var',
        {}
    )

    # =====================================
    # PRE nb Constraints
    # =====================================

    for (k, i, j), val in pre_nb.items():

        if (k, i, j) in nb_dict:

            var = nb_dict[(k, i, j)]

            # ================================
            # Upper Bound
            # ================================

            if mode == 'ub':

                cname = (
                    f"PRE_nb_ub_{k}_{i}_{j}"
                )

                model.addConstr(

                    var <= val + 1e-6,

                    name=cname

                )

            # ================================
            # Fixed
            # ================================

            elif mode == 'fix':

                cname = (
                    f"PRE_nb_eq_{k}_{i}_{j}"
                )

                model.addConstr(

                    var == val,

                    name=cname

                )

            # ================================
            # Soft Constraint
            # ================================

            elif mode == 'soft':

                s = model.addVar(

                    lb=0.0,

                    name=(
                        f"PRE_s_nb_"
                        f"{k}_{i}_{j}"
                    )

                )

                model.addConstr(

                    var <= val + s,

                    name=(
                        f"PRE_nb_soft_"
                        f"{k}_{i}_{j}"
                    )

                )

                created_slacks.append(s)

    # =====================================
    # PRE nb_total Constraints
    # =====================================

    for (k, i), val in pre_nb_total.items():

        if (k, i) in nbtot_dict:

            var = nbtot_dict[(k, i)]

            # ================================
            # Upper Bound
            # ================================

            if mode == 'ub':

                cname = (
                    f"PRE_nbtotal_ub_"
                    f"{k}_{i}"
                )

                model.addConstr(

                    var <= val + 1e-6,

                    name=cname

                )

            # ================================
            # Fixed
            # ================================

            elif mode == 'fix':

                cname = (
                    f"PRE_nbtotal_eq_"
                    f"{k}_{i}"
                )

                model.addConstr(

                    var == val,

                    name=cname

                )

            # ================================
            # Soft Constraint
            # ================================

            elif mode == 'soft':

                s = model.addVar(

                    lb=0.0,

                    name=(
                        f"PRE_s_nbtotal_"
                        f"{k}_{i}"
                    )

                )

                model.addConstr(

                    var <= val + s,

                    name=(
                        f"PRE_nbtotal_soft_"
                        f"{k}_{i}"
                    )

                )

                created_slacks.append(s)

    # =====================================
    # Update Model
    # =====================================

    model.update()

    # =====================================
    # Add Penalty
    # =====================================

    if created_slacks:

        try:

            cur_obj = model.getObjective()

        except Exception:

            cur_obj = None

        penalty_expr = (

            pre_penalty
            *
            quicksum(created_slacks)

        )

        if cur_obj is not None:

            model.setObjective(

                cur_obj - penalty_expr

            )

        else:

            model.setObjective(

                -penalty_expr

            )

        model.update()

    return created_slacks
