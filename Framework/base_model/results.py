# results.py

from config import SOLUTION_FILE


def log_message(file, message):

    print(message)

    file.write(message + "\n")


def print_results(
    mdl,
    vars,
    sets
):

    with open(SOLUTION_FILE, "w") as file:

        log_message(file, "=================================")
        log_message(file, "OPTIMIZATION RESULTS")
        log_message(file, "=================================")

        runtime = mdl.Runtime

        log_message(file, "")
        log_message(file, f"Run Time: {runtime:.4f} seconds")

        log_message(file, "")
        log_message(file, f"Objective Function Value = {mdl.ObjVal}")

        # =========================
        # VARIABLES
        # =========================

        tau = vars['tau']
        z = vars['z']
        y = vars['y']
        x = vars['x']

        a1 = vars['a1']
        d = vars['d']
        h = vars['h']

        alpha_up1 = vars['alpha_up1']
        alpha_dn2 = vars['alpha_dn2']
        alpha_up3 = vars['alpha_up3']
        alpha_dn4 = vars['alpha_dn4']

        beta_up2 = vars['beta_up2']
        beta_dn1 = vars['beta_dn1']
        beta_up4 = vars['beta_up4']
        beta_dn3 = vars['beta_dn3']

        RS1 = vars['RS1']
        RS2 = vars['RS2']
        RS3 = vars['RS3']
        RS4 = vars['RS4']

        w = vars['w']
        w_b = vars['w_b']
        w_b1 = vars['w_b1']

        n_b = vars['n_b']
        n_b1 = vars['n_b1']
        n_a = vars['n_a']

        n1 = vars['n1']
        v = vars['v']
        sai = vars['sai']

        # =========================
        # SETS
        # =========================

        K_u = sets['K_u']
        K_d = sets['K_d']

        B = sets['B']
        A = sets['A']
        C1 = sets['C1']
        D = sets['D']

        h_min = sets.get('h_min', 90)
        cap = sets.get('C', 250)

        # =========================
        # tau
        # =========================

        log_message(file, "\nValues of tau:")

        for k in K_u:

            if tau[k].X > 0:

                log_message(file, f"tau[{k}] = {tau[k].X}")

        for k in K_d:

            if tau[k].X > 0:

                log_message(file, f"tau[{k}] = {tau[k].X}")

        # =========================
        # RS Variables
        # =========================

        log_message(file, "\nRolling Stock Variables:")

        if RS1.X > 0:
            log_message(file, f"RS1 = {RS1.X}")

        if RS2.X > 0:
            log_message(file, f"RS2 = {RS2.X}")

        if RS3.X > 0:
            log_message(file, f"RS3 = {RS3.X}")

        if RS4.X > 0:
            log_message(file, f"RS4 = {RS4.X}")

        # =========================
        # h
        # =========================

        log_message(file, "\nValues of h:")

        for d_val in D:

            if h[d_val].X > 0:

                log_message(file, f"h[{d_val}] = {h[d_val].X}")

        # =========================
        # z
        # =========================

        log_message(file, "\nValues of z:")

        for b in B:

            if z[b].X > 0:

                log_message(file, f"z[{b}] = {z[b].X}")

        # =========================
        # y
        # =========================

        log_message(file, "\nValues of y:")

        for a in A:

            if y[a].X > 0:

                log_message(file, f"y[{a}] = {y[a].X}")

        # =========================
        # x
        # =========================

        log_message(file, "\nValues of x:")

        for c in C1:

            if x[c].X > 0:

                log_message(file, f"x[{c}] = {x[c].X}")

        # =========================
        # a1
        # =========================

        log_message(file, "\nValues of a1:")

        for c in C1:

            if a1[c].X > 0:

                log_message(file, f"a1[{c}] = {a1[c].X}")

        # =========================
        # d
        # =========================

        log_message(file, "\nValues of d:")

        for c in C1:

            if d[c].X > 0:

                log_message(file, f"d[{c}] = {d[c].X}")

        # =========================
        # Alpha Variables
        # =========================

        log_message(file, "\nAlpha Variables:")

        for k in K_u:

            if alpha_up1[k].X > 0:

                log_message(file, f"alpha_up1[{k}] = {alpha_up1[k].X}")

        for k in K_d:

            if alpha_dn2[k].X > 0:

                log_message(file, f"alpha_dn2[{k}] = {alpha_dn2[k].X}")

        for k in K_u:

            if alpha_up3[k].X > 0:

                log_message(file, f"alpha_up3[{k}] = {alpha_up3[k].X}")

        for k in K_d:

            if alpha_dn4[k].X > 0:

                log_message(file, f"alpha_dn4[{k}] = {alpha_dn4[k].X}")

        # =========================
        # Beta Variables
        # =========================

        log_message(file, "\nBeta Variables:")

        for k in K_u:

            if beta_up2[k].X > 0:

                log_message(file, f"beta_up2[{k}] = {beta_up2[k].X}")

        for k in K_d:

            if beta_dn1[k].X > 0:

                log_message(file, f"beta_dn1[{k}] = {beta_dn1[k].X}")

        for k in K_u:

            if beta_up4[k].X > 0:

                log_message(file, f"beta_up4[{k}] = {beta_up4[k].X}")

        for k in K_d:

            if beta_dn3[k].X > 0:

                log_message(file, f"beta_dn3[{k}] = {beta_dn3[k].X}")

        # =========================
        # Passenger Variables
        # =========================

        log_message(file, "\nValues of w:")

        for b in B:

            if w[b].X > 0:

                log_message(file, f"w[{b}] = {w[b].X}")

        log_message(file, "\nValues of w_b:")

        for c in C1:

            if w_b[c].X > 0:

                log_message(file, f"w_b[{c}] = {w_b[c].X}")

        log_message(file, "\nValues of w_b1:")

        for b in B:

            if w_b1[b].X > 0:

                log_message(file, f"w_b1[{b}] = {w_b1[b].X}")

        log_message(file, "\nValues of n_b:")

        for c in C1:

            if n_b[c].X > 0:

                log_message(file, f"n_b[{c}] = {n_b[c].X}")

        log_message(file, "\nValues of n_b1:")

        for b in B:

            if n_b1[b].X > 0:

                log_message(file, f"n_b1[{b}] = {n_b1[b].X}")

        log_message(file, "\nValues of n_a:")

        for c in C1:

            if n_a[c].X > 0:

                log_message(file, f"n_a[{c}] = {n_a[c].X}")

        log_message(file, "\nValues of n1:")

        for c in C1:

            if n1[c].X > 0:

                log_message(file, f"n1[{c}] = {n1[c].X}")

        log_message(file, "\nValues of v:")

        for b in B:

            if v[b].X > 0:

                log_message(file, f"v[{b}] = {v[b].X}")

        log_message(file, "\nValues of sai:")

        for c in C1:

            if sai[c].X > 0:

                log_message(file, f"sai[{c}] = {sai[c].X}")

        # =========================
        # PERFORMANCE METRICS
        # =========================

        eps = 1e-6

        plf_by_ki = {}

        for (k, i) in n1.keys():

            try:

                if tau[k].X > 0.5:

                    val = n1[(k, i)].X

                    plf = val / float(cap)

                    plf_by_ki[(k, i)] = plf

            except Exception:
                continue

        max_plf = max(plf_by_ki.values()) if plf_by_ki else 0.0

        overloaded = [

            (k, i, p)

            for (k, i), p in plf_by_ki.items()

            if p > 0.85
        ]

        cong_seconds = 0.0

        services_with_cong = set([k for (k, i, p) in overloaded])

        for k in services_with_cong:

            idx = (k - 1, k)

            if idx in h.keys():

                cong_seconds += h[idx].X

            else:

                cong_seconds += 120.0

        cong_minutes = cong_seconds / 60.0

        sum_v = 0.0
        sum_w = 0.0

        for (k, i, j) in B:

            try:

                if (

                    x[(k, i)].X > 0.5
                    and
                    x[(k, j)].X > 0.5

                ):

                    sum_v += float(v[(k, i, j)].X)

                    sum_w += float(w[(k, i, j)].X)

            except Exception:
                continue

        left_behind_rate = (

            sum_v / sum_w

            if sum_w > eps else 0.0
        )

        dcr_by_ki = {}

        for (k, i) in w_b.keys():

            try:

                if tau[k].X > 0.5:

                    arrivals = float(w_b[(k, i)].X)

                    idx = (k - 1, k)

                    if idx in h.keys():

                        hw = max(h[idx].X, eps)

                    else:

                        hw = h_min

                    cap_per_time = float(cap) / hw

                    dcr_val = arrivals / cap_per_time

                    dcr_by_ki[(k, i)] = dcr_val

            except Exception:
                continue

        peak_dcr = max(dcr_by_ki.values()) if dcr_by_ki else 0.0

        log_message(file, "\n=================================")
        log_message(file, "PERFORMANCE METRICS")
        log_message(file, "=================================")

        log_message(
            file,
            f"Peak Load Factor (PLF) = {max_plf:.4f}"
        )

        if overloaded:

            log_message(
                file,
                "Overloaded Services (PLF > 0.85):"
            )

            for (k, i, p) in overloaded:

                log_message(
                    file,
                    f"({k},{i}) -> {p:.4f}"
                )

        else:

            log_message(
                file,
                "No overloaded services."
            )

        log_message(
            file,
            f"Congestion Duration = {cong_seconds:.2f} sec ({cong_minutes:.2f} min)"
        )

        log_message(
            file,
            f"Left-Behind Passenger Rate = {left_behind_rate:.4f}"
        )

        log_message(
            file,
            f"Peak Demand-to-Capacity Ratio = {peak_dcr:.4f}"
        )