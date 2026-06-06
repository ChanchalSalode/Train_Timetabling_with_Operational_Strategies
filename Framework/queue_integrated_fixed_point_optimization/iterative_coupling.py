# iterative_coupling.py

import traceback

from gurobipy import GRB

from precompute_params import (
    compute_precomputed_params,
    add_precomputed_constraints_to_model,
    compute_pre_penalty_from_model
)

from helpers import (
    average_absolute_difference
)


# =========================================
# GLOBAL HISTORY STORAGE
# =========================================

PREPROC_HISTORY = []


# =========================================
# COMPUTE ADAPTIVE PRE-PENALTY
# =========================================

def compute_and_apply_pre_penalty(
    model,
    preproc_history,
    alpha=0.1,
    min_penalty=0.1,
    beta=1.0
):

    try:

        pre = compute_pre_penalty_from_model(
            model,
            preproc_history,
            alpha=alpha,
            min_penalty=min_penalty,
            beta=beta
        )

    except Exception as e:

        print(
            '[Preproc] Failed to compute '
            'adaptive pre_penalty, '
            'using fallback 1.0:',
            e
        )

        pre = 1.0

    try:

        print(
            f"[Preproc] adaptive "
            f"pre_penalty = {pre:.6f}"
        )

    except Exception:

        print(
            '[Preproc] adaptive '
            'pre_penalty computed.'
        )

    return float(pre)


# =========================================
# ITERATIVE PREPROCESSOR + SOLVE
# =========================================

def run_iterative_with_preprocessor(

    model,

    all_stations,

    arrivals_by_station_od,

    services_csv_or_df,

    var_maps,

    max_iters=3,

    tol=1e-3,

    mu_board=0.5,

    dwell_base_sec=20.0,

    dwell_per_pax_sec=0.25,

    initial_backlogs=None,

    mode='soft',

    pre_penalty=10000.0,

    update_services_function=None

):

    prev_pre_nb = None

    services_src = services_csv_or_df

    # =====================================
    # Iterative Loop
    # =====================================

    for itr in range(

        1,
        max_iters + 1

    ):

        print(
            f"[Preproc Iter {itr}] "
            f"Computing precomputed params..."
        )

        # =================================
        # Compute Queue Approximation
        # =================================

        (
            pre_w,
            pre_nb,
            pre_v,
            pre_nb_total

        ) = compute_precomputed_params(

            all_stations,

            arrivals_by_station_od,

            services_src,

            mu_board=mu_board,

            dwell_base_sec=dwell_base_sec,

            dwell_per_pax_sec=dwell_per_pax_sec,

            initial_backlogs=initial_backlogs

        )

        # =================================
        # Store Iteration History
        # =================================

        try:

            PREPROC_HISTORY.append({

                'iter': itr,

                'pre_w': pre_w.copy(),

                'pre_nb': pre_nb.copy(),

                'pre_v': pre_v.copy(),

                'pre_nb_total': pre_nb_total.copy()

            })

        except Exception:

            PREPROC_HISTORY.append({

                'iter': itr

            })

        # =================================
        # Add PRE Constraints
        # =================================

        print(
            f"[Preproc Iter {itr}] "
            f"Adding constraints to model..."
        )

        add_precomputed_constraints_to_model(

            model,

            pre_w,

            pre_nb,

            pre_nb_total,

            var_maps,

            mode=mode

        )

        # =================================
        # Optimize
        # =================================

        print(
            f"[Preproc Iter {itr}] "
            f"Solving MIP..."
        )

        model.optimize()

        # =================================
        # Diagnostics
        # =================================

        try:

            n_pre = sum(
                1
                for c in model.getConstrs()
                if (
                    c.ConstrName
                    and
                    c.ConstrName.startswith('PRE_')
                )
            )

            n_slack_vars = sum(
                1
                for v in model.getVars()
                if (
                    v.VarName
                    and
                    v.VarName.startswith('PRE_s_')
                )
            )

            penalty_contrib = 0.0

            for v in model.getVars():

                if (
                    v.VarName
                    and
                    v.VarName.startswith('PRE_s_')
                ):

                    xv = getattr(v, 'X', None)

                    if xv is None:

                        xv = getattr(v, 'x', 0.0)

                    try:

                        penalty_contrib += float(xv)

                    except Exception:

                        pass

            penalty_contrib *= float(pre_penalty)

            print(
                f"[Preproc Diagnostics] "
                f"PRE constraints={n_pre}, "
                f"PRE slacks={n_slack_vars}, "
                f"pre_penalty={pre_penalty}, "
                f"penalty_contrib={penalty_contrib}"
            )

        except Exception as diag_error:

            print(
                "Failed to compute preproc diagnostics:",
                diag_error
            )

        status = model.Status
        # =================================
        # Infeasible Case
        # =================================

        if status == GRB.INFEASIBLE:

            print(

                f"[Preproc Iter {itr}] "

                f"Model became INFEASIBLE "

                f"after PRE constraints."

            )

            print(
                "Computing IIS..."
            )

            try:

                iis_name = (
                    f"mdl_preproc_iter{itr}.ilp"
                )

                model.computeIIS()

                model.write(iis_name)

                print(
                    f"[Preproc Iter {itr}] "
                    f"IIS written to "
                    f"{iis_name}."
                )

            except Exception as e:

                print(
                    f"[Preproc Iter {itr}] "
                    f"Failed IIS: {e}"
                )

            # =============================
            # Remove PRE Constraints
            # =============================

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

            print(

                f"[Preproc Iter {itr}] "

                f"PRE constraints removed."

            )

            break

        # =================================
        # Bad Solver Status
        # =================================

        if status not in (

            GRB.OPTIMAL,

            GRB.TIME_LIMIT,

            GRB.SUBOPTIMAL,

            GRB.INTERRUPTED

        ):

            print(

                f"[Preproc Iter {itr}] "

                f"Solver returned status "

                f"{status}. "

                f"Stopping iterative preprocessor."

            )

            break

        # =================================
        # Update Services
        # =================================

        if update_services_function is not None:

            try:

                update_services_function(
                    model,
                    services_src
                )

            except Exception:

                print(

                    "NOTE: "

                    "update_services_from_solution "

                    "failed."

                )

                traceback.print_exc()

        # =================================
        # Convergence Check
        # =================================

        if prev_pre_nb is not None:

            diff = average_absolute_difference(

                pre_nb,
                prev_pre_nb

            )

            print(

                f"[Preproc Iter {itr}] "

                f"avg change in pre_nb "

                f"= {diff:.6f}"

            )

            if diff < tol:

                print(

                    f"[Preproc Iter {itr}] "

                    f"Converged "

                    f"(tol={tol}). "
                    f"Stopping iter."

                )

                break

        # =================================
        # Store Current Iteration
        # =================================

        prev_pre_nb = pre_nb

    # =====================================
    # Completed
    # =====================================

    print(
        "Iterative preprocessor+solve finished."
    )


# =========================================
# GET PREPROCESS HISTORY
# =========================================

def get_preprocessing_history():

    return PREPROC_HISTORY
