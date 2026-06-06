# solve.py

from config import (
    MU_BOARD,
    DWELL_BASE_SEC,
    DWELL_PER_PAX_SEC,
    PREPROC_MAX_ITERS,
    PREPROC_TOL,
    PREPROC_MODE
)

from iterative_coupling import (
    PREPROC_HISTORY,
    compute_and_apply_pre_penalty,
    run_iterative_with_preprocessor
)


def optimize_model(model, label='model'):
    """Run a Gurobi optimization and keep solve calls out of orchestration code."""

    print(f"\nSolving {label}...")
    model.optimize()
    print(f"{label.capitalize()} solve completed with status {model.Status}.")
    return model.Status


def run_preprocessor_solve(
    model,
    arrivals_by_station_od,
    services_df,
    var_maps,
    initial_backlogs,
    pre_penalty=None,
    max_iters=PREPROC_MAX_ITERS,
    tol=PREPROC_TOL,
    mode=PREPROC_MODE,
    update_services_function=None
):
    """Run the queue preprocessor and repeated MIP solves."""

    if pre_penalty is None:

        pre_penalty = compute_and_apply_pre_penalty(
            model,
            PREPROC_HISTORY,
            alpha=0.1,
            min_penalty=0.1,
            beta=1.0
        )

    return run_iterative_with_preprocessor(
        model=model,
        all_stations=sorted(arrivals_by_station_od.keys()),
        arrivals_by_station_od=arrivals_by_station_od,
        services_csv_or_df=services_df,
        var_maps=var_maps,
        max_iters=max_iters,
        tol=tol,
        mu_board=MU_BOARD,
        dwell_base_sec=DWELL_BASE_SEC,
        dwell_per_pax_sec=DWELL_PER_PAX_SEC,
        initial_backlogs=initial_backlogs,
        mode=mode,
        pre_penalty=pre_penalty,
        update_services_function=update_services_function
    )
