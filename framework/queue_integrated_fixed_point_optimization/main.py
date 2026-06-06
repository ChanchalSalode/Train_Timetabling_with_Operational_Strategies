# main.py

from config import (
    S_U,
    S_D,
    C,
    H,
    T_WINDOW,
    _times,
    _n_bins
)

from input_loader import (
    load_running_time,
    load_dwell_time,
    load_demand_matrix
)

from backlog_loader import (
    get_backlog_matrix
)

from demand_processing import (
    build_phi_from_excel,
    build_pwl_breaks_for_dir,
    build_arrivals_by_station,
    load_od_share_per_bin
)

from model_builder import (
    build_trb_model
)

from metrics import (
    collect_performance_metrics_from_model
)

from service_generator import (
    generate_initial_services,
    update_services_from_solution
)

from solve import (
    optimize_model,
    run_preprocessor_solve
)

from diagnostics import (
    write_metrics_comparison
)

# ============================================================
# MAIN EXECUTION
# ============================================================

def main():

    print("\n====================================")
    print("QFPO MODULAR FRAMEWORK")
    print("====================================")

    # ========================================================
    # LOAD INPUT DATA
    # ========================================================

    print("\nLoading input data...")

    r = load_running_time()

    e = load_dwell_time()

    p = load_demand_matrix()

    print("Input data loaded successfully.")

    # ========================================================
    # LOAD ARRIVAL PROFILE
    # ========================================================

    print("\nLoading arrival profiles...")

    phi = build_phi_from_excel()

    od_share_per_bin = load_od_share_per_bin(
        demand=p,
        n_stations=S_D,
        n_bins=_n_bins
    )

    print("Arrival profiles loaded successfully.")

    # ========================================================
    # LOAD BACKLOG
    # ========================================================

    print("\nLoading backlog matrix...")

    P_backlog = get_backlog_matrix(p)

    print("Backlog matrix ready.")

    # ========================================================
    # BUILD PWL ARRIVALS
    # ========================================================

    print("\nBuilding PWL demand profiles...")

    # pwl_up = build_pwl_breaks_for_dir(
    #     p,
    #     phi,
    #     list(range(1, S_U + 1))
    # )

    # pwl_down = build_pwl_breaks_for_dir(
    #     p,
    #     phi,
    #     list(range(9, S_D + 1))
    # )

    # =========================================
    # PWL SETTINGS
    # =========================================

    times = _times

    n_bins = _n_bins

    # =========================================
    # UP DIRECTION PWL
    # =========================================

    pwl_up = build_pwl_breaks_for_dir(

        p,

        phi,

        list(range(1, S_U + 1)),

        times,

        n_bins,

        od_share_per_bin=od_share_per_bin

    )

    # =========================================
    # DOWN DIRECTION PWL
    # =========================================

    pwl_down = build_pwl_breaks_for_dir(

        p,

        phi,

        list(range(9, S_D + 1)),

        times,

        n_bins,

        od_share_per_bin=od_share_per_bin

    )

    print("PWL demand profiles created.")

    # ========================================================
    # ARRIVAL SERIES
    # ========================================================

    arrivals_by_station_od = build_arrivals_by_station(
        pwl_up,
        pwl_down
    )

    # ========================================================
    # BUILD MODEL
    # ========================================================

    print("\nBuilding optimization model...")

    mdl, var_dicts = build_trb_model(
        p,
        r,
        e,
        backlog=P_backlog,
        pwl_up=pwl_up,
        pwl_down=pwl_down
    )

    print("Model built successfully.")

    # ========================================================
    # INITIAL OPTIMIZATION
    # ========================================================

    optimize_model(
        mdl,
        label='initial model'
    )

    # ========================================================
    # EXTRACT VARIABLES
    # ========================================================

    # ========================================================
    # CREATE SERVICE TABLE
    # ========================================================

    print("\nCreating service table...")

    services_df = generate_initial_services(
        var_dicts,
        H,
        C
    )

    print("Service table created.")

    # ========================================================
    # INITIAL BACKLOGS
    # ========================================================

    initial_backlogs_per_station = {}

    for st in sorted(
        arrivals_by_station_od.keys()
    ):

        sub = {

            od: val

            for od, val
            in P_backlog.items()

            if od[0] == st

        }

        initial_backlogs_per_station[st] = sub

    # ========================================================
    # VARIABLE MAPS
    # ========================================================

    var_maps = {

        'w_var': var_dicts.get(
            'w',
            {}
        ),

        'nb_var': var_dicts.get(
            'n_b1',
            {}
        ),

        'nb_total_var': var_dicts.get(
            'n_b',
            {}
        )

    }

    # ========================================================
    # METRICS BEFORE ITERATION
    # ========================================================

    print("\nCollecting baseline metrics...")

    metrics_before = (

        collect_performance_metrics_from_model(
            mdl,
            var_dicts
        )

    )

    # ========================================================
    # ITERATIVE PREPROCESSOR
    # ========================================================

    print("\n[Preprocessor] Running iterative coupling (max_iters=3)...")

    def update_services(_model, services_src):

        update_services_from_solution(
            services_src,
            var_dicts,
            H,
            T_WINDOW
        )

    run_preprocessor_solve(
        model=mdl,
        arrivals_by_station_od=arrivals_by_station_od,
        services_df=services_df,
        var_maps=var_maps,
        initial_backlogs=initial_backlogs_per_station,
        update_services_function=update_services
    )

    print("[Preprocessor] Iterative coupling complete.")

    # ========================================================
    # METRICS AFTER ITERATION
    # ========================================================

    print("\nCollecting final metrics...")

    metrics_after = (

        collect_performance_metrics_from_model(
            mdl,
            var_dicts
        )

    )

    # ========================================================
    # EXPORT METRICS
    # ========================================================

    print("\nWriting metrics comparison...")

    write_metrics_comparison(

        metrics_before=metrics_before,

        metrics_after=metrics_after,

        output_file=(
            "Qmetrics_comparison.txt"
        )

    )

    # ========================================================
    # FINAL MESSAGE
    # ========================================================

    print("\n====================================")
    print("QFPO MODULAR FRAMEWORK COMPLETED")
    print("====================================")


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    main()
