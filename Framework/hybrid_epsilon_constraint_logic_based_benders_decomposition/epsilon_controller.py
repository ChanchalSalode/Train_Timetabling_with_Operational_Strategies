# epsilon_controller.py

from gurobipy import GRB

from callback import benders_callback

from subproblem import solve_passenger_subproblem

from performance_metrics import analyze_passenger_solution

from config import (
    TIME_LIMIT,
    MAX_Y
)


# =========================================
# EPSILON CONSTRAINT CONTROLLER
# =========================================

def run_epsilon_loop(
    mdl,
    vars,
    sets,
    params,
    data
):

    # =====================================
    # Variables
    # =====================================

    x = vars['x']

    d = vars['d']

    z = vars['z']

    h = vars['h']

    tau = vars['tau']

    y = vars['y']

    y_level = vars['y_level']

    # =====================================
    # Store callback data inside model
    # =====================================

    mdl._params_data = params
    mdl._input_data = data

    # =====================================
    # Initial Optimization
    # =====================================

    print("\n=================================")
    print("INITIAL HYBRID SOLUTION")
    print("=================================")

    mdl.optimize(

        lambda model, where:

        benders_callback(
            model,
            where,
            vars,
            sets
        )

    )

    mdl.update()

    # =====================================
    # Extract Initial y*
    # =====================================

    y_vals = mdl.getAttr(
        "X",
        y
    )

    y_star = int(

        round(

            sum(
                y_vals.values()
            )

        )

    )

    level = y_star

    print(f"\nInitial y* = {y_star}")

    # =====================================
    # Start ε-loop
    # =====================================

    while level <= MAX_Y:

        print("\n====================================")
        print(f"SOLUTION FOR ε = {level}")
        print("====================================")

        # =================================
        # Update epsilon lower bound
        # =================================

        y_level.lb = level

        mdl.update()

        # =================================
        # Set time limit
        # =================================

        mdl.setParam(
            GRB.Param.TimeLimit,
            TIME_LIMIT
        )

        # =================================
        # Optimize
        # =================================

        mdl.optimize(

            lambda model, where:

            benders_callback(
                model,
                where,
                vars,
                sets
            )

        )

        # =================================
        # TIME LIMIT CASE
        # =================================

        if mdl.Status == GRB.TIME_LIMIT:

            if mdl.SolCount > 0:

                print(
                    f"[TIME LIMIT] "
                    f"Extracting incumbent for ε = {level}"
                )

                mdl.update()

                # =========================
                # Extract Master Solution
                # =========================

                x_val = dict(
                    mdl.getAttr("X", x)
                )

                d_val = dict(
                    mdl.getAttr("X", d)
                )

                z_val = dict(
                    mdl.getAttr("X", z)
                )

                h_val = dict(
                    mdl.getAttr("X", h)
                )

                tau_val = dict(
                    mdl.getAttr("X", tau)
                )

                # =========================
                # Solve Passenger Subproblem
                # =========================

                sub_vars = solve_passenger_subproblem(

                    sets=sets,

                    params=params,

                    data=data,

                    x_val=x_val,

                    d_val=d_val,

                    z_val=z_val

                )

                # =========================
                # Analyze Solution
                # =========================

                if sub_vars is not None:

                    analyze_passenger_solution(

                        solution_id=level,

                        tau_val=tau_val,

                        x_val=x_val,

                        d_val=d_val,

                        z_val=z_val,

                        h_val=h_val,

                        sub_vars=sub_vars,

                        output_prefix="TIME_LIMIT"

                    )

            else:

                print(
                    f"[TIME LIMIT] "
                    f"No incumbent for ε = {level}"
                )

            level += 1

            continue

        # =================================
        # Infeasible Case
        # =================================

        if (
            mdl.Status == GRB.INFEASIBLE
            or
            mdl.SolCount == 0
        ):

            print(
                f"No feasible solution "
                f"for ε = {level}"
            )

            level += 1

            continue

        # =================================
        # Feasible Solution Found
        # =================================

        print(
            f"Hybrid ε+Benders solution "
            f"found for y_level = {level}"
        )

        # =================================
        # Extract Master Variables
        # =================================

        mdl.update()

        x_val = dict(
            mdl.getAttr("X", x)
        )

        d_val = dict(
            mdl.getAttr("X", d)
        )

        z_val = dict(
            mdl.getAttr("X", z)
        )

        h_val = dict(
            mdl.getAttr("X", h)
        )

        tau_val = dict(
            mdl.getAttr("X", tau)
        )

        # =================================
        # Solve Passenger Subproblem
        # =================================

        sub_vars = solve_passenger_subproblem(
            sets=sets,
            params=params,
            data=data,
            x_val=x_val,
            d_val=d_val,
            z_val=z_val
        )

        # =================================
        # Analyze Solution
        # =================================

        if sub_vars is not None:

            analyze_passenger_solution(

                solution_id=level,

                tau_val=tau_val,

                x_val=x_val,

                d_val=d_val,

                z_val=z_val,

                h_val=h_val,

                sub_vars=sub_vars,

                output_prefix="FINAL"

            )

        # =================================
        # Next epsilon level
        # =================================

        level += 1

    # =====================================
    # Completed
    # =====================================

    print("\n=================================")
    print("EPSILON LOOP COMPLETED")
    print("=================================")