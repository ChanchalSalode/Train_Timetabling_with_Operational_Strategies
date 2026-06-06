# callback.py

from gurobipy import GRB, quicksum

from subproblem import solve_passenger_subproblem


# =========================================
# FEASIBILITY CACHE
# =========================================

feasibility_cache = {}


# =========================================
# BENDERS CALLBACK
# =========================================

def benders_callback(
    model,
    where,
    vars,
    sets
):

    # =====================================
    # Trigger only at integer solutions
    # =====================================

    if where == GRB.Callback.MIPSOL:

        # =================================
        # Variables
        # =================================

        x = vars['x']

        d = vars['d']

        z = vars['z']

        # =================================
        # Sets
        # =================================

        C1 = sets['C1']

        B = sets['B']

        # =================================
        # Extract incumbent solution
        # =================================

        x_val = {

            (k, i):

            model.cbGetSolution(
                x[k, i]
            )

            for (k, i)
            in C1

        }

        d_val = {

            (k, i):

            model.cbGetSolution(
                d[k, i]
            )

            for (k, i)
            in C1

        }

        z_val = {

            (k, m, n):

            model.cbGetSolution(
                z[k, m, n]
            )

            for (k, m, n)
            in B

        }

        # =================================
        # Create pattern key
        # =================================

        pattern_key = tuple(

            (k, i)

            for (k, i)
            in C1

            if x_val[(k, i)] > 0.5

        )

        # =================================
        # Cache Check
        # =================================

        if pattern_key in feasibility_cache:

            # =============================
            # Infeasible pattern
            # =============================

            if not feasibility_cache[pattern_key]:

                violated_keys = list(
                    pattern_key
                )

                model.cbLazy(

                    quicksum(

                        x[k, i]

                        for (k, i)
                        in violated_keys

                    )

                    <=

                    len(violated_keys) - 1

                )

            return

        # =================================
        # Solve Passenger Subproblem
        # =================================

        sub_vars = solve_passenger_subproblem(
            sets=sets,
            params=model._params_data,
            data=model._input_data,
            x_val=x_val,
            d_val=d_val,
            z_val=z_val
        )

        # =================================
        # Store feasibility result
        # =================================

        feasibility_cache[pattern_key] = (

            sub_vars is not None

        )

        # =================================
        # If feasible → no cut
        # =================================

        if sub_vars is not None:

            return

        # =================================
        # Generate no-good cut
        # =================================

        violated_keys = [

            (k, i)

            for (k, i)
            in C1

            if x_val[(k, i)] > 0.5

        ]

        model.cbLazy(

            quicksum(

                x[k, i]

                for (k, i)
                in violated_keys

            )

            <=

            len(violated_keys) - 1

        )