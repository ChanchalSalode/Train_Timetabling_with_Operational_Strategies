# epsilon_loop.py

from gurobipy import GRB, quicksum

from performance_metrics import compute_performance_metrics


def run_epsilon_loop(
    mdl,
    vars,
    sets,
    params
):

    # ============================================
    # Variables
    # ============================================

    y = vars['y']
    tau = vars['tau']
    h = vars['h']
    x = vars['x']
    w = vars['w']
    w_b = vars['w_b']
    v = vars['v']
    n1 = vars['n1']

    # ============================================
    # Sets
    # ============================================

    A = sets['A']
    B = sets['B']

    K_u = sets['K_u']
    K_d = sets['K_d']

    # ============================================
    # Parameters
    # ============================================

    C = params['C']
    h_min = params['h_min']

    # ============================================
    # Initial Solve
    # ============================================

    mdl.setParam(
        'OutputFlag',
        params['OUTPUT_FLAG']
    )

    mdl.optimize()

    solution_count = 1

    # ============================================
    # Epsilon Loop
    # ============================================

    while mdl.Status == GRB.OPTIMAL:

        print("\n====================================")
        print(f"SOLUTION {solution_count}")
        print("====================================")

        # =====================================
        # Store y values
        # =====================================

        y_val = {}

        for a in A:

            y_val[a] = y[a].X

            if y[a].X > 0:

                print(
                    f"y[{a}] = {y[a].X}"
                )

        # =====================================
        # Objective Value
        # =====================================

        print("\nObjective Function Value:")

        print(
            f"ObjVal: {mdl.ObjVal}"
        )

        # =====================================
        # Epsilon Objective Expression
        # =====================================

        kk_expr = (
            quicksum(
                y_val[l,k,6]
                for l in K_u
                for k in K_d
            )
            +
            quicksum(
                y_val[l,k,8]
                for l in K_u
                for k in K_d
            )
            +
            quicksum(
                y_val[l,k,14]
                for l in K_d
                for k in K_u
            )
            +
            quicksum(
                y_val[l,k,16]
                for l in K_d
                for k in K_u
            )
        )

        print(
            f"Value of kk: {kk_expr.getValue()}"
        )

        # =====================================
        # Performance Metrics
        # =====================================

        compute_performance_metrics(
            vars=vars,
            sets=sets,
            params=params
        )

        # =====================================
        # Add Epsilon Constraint
        # =====================================

        mdl.addConstr(

            quicksum(
                y[l,k,6]
                for l in K_u
                for k in K_d
            )
            +
            quicksum(
                y[l,k,8]
                for l in K_u
                for k in K_d
            )
            +
            quicksum(
                y[l,k,14]
                for l in K_d
                for k in K_u
            )
            +
            quicksum(
                y[l,k,16]
                for l in K_d
                for k in K_u
            )

            >= kk_expr + 1

        )

        # =====================================
        # Resolve
        # =====================================

        mdl.optimize()

        solution_count += 1
