# model_builder.py

from gurobipy import (
    Model,
    GRB,
    quicksum
)

from variables import create_variables

from master_constraints import (
    add_master_constraints
)

from linearization_constraints import (
    add_linearization_constraints
)


# =========================================
# BUILD HYBRID EPSILON MODEL
# =========================================

def build_model(
    sets,
    data,
    params
):

    # =====================================
    # Create Model
    # =====================================

    mdl = Model(
        params['MODEL_NAME']
    )

    # =====================================
    # Enable Lazy Constraints
    # =====================================

    mdl.Params.LazyConstraints = 1

    # =====================================
    # Create Variables
    # =====================================

    # vars = create_variables(
    #     mdl,
    #     sets
    # )
    vars = create_variables(
        mdl=mdl,
        sets=sets,
        params=params
    )

    # =====================================
    # Add Master Constraints
    # =====================================

    add_master_constraints(
        mdl=mdl,
        vars=vars,
        sets=sets,
        data=data,
        params=params
    )

    mdl.setParam(
        'OutputFlag',
        params['OUTPUT_FLAG']
    )

    mdl.setParam(
        'LogFile',
        params['GUROBI_LOG_FILE']
    )
    # =====================================
    # Add Linearization Constraints
    # =====================================

    add_linearization_constraints(
        mdl,
        vars,
        sets,
        params
    )

    # =====================================
    # Variables
    # =====================================

    zz = vars['zz']

    d = vars['d']

    # =====================================
    # Sets
    # =====================================

    K_u = sets['K_u']

    K_d = sets['K_d']

    S_u = sets['S_u']

    S_d = sets['S_d']

    # =====================================
    # Objective Function
    # =====================================

    objective = (

        # ================================
        # Up Direction
        # ================================

        quicksum(
            zz[k, 1, 6]
            for k in K_u
        )

        +

        quicksum(
            zz[k, 1, 8]
            for k in K_u
        )

        +

        quicksum(
            zz[k, 3, 6]
            for k in K_u
        )

        +

        quicksum(
            zz[k, 3, 8]
            for k in K_u
        )

        +

        # ================================
        # Down Direction
        # ================================

        quicksum(
            zz[k, 9, 14]
            for k in K_d
        )

        +

        quicksum(
            zz[k, 9, 16]
            for k in K_d
        )

        +

        quicksum(
            zz[k, 11, 14]
            for k in K_d
        )

        +

        quicksum(
            zz[k, 11, 16]
            for k in K_d
        )

        +

        # ================================
        # Headway Terms
        # ================================

        quicksum(

            d[k, i] - d[k - 1, i]

            for k in K_u
            for i in S_u

            if k != 1

        )

        +

        quicksum(

            d[k, i] - d[k - 1, i]

            for k in K_d
            for i in S_d

            if k != 10

        )

    )

    # =====================================
    # Set Objective
    # =====================================

    mdl.setObjective(
        objective,
        GRB.MINIMIZE
    )

    # =====================================
    # Return Model
    # =====================================

    return mdl, vars
