# model_builder.py

from gurobipy import Model

from variables import create_variables

from constraints import add_operation_planning_constraints

from linearization_constraints import add_linearization_constraints

from objective import set_objective


def build_model(
    sets,
    params,
    data
):

    # =========================================
    # Create Model
    # =========================================

    mdl = Model(
        params['MODEL_NAME']
    )

    # =========================================
    # Create Variables
    # =========================================

    vars = create_variables(
        mdl=mdl,
        sets=sets,
        params=params
    )

    # =========================================
    # Add Constraints
    # =========================================

    add_operation_planning_constraints(
        mdl=mdl,
        vars=vars,
        sets=sets,
        params=params,
        data=data
    )

    add_linearization_constraints(
        mdl=mdl,
        vars=vars,
        sets=sets,
        params=params
    )

    set_objective(
        mdl=mdl,
        vars=vars,
        sets=sets
    )

    # =========================================
    # Silent Output
    # =========================================

    mdl.setParam(
        'OutputFlag',
        params['OUTPUT_FLAG']
    )

    mdl.setParam(
        'LogToConsole',
        params['LOG_TO_CONSOLE']
    )

    mdl.setParam(
        'LogFile',
        params['GUROBI_LOG_FILE']
    )

    # =========================================
    # Return Model and Variables
    # =========================================

    return mdl, vars
