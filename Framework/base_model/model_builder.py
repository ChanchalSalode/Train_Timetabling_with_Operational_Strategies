# model_builder.py

from gurobipy import Model

from sets import create_sets
from variables import add_variables
from constraints import add_constraints
from objective import set_objective


def build_model(params):

    mdl = Model(
        'Model_1_Minimizing_Opearation_Cost'
    )

    sets = create_sets()

    vars = add_variables(
        mdl,
        sets
    )

    add_constraints(
        mdl,
        vars,
        sets,
        params
    )

    set_objective(
        mdl,
        vars,
        sets
    )

    return mdl, vars, sets
