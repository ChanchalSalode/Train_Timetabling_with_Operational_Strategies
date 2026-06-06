# # solver.py

from config import TIME_LIMIT
from config import GUROBI_LOG_FILE


def solve_model(mdl):

    mdl.Params.TimeLimit = TIME_LIMIT

    # Save Gurobi log
    mdl.Params.LogFile = GUROBI_LOG_FILE

    mdl.optimize()