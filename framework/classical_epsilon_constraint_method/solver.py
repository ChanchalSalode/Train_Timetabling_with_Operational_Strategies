# solver.py

from epsilon_loop import run_epsilon_loop


def solve_model(
    mdl,
    vars,
    sets,
    params
):

    run_epsilon_loop(
        mdl=mdl,
        vars=vars,
        sets=sets,
        params=params
    )
