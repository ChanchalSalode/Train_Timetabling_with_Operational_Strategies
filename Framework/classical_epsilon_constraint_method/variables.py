# variables.py

from gurobipy import GRB


def create_variables(
    mdl,
    sets,
    params
):

    # =========================================
    # Sets
    # =========================================

    K_u = sets['K_u']
    K_d = sets['K_d']

    K_dd = sets['K_dd']

    B = sets['B']
    A = sets['A']
    C1 = sets['C1']
    D = sets['D']

    # =========================================
    # Parameters
    # =========================================

    C = params['C']

    # =========================================
    # Binary Variables
    # =========================================

    tau = mdl.addVars(
        K_dd,
        vtype=GRB.BINARY,
        name='tau'
    )

    z = mdl.addVars(
        B,
        vtype=GRB.BINARY,
        name='z'
    )

    y = mdl.addVars(
        A,
        vtype=GRB.BINARY,
        name='y'
    )

    x = mdl.addVars(
        C1,
        vtype=GRB.BINARY,
        name='x'
    )

    sai = mdl.addVars(
        C1,
        vtype=GRB.BINARY,
        name='sai'
    )

    alpha_u1 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='alpha_u1'
    )

    alpha_d2 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='alpha_d2'
    )

    alpha_u3 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='alpha_u3'
    )

    alpha_d4 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='alpha_d4'
    )

    beta_u2 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='beta_u2'
    )

    beta_d1 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='beta_d1'
    )

    beta_u4 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='beta_u4'
    )

    beta_d3 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='beta_d3'
    )

    # =========================================
    # Continuous Variables
    # =========================================

    a1 = mdl.addVars(
        C1,
        vtype=GRB.CONTINUOUS,
        name='a1'
    )

    d = mdl.addVars(
        C1,
        vtype=GRB.CONTINUOUS,
        name='d'
    )

    h = mdl.addVars(
        D,
        vtype=GRB.CONTINUOUS,
        name='h'
    )

    w = mdl.addVars(
        B,
        vtype=GRB.CONTINUOUS,
        lb=0.0,
        name='w'
    )

    w_b = mdl.addVars(
        C1,
        vtype=GRB.CONTINUOUS,
        lb=0.0,
        name='w_b'
    )

    w_b1 = mdl.addVars(
        B,
        vtype=GRB.CONTINUOUS,
        lb=0.0,
        name='w_b1'
    )

    n_b = mdl.addVars(
        C1,
        vtype=GRB.CONTINUOUS,
        lb=0.0,
        name='n_b'
    )

    n_b1 = mdl.addVars(
        B,
        vtype=GRB.CONTINUOUS,
        lb=0.0,
        name='n_b1'
    )

    n_a = mdl.addVars(
        C1,
        vtype=GRB.CONTINUOUS,
        lb=0.0,
        name='n_a'
    )

    n1 = mdl.addVars(
        C1,
        vtype=GRB.CONTINUOUS,
        lb=0,
        ub=C,
        name='n1'
    )

    v = mdl.addVars(
        B,
        vtype=GRB.CONTINUOUS,
        lb=0.0,
        name='v'
    )

    zz = mdl.addVars(
        B,
        vtype=GRB.CONTINUOUS,
        name='zz'
    )

    # =========================================
    # Integer Variables
    # =========================================

    RS1 = mdl.addVar(
        vtype=GRB.INTEGER,
        lb=0,
        name='RS1'
    )

    RS2 = mdl.addVar(
        vtype=GRB.INTEGER,
        lb=0,
        name='RS2'
    )

    RS3 = mdl.addVar(
        vtype=GRB.INTEGER,
        lb=0,
        name='RS3'
    )

    RS4 = mdl.addVar(
        vtype=GRB.INTEGER,
        lb=0,
        name='RS4'
    )

    # =========================================
    # Return Variables
    # =========================================

    return {

        'tau': tau,
        'z': z,
        'y': y,
        'x': x,

        'a1': a1,
        'd': d,
        'h': h,

        'alpha_u1': alpha_u1,
        'alpha_d2': alpha_d2,

        'alpha_u3': alpha_u3,
        'alpha_d4': alpha_d4,

        'beta_u2': beta_u2,
        'beta_d1': beta_d1,

        'beta_u4': beta_u4,
        'beta_d3': beta_d3,

        'RS1': RS1,
        'RS2': RS2,
        'RS3': RS3,
        'RS4': RS4,

        'w': w,
        'w_b': w_b,
        'w_b1': w_b1,

        'n_b': n_b,
        'n_b1': n_b1,

        'n_a': n_a,
        'n1': n1,

        'v': v,

        'sai': sai,

        'zz': zz

    }