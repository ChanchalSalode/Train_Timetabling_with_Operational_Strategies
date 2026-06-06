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

    # =========================================
    # Alpha Variables
    # =========================================

    alpha_up1 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='alpha_up1'
    )

    alpha_dn2 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='alpha_dn2'
    )

    alpha_up3 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='alpha_up3'
    )

    alpha_dn4 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='alpha_dn4'
    )

    # =========================================
    # Beta Variables
    # =========================================

    beta_up2 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='beta_up2'
    )

    beta_dn1 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='beta_dn1'
    )

    beta_up4 = mdl.addVars(
        K_u,
        vtype=GRB.BINARY,
        name='beta_up4'
    )

    beta_dn3 = mdl.addVars(
        K_d,
        vtype=GRB.BINARY,
        name='beta_dn3'
    )

    # =========================================
    # Rolling Stock Variables
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
    # Epsilon Variable
    # =========================================

    y_level = mdl.addVar(
        lb=0,
        vtype=GRB.INTEGER,
        name='y_level'
    )

    # =========================================
    # Linearization Variable
    # =========================================

    zz = mdl.addVars(
        B,
        vtype=GRB.CONTINUOUS,
        name='zz'
    )

    # =========================================
    # Monotonicity Constraints
    # =========================================

    mdl.addConstrs(

        tau[k] >= tau[k + 1]

        for k in K_dd

        if (k + 1) in K_dd

    )

    # =========================================
    # Pareto Constraint
    # =========================================

    mdl.addConstr(

        sum(
            y[k, l, 6]

            for k in K_u
            for l in K_d

        )

        +

        sum(
            y[k, l, 8]

            for k in K_u
            for l in K_d

        )

        +

        sum(
            y[k, l, 14]

            for k in K_d
            for l in K_u

        )

        +

        sum(
            y[k, l, 16]

            for k in K_d
            for l in K_u

        )

        >= y_level,

        name='Pareto_y_level'

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

        'alpha_up1': alpha_up1,
        'alpha_dn2': alpha_dn2,
        'alpha_up3': alpha_up3,
        'alpha_dn4': alpha_dn4,

        'beta_up2': beta_up2,
        'beta_dn1': beta_dn1,
        'beta_up4': beta_up4,
        'beta_dn3': beta_dn3,

        'RS1': RS1,
        'RS2': RS2,
        'RS3': RS3,
        'RS4': RS4,

        'y_level': y_level,

        'zz': zz,

    }
