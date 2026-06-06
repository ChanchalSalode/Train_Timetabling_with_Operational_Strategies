# variables.py

from gurobipy import GRB
from config import *


def add_variables(mdl, sets):

    K_dd = sets['K_dd']
    K_u = sets['K_u']
    K_d = sets['K_d']

    B = sets['B']
    A = sets['A']
    C1 = sets['C1']
    D = sets['D']

    tau = mdl.addVars(K_dd, vtype=GRB.BINARY, name='tau')

    z = mdl.addVars(B, vtype=GRB.BINARY, name='z')

    y = mdl.addVars(A, vtype=GRB.BINARY, name='y')

    x = mdl.addVars(C1, vtype=GRB.BINARY, name='x')

    a1 = mdl.addVars(C1, vtype=GRB.CONTINUOUS, name='a1')

    d = mdl.addVars(C1, vtype=GRB.CONTINUOUS, name='d')

    h = mdl.addVars(D, vtype=GRB.CONTINUOUS, name='h')

    alpha_up1 = mdl.addVars(K_u, vtype=GRB.BINARY)

    alpha_dn2 = mdl.addVars(K_d, vtype=GRB.BINARY)

    alpha_up3 = mdl.addVars(K_u, vtype=GRB.BINARY)

    alpha_dn4 = mdl.addVars(K_d, vtype=GRB.BINARY)

    beta_up2 = mdl.addVars(K_u, vtype=GRB.BINARY)

    beta_dn1 = mdl.addVars(K_d, vtype=GRB.BINARY)

    beta_up4 = mdl.addVars(K_u, vtype=GRB.BINARY)

    beta_dn3 = mdl.addVars(K_d, vtype=GRB.BINARY)

    RS1 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS1')

    RS2 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS2')

    RS3 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS3')

    RS4 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS4')

    w = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='w')

    w_b = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='w_b')

    w_b1 = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='w_b1')

    n_b = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='n_b')

    n_b1 = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='n_b1')

    n_a = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='n_a')

    n1 = mdl.addVars(
        C1,
        vtype=GRB.CONTINUOUS,
        lb=0,
        ub=C,
        name='n1'
    )

    v = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='v')

    sai = mdl.addVars(C1, vtype=GRB.BINARY, name='sai')

    return locals()