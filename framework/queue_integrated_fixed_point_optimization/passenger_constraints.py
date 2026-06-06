# passenger_constraints.py

from gurobipy import quicksum


def add_passenger_constraints(
    mdl,
    vars,
    sets,
    params
):

    # =========================================
    # Variables
    # =========================================

    x = vars['x']

    z = vars['z']

    w = vars['w']

    w_b = vars['w_b']
    w_b1 = vars['w_b1']

    n_b = vars['n_b']
    n_b1 = vars['n_b1']

    n_a = vars['n_a']

    n1 = vars['n1']

    v = vars['v']

    sai = vars['sai']

    # =========================================
    # Sets
    # =========================================

    K_u = sets['K_u']
    K_d = sets['K_d']

    S_u = sets['S_u']
    S_d = sets['S_d']

    # =========================================
    # Parameters
    # =========================================

    C = params['C']

    M = params['M']

    # =========================================
    # Boarding Constraints
    # =========================================

    mdl.addConstrs(w_b[k,i] == quicksum(w_b1[k,i,j] for j in S_u if i < j) for k in K_u for i in S_u)
    mdl.addConstrs(w_b[k,i] == quicksum(w_b1[k,i,j] for j in S_d if i < j) for k in K_d for i in S_d)
    mdl.addConstrs(n_a[k,i] == quicksum(n_b1[k,j,i] for j in S_u if j < i) for k in K_u for i in S_u)
    mdl.addConstrs(n_a[k,i] == quicksum(n_b1[k,j,i] for j in S_d if j < i) for k in K_d for i in S_d)
    mdl.addConstrs(n_b[k,i] - n_b1[k,i,j] == w_b[k,i] - w_b1[k,i,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(n_b[k,i] - n_b1[k,i,j] == w_b[k,i] - w_b1[k,i,j] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= 0 for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= w[k,i,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,i] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M*(2-x[k,i]-x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= 0 for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= w[k,i,j] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,i] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] <= M*x[k,j] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M*(2-x[k,i]-x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_u for i in S_u)
    mdl.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_u for i in S_u if i!= 1)
    mdl.addConstrs(n_b[k,i] >= w_b[k,i] - M*(1-sai[k,i]) for k in K_u for i in S_u)
    mdl.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M*sai[k,i] for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_d for i in S_d)
    mdl.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_d for i in S_d if i!= 9)
    mdl.addConstrs(n_b[k,i] >= w_b[k,i] - M*(1-sai[k,i]) for k in K_d for i in S_d)
    mdl.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M*sai[k,i] for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(n1[k,6] <= C * (1 - z[k,1,6]) for k in K_u)
    mdl.addConstrs(n1[k,6] <= C * (1 - z[k,3,6]) for k in K_u)
    mdl.addConstrs(n1[k,8] <= C * (1 - z[k,1,8]) for k in K_u)
    mdl.addConstrs(n1[k,8] <= C * (1 - z[k,3,8]) for k in K_u)
    mdl.addConstrs(n1[k,14] <= C * (1 - z[k,9,14]) for k in K_d)
    mdl.addConstrs(n1[k,14] <= C * (1 - z[k,11,14]) for k in K_d)
    mdl.addConstrs(n1[k,16] <= C * (1 - z[k,9,16]) for k in K_d)
    mdl.addConstrs(n1[k,16] <= C * (1 - z[k,11,16]) for k in K_d)
    mdl.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(n1[k,1] == n_b[k,1] for k in K_u)
    mdl.addConstrs(n1[k,9] == n_b[k,9] for k in K_d)
    mdl.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M*(2 - x[k,i] - x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M*(2 - x[k,i] - x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,i] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,j] for k in K_u for i in S_u for j in S_u if i < j)
    mdl.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M*(2 - x[k,i] - x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M*(2 - x[k,i] - x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,i] for k in K_d for i in S_d for j in S_d if i < j)
    mdl.addConstrs(v[k,i,j] <= M*x[k,j] for k in K_d for i in S_d for j in S_d if i < j)