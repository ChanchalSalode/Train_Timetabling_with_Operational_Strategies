# constraints.py

from gurobipy import quicksum


def add_operation_planning_constraints(
    mdl,
    vars,
    sets,
    params,
    data
):

    # -----------------------------
    # Variables
    # -----------------------------

    tau = vars['tau']
    z = vars['z']
    y = vars['y']
    x = vars['x']
    a1 = vars['a1']
    d = vars['d']
    h = vars['h']

    alpha_u1 = vars['alpha_u1']
    alpha_d2 = vars['alpha_d2']
    alpha_u3 = vars['alpha_u3']
    alpha_d4 = vars['alpha_d4']

    beta_u2 = vars['beta_u2']
    beta_d1 = vars['beta_d1']
    beta_u4 = vars['beta_u4']
    beta_d3 = vars['beta_d3']

    RS1 = vars['RS1']
    RS2 = vars['RS2']
    RS3 = vars['RS3']
    RS4 = vars['RS4']

    w = vars['w']
    w_b = vars['w_b']
    w_b1 = vars['w_b1']

    n_b = vars['n_b']
    n_b1 = vars['n_b1']
    n_a = vars['n_a']
    n1 = vars['n1']

    v = vars['v']
    sai = vars['sai']
    # -----------------------------
    # Sets
    # -----------------------------

    K_u = sets['K_u']
    K_d = sets['K_d']

    S_u = sets['S_u']
    S_d = sets['S_d']

    B = sets['B']

    # -----------------------------
    # Parameters
    # -----------------------------

    h_min = params['h_min']
    h_max = params['h_max']

    delta_min = params['delta_min']

    RS = params['RS']

    C = params['C']

    M = params['M']
    M_time = params['M_time']

    # -----------------------------
    # Data
    # -----------------------------

    r = data['r']
    e = data['e']
    p = data['p']

    # =========================================================
    # OPERATION PLANNING CONSTRAINTS
    # =========================================================

    #Operation Planning Constraints
    mdl.addConstrs(z[k,1,6] + z[k,1,8] + z[k,3,6]+ z[k,3,8]  == tau[k] for k in K_u)
    mdl.addConstrs(z[k,9,14] + z[k,9,16] + z[k,11,14] + z[k,11,16]  == tau[k] for k in K_d)
    mdl.addConstrs(x[k,i] <= tau[k] for i in S_u for k in K_u)
    mdl.addConstrs(x[k,i] <= tau[k] for i in S_d for k in K_d)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M*(1-z[k,1,6])for k in K_u)
    mdl.addConstrs(x[k,i] >= z[k,1,6] for k in K_u for i in S_u if i >= 1 if i <= 6)
    mdl.addConstrs(x[k,i] >= z[k,1,8] for k in K_u for i in S_u if i >= 1 if i <= 8)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M*(1-z[k,3,6])for k in K_u)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M*(1-z[k,3,6])for k in K_u)
    mdl.addConstrs(x[k,i] >= z[k,3,6] for k in K_u for i in S_u if i >= 3 if i <= 6)
    mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M*(1-z[k,3,8])for k in K_u)
    mdl.addConstrs(x[k,i] >= z[k,3,8] for k in K_u for i in S_u if i >= 3 if i <= 8)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M*(1-z[k,9,14])for k in K_d)
    mdl.addConstrs(x[k,i] >= z[k,9,14] for k in K_d for i in S_d if i >= 9 if i <= 14)
    mdl.addConstrs(x[k,i] >= z[k,9,16] for k in K_d for i in S_d if i >= 9 if i <= 16)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M*(1-z[k,11,14])for k in K_d)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M*(1-z[k,11,14])for k in K_d)
    mdl.addConstrs(x[k,i] >= z[k,11,14] for k in K_d for i in S_d if i >= 11 if i <= 14)
    mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M*(1-z[k,11,16])for k in K_d)
    mdl.addConstrs(x[k,i] >= z[k,11,16] for k in K_d for i in S_d if i >= 11 if i <= 16)
    mdl.addConstrs(x[k-1,i] + x[k,i] >= 1 for k in K_u for i in S_u if k != 1)
    mdl.addConstrs(x[k-1,i] + x[k,i] >= 1 for k in K_d for i in S_d if k != 7)
    mdl.addConstr(d[1,1] == 27000)
    mdl.addConstr(d[7,9] == 27000)
    mdl.addConstrs(d[k,1] <= 28800 for k in K_u)
    mdl.addConstrs(d[k,3] <= 28800 for k in K_u)
    mdl.addConstrs(d[k,9] <= 28800 for k in K_d)
    mdl.addConstrs(d[k,11] <= 28800 for k in K_d)
    mdl.addConstrs(a1[k,i] - d[k,i-1] == r[i-1,i] for k in K_u for i in S_u if i != 1)
    mdl.addConstrs(a1[k,i] - d[k,i-1] == r[i-1,i] for k in K_d for i in S_d if i != 9)
    mdl.addConstrs(d[k,i] - a1[k,i] == e[i] for k in K_u for i in S_u)
    mdl.addConstrs(d[k,i] - a1[k,i] == e[i] for k in K_d for i in S_d)
    mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_u if k!= 1)
    mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_u if k!= 1)
    mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_d if k!= 7)
    mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_d if k!= 7)
    mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_u for i in S_u if k != 1)
    mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_d for i in S_d if k != 7)
    mdl.addConstrs(2*y[k,l,6] <= z[k,1,6] +z[l,11,16] + z[k,3,6] + z[l,11,14] for k in K_u for l in K_d)
    mdl.addConstrs(2*y[k,l,8] <= z[k,1,8] + z[l,9,16] +z[k,3,8] + z[l,9,14] for k in K_u for l in K_d)
    mdl.addConstrs(2*y[k,l,14] <= z[k,9,14] + z[l,3,8]+z[k,11,14] + z[l,3,6] for k in K_d for l in K_u)
    mdl.addConstrs(2*y[k,l,16] <= z[k,9,16] + z[l,1,8]+z[k,11,16] + z[l,1,6] for k in K_d for l in K_u)
    mdl.addConstrs(a1[l,9] - d[k,8] >= delta_min + M_time*(y[k,l,8]-1) for k in K_u for l in K_d)
    mdl.addConstrs(a1[l,11] - d[k,6] >= delta_min + M_time*(y[k,l,6]-1) for k in K_u for l in K_d)
    mdl.addConstrs(a1[l,1] - d[k,16] >= delta_min + M_time*(y[k,l,16]-1) for k in K_d for l in K_u)
    mdl.addConstrs(a1[l,3] - d[k,14] >= delta_min + M_time*(y[k,l,14]-1) for k in K_d for l in K_u)
    mdl.addConstrs(quicksum(y[l,k,16]for l in K_d)+quicksum(y[l,k,14]for l in K_d) + alpha_u1[k] + alpha_u3[k] == tau[k] for k in K_u)
    mdl.addConstrs(quicksum(y[l,k,6]for l in K_u)+quicksum(y[l,k,8]for l in K_u) + alpha_d2[k] + alpha_d4[k] == tau[k] for k in K_d)
    mdl.addConstrs(quicksum(y[k,l,6]for l in K_d)+quicksum(y[k,l,8]for l in K_d) + beta_u2[k]+ beta_u4[k] == tau[k] for k in K_u)
    mdl.addConstrs(quicksum(y[k,l,16]for l in K_u)+quicksum(y[k,l,14]for l in K_u) + beta_d1[k]+ beta_d3[k] == tau[k] for k in K_d)
    mdl.addConstr(quicksum(alpha_u1[k] for k in K_u) + quicksum(alpha_d2[k] for k in K_d)+ quicksum(alpha_u3[k] for k in K_u) + quicksum(alpha_d4[k] for k in K_d)<= RS)
    mdl.addConstr(quicksum(alpha_u1[k] for k in K_u) <= RS1)
    mdl.addConstr(quicksum(alpha_d2[k] for k in K_d) <= RS2)
    mdl.addConstr(quicksum(alpha_u3[k] for k in K_u) <= RS3)
    mdl.addConstr(quicksum(alpha_d4[k] for k in K_d) <= RS4)
    mdl.addConstrs(alpha_u1[k] <= z[k,1,8] + z[k,1,6] for k in K_u)
    mdl.addConstrs(alpha_d2[k] <= z[k,9,16] + z[k,9,14] for k in K_d)
    mdl.addConstrs(beta_d1[k] <= z[k,9,16]+ z[k,11,16]for k in K_d)
    mdl.addConstrs(beta_u2[k] <= z[k,1,8] + z[k,3,8] for k in K_u)
    mdl.addConstrs(alpha_u3[k] <= z[k,3,8] + z[k,3,6] for k in K_u)
    mdl.addConstrs(alpha_d4[k] <= z[k,11,16] + z[k,11,14] for k in K_d)
    mdl.addConstrs(beta_d3[k] <= z[k,9,14]+z[k,11,14] for k in K_d)
    mdl.addConstrs(beta_u4[k] <= z[k,1,6] + z[k,3,6] for k in K_u)

    #Passenger Demand Constraints (Time Invarient)
    mdl.addConstrs(w[1,i,j] == (p[i,j]/1800)*120 for i in S_u for j in S_u if i < j)
    mdl.addConstrs(w[k,i,j] == v[k-1,i,j] + (p[i,j]/1800)*(d[k,i] - d[k-1,i]) for k in K_u for i in S_u for j in S_u if i < j if k != 1)
    mdl.addConstrs(w[7,i,j] == (p[i,j]/1800)*120 for i in S_d for j in S_d if i < j)
    mdl.addConstrs(w[k,i,j] == v[k-1,i,j] + (p[i,j]/1800)*(d[k,i] - d[k-1,i]) for k in K_d for i in S_d for j in S_d if i < j if k != 7)
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
    mdl.addConstrs(n1[k, 6] <= C * (1 - z[k, 1, 6]) for k in K_u)
    mdl.addConstrs(n1[k, 6] <= C * (1 - z[k, 3, 6]) for k in K_u)
    mdl.addConstrs(n1[k, 8] <= C * (1 - z[k, 1, 8]) for k in K_u)
    mdl.addConstrs(n1[k, 8] <= C * (1 - z[k, 3, 8]) for k in K_u)
    mdl.addConstrs(n1[k, 14] <= C * (1 - z[k, 9, 14]) for k in K_d)
    mdl.addConstrs(n1[k, 14] <= C * (1 - z[k, 11, 14]) for k in K_d)
    mdl.addConstrs(n1[k, 16] <= C * (1 - z[k, 9, 16]) for k in K_d)
    mdl.addConstrs(n1[k, 16] <= C * (1 - z[k, 11, 16]) for k in K_d)
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
