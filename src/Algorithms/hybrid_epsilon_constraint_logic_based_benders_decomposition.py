import pandas as pd
# Model Parameters
S_U = 8
S_D = 16
K_U = 9
K_D = 17
h_max = 360
h_min = 90
delta_min = 135
RS = 14
C = 250
M_time = 100000.0
M_flow = 5000.0
M = M_flow

#Reading data from excel
df1= pd.read_excel('Input_Parameters.xlsx', 'Pure_running_time')
my_tuple = [tuple(x) for x in df1.values]
pr =dict((tuple((a, b)), c) for a,b,c in df1.values)
df1= pd.read_excel('Input_Parameters.xlsx', 'Acceleration_time')
my_tuple = [tuple(x) for x in df1.values]
ac =dict((tuple((a, b)), c) for a,b,c in df1.values)
df1= pd.read_excel('Input_Parameters.xlsx', 'Deceleration_time')
my_tuple = [tuple(x) for x in df1.values]
dc =dict((tuple((a, b)), c) for a,b,c in df1.values)
df1= pd.read_excel('Input_Parameters.xlsx', 'Demand_30(7.30-8.00am)mint_peak')
my_tuple = [tuple(x) for x in df1.values]
p =dict((tuple((a, b)), c) for a,b,c in df1.values)
df1 = pd.read_excel('Input_Parameters.xlsx', 'Dwelling_time')
my_tuple = [tuple(x) for x in df1.values]
e = dict((a, c) for a, c in my_tuple)

#Index Sets
S_d = [i for i in range(9,S_D+1)]
S_dd = [i for i in range(1,S_D+1)]
K_d = [i for i in range(10,K_D+1)]
K_dd = [i for i in range(1,K_D+1)]
K_u = [i for i in range(1,K_U+1)]
S_u = [i for i in range(1,S_U+1)]
B = [(k,m,n) for k in K_dd for m in S_dd for n in S_dd]
A = [(k,l,m) for k in K_dd for l in K_dd for m in S_dd]
C1 = [(k,i) for k in K_dd for i in S_dd]
D = [(k,l) for k in K_dd for l in K_dd]

from gurobipy import Model, GRB, quicksum
mdl =Model('TRB_Santiago_16station')
mdl.Params.LazyConstraints = 1

#Decision Variables
tau = mdl.addVars(K_dd, vtype=GRB.BINARY)
mdl.addConstrs(tau[k] >= tau[k+1] for k in K_dd if k+1 in K_dd)
z = mdl.addVars(B, vtype=GRB.BINARY)
y = mdl.addVars(A, vtype=GRB.BINARY)
x = mdl.addVars(C1, vtype=GRB.BINARY)
a1 = mdl.addVars(C1, vtype=GRB.CONTINUOUS)
d = mdl.addVars(C1, vtype=GRB.CONTINUOUS)
h = mdl.addVars(D, vtype=GRB.CONTINUOUS)
alpha_up1 = mdl.addVars(K_u, vtype=GRB.BINARY)
alpha_dn2 = mdl.addVars(K_d, vtype=GRB.BINARY)
alpha_up3 = mdl.addVars(K_u, vtype=GRB.BINARY)
alpha_dn4 = mdl.addVars(K_d, vtype=GRB.BINARY)
beta_up2 = mdl.addVars(K_u, vtype=GRB.BINARY)
beta_dn1 = mdl.addVars(K_d, vtype=GRB.BINARY)
beta_up4 = mdl.addVars(K_u, vtype=GRB.BINARY)
beta_dn3 = mdl.addVars(K_d, vtype=GRB.BINARY)
RS1 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS1')
RS2 = mdl.addVar(vtype=GRB.INTEGER,lb=0, name='RS2')
RS3 = mdl.addVar(vtype=GRB.INTEGER, lb=0, name='RS3')
RS4 = mdl.addVar(vtype=GRB.INTEGER, lb=0,  name='RS4')
y_level = mdl.addVar(lb=0, vtype=GRB.INTEGER, name="y_level")
mdl.addConstr(quicksum(y[k,l,6]  for k in K_u for l in K_d) + quicksum(y[k,l,8]  for k in K_u for l in K_d) + quicksum(y[k,l,14] for k in K_d for l in K_u) + quicksum(y[k,l,16] for k in K_d for l in K_u) >= y_level, name="Pareto_y_level")
zz = mdl.addVars(B, vtype=GRB.CONTINUOUS)

#Master Problem: Operation Planning Constraints
mdl.addConstrs(z[k,1,6] + z[k,1,8] + z[k,3,6]+ z[k,3,8]  == tau[k] for k in K_u)
mdl.addConstrs(z[k,9,14] + z[k,9,16] + z[k,11,14] + z[k,11,16]  == tau[k] for k in K_d)
mdl.addConstrs(x[k,i] <= tau[k] for i in S_u for k in K_u)
mdl.addConstrs(x[k,i] <= tau[k] for i in S_d for k in K_d)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M*(1-z[k,1,6])for k in K_u)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M*(1-z[k,3,6])for k in K_u)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M*(1-z[k,3,6])for k in K_u)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M*(1-z[k,3,8])for k in K_u)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M*(1-z[k,9,14])for k in K_d)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M*(1-z[k,11,14])for k in K_d)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M*(1-z[k,11,14])for k in K_d)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M*(1-z[k,11,16])for k in K_d)
mdl.addConstrs(x[k-1,i] + x[k,i] >= tau[k] for k in K_u for i in S_u if k != 1)
mdl.addConstrs(x[k-1,i] + x[k,i] >= tau[k] for k in K_d for i in S_d if k != 10)
mdl.addConstr(d[1,1] == 27000)
mdl.addConstr(d[10,9] == 27000)
mdl.addConstrs(d[k,1] <= 28800 for k in K_u)
mdl.addConstrs(d[k,3] <= 28800 for k in K_u)
mdl.addConstrs(d[k,9] <= 28800 for k in K_d)
mdl.addConstrs(d[k,11] <= 28800 for k in K_d)
mdl.addConstrs(a1[k,i] - d[k,i-1] == pr[i-1,i] + x[k,i-1]*ac[i-1,i] +x[k,i]*dc[i-1,i] for k in K_u for i in S_u if i != 1)
mdl.addConstrs(a1[k,i] - d[k,i-1] == pr[i-1,i] + x[k,i-1]*ac[i-1,i] +x[k,i]*dc[i-1,i] for k in K_d for i in S_d if i != 9)
mdl.addConstrs(d[k,i] - a1[k,i] >= e[i]*x[k,i] for k in K_u for i in S_u)
mdl.addConstrs(d[k,i] - a1[k,i] >= e[i]*x[k,i] for k in K_d for i in S_d)
mdl.addConstrs(quicksum(x[k,i]for i in S_u) >= 4 for k in K_u) 
mdl.addConstrs(quicksum(x[k,i]for i in S_d) >= 4 for k in K_d)
mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_u if k!= 1)
mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_u if k!= 1)
mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_d if k!= 10)
mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_d if k!= 10)
mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_u for i in S_u if k != 1)
mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_d for i in S_d if k != 10)
mdl.addConstrs(2*y[k,l,6] <= z[k,1,6] +z[l,11,16] + z[k,3,6] + z[l,11,14] for k in K_u for l in K_d)
mdl.addConstrs(2*y[k,l,8] <= z[k,1,8] + z[l,9,16] +z[k,3,8] + z[l,9,14] for k in K_u for l in K_d)
mdl.addConstrs(2*y[k,l,14] <= z[k,9,14] + z[l,3,8]+z[k,11,14] + z[l,3,6] for k in K_d for l in K_u)
mdl.addConstrs(2*y[k,l,16] <= z[k,9,16] + z[l,1,8]+z[k,11,16] + z[l,1,6] for k in K_d for l in K_u)
mdl.addConstrs(a1[l,9] - d[k,8] >= delta_min + M_time*(y[k,l,8]-1) for k in K_u for l in K_d)
mdl.addConstrs(a1[l,11] - d[k,6] >= delta_min + M_time*(y[k,l,6]-1) for k in K_u for l in K_d)
mdl.addConstrs(a1[l,1] - d[k,16] >= delta_min + M_time*(y[k,l,16]-1) for k in K_d for l in K_u)
mdl.addConstrs(a1[l,3] - d[k,14] >= delta_min + M_time*(y[k,l,14]-1) for k in K_d for l in K_u)
mdl.addConstrs(quicksum(y[l,k,16]for l in K_d)+quicksum(y[l,k,14]for l in K_d) + alpha_up1[k] + alpha_up3[k] == tau[k] for k in K_u)
mdl.addConstrs(quicksum(y[l,k,6]for l in K_u)+quicksum(y[l,k,8]for l in K_u) + alpha_dn2[k] + alpha_dn4[k] == tau[k] for k in K_d)
mdl.addConstrs(quicksum(y[k,l,6]for l in K_d)+quicksum(y[k,l,8]for l in K_d) + beta_up2[k]+ beta_up4[k] == tau[k] for k in K_u)
mdl.addConstrs(quicksum(y[k,l,16]for l in K_u)+quicksum(y[k,l,14]for l in K_u) + beta_dn1[k]+ beta_dn3[k] == tau[k] for k in K_d)
mdl.addConstr(quicksum(alpha_up1[k] for k in K_u) + quicksum(alpha_dn2[k] for k in K_d)+ quicksum(alpha_up3[k] for k in K_u) + quicksum(alpha_dn4[k] for k in K_d)<= RS)
mdl.addConstr(quicksum(alpha_up1[k] for k in K_u) <= RS1)
mdl.addConstr(quicksum(alpha_dn2[k] for k in K_d) <= RS2)
mdl.addConstr(quicksum(alpha_up3[k] for k in K_u) <= RS3)
mdl.addConstr(quicksum(alpha_dn4[k] for k in K_d) <= RS4)
mdl.addConstrs(alpha_up1[k] <= z[k,1,8] + z[k,1,6] for k in K_u)
mdl.addConstrs(alpha_dn2[k] <= z[k,9,16] + z[k,9,14] for k in K_d)
mdl.addConstrs(beta_dn1[k] <= z[k,9,16]+ z[k,11,16]for k in K_d)
mdl.addConstrs(beta_up2[k] <= z[k,1,8] + z[k,3,8] for k in K_u)
mdl.addConstrs(alpha_up3[k] <= z[k,3,8] + z[k,3,6] for k in K_u)
mdl.addConstrs(alpha_dn4[k] <= z[k,11,16] + z[k,11,14] for k in K_d)
mdl.addConstrs(beta_dn3[k] <= z[k,9,14]+z[k,11,14] for k in K_d)
mdl.addConstrs(beta_up4[k] <= z[k,1,6] + z[k,3,6] for k in K_u)

#Constraints used to linearize objective function
mdl.addConstrs(zz[k,1,6] >= (d[k,6] - d[k,1]) - M*(1-z[k,1,6])for k in K_u)
mdl.addConstrs(zz[k,1,8] >= (d[k,8] - d[k,1]) - M*(1-z[k,1,8])for k in K_u)
mdl.addConstrs(zz[k,3,6] >= (d[k,6] - d[k,3]) - M*(1-z[k,3,6])for k in K_u)
mdl.addConstrs(zz[k,3,8] >= (d[k,8] - d[k,3]) - M*(1-z[k,3,8])for k in K_u)
mdl.addConstrs(zz[k,9,14] >= (d[k,14] - d[k,9]) - M*(1-z[k,9,14])for k in K_d)
mdl.addConstrs(zz[k,9,16] >= (d[k,16] - d[k,9]) - M*(1-z[k,9,16])for k in K_d)
mdl.addConstrs(zz[k,11,14] >= (d[k,14] - d[k,11]) - M*(1-z[k,11,14])for k in K_d)
mdl.addConstrs(zz[k,11,16] >= (d[k,16] - d[k,11]) - M*(1-z[k,11,16])for k in K_d)
mdl.addConstrs(zz[k,1,6] <= M*(z[k,1,6])for k in K_u)
mdl.addConstrs(zz[k,1,8] <= M*(z[k,1,8])for k in K_u)
mdl.addConstrs(zz[k,3,6] <= M*(z[k,3,6])for k in K_u)
mdl.addConstrs(zz[k,3,8] <= M*(z[k,3,8])for k in K_u)
mdl.addConstrs(zz[k,9,14] <= M*(z[k,9,14])for k in K_d)
mdl.addConstrs(zz[k,9,16] <= M*(z[k,9,16])for k in K_d)
mdl.addConstrs(zz[k,11,14] <= M*(z[k,11,14])for k in K_d)
mdl.addConstrs(zz[k,11,16] <= M*(z[k,11,16])for k in K_d)
               
mdl.setObjective(quicksum(zz[k,1,6]for k in K_u) +quicksum(zz[k,1,8]for k in K_u) +quicksum(zz[k,3,6]for k in K_u) +quicksum(zz[k,3,8]for k in K_u)
                 +quicksum(zz[k,9,14]for k in K_d) +quicksum(zz[k,9,16]for k in K_d) +quicksum(zz[k,11,14]for k in K_d) +quicksum(zz[k,11,16]for k in K_d)
                 +quicksum(d[k,i]-d[k-1,i] for k in K_u for i in S_u if k!=1) +quicksum(d[k,i]-d[k-1,i] for k in K_d for i in S_d if k!=10)
)
#Sub-Problem
def solve_passenger_subproblem(x_val, d_val, z_val):
    sub = Model("Passenger_Subproblem")
    w   = sub.addVars(B, lb=0)
    w_b = sub.addVars(C1, lb=0)
    w_b1= sub.addVars(B, lb=0)
    n_b = sub.addVars(C1, lb=0)
    n_b1= sub.addVars(B, lb=0)
    n_a = sub.addVars(C1, lb=0)
    n1  = sub.addVars(C1, lb=0, ub=C)
    v   = sub.addVars(B, lb=0)
    sai = sub.addVars(C1, vtype=GRB.BINARY)
    
    #Sub-Problem: Passenger Demand Constraints (Time Invarient)
    sub.addConstrs(w[1,i,j] == (p[i,j]/1800)*120 for i in S_u for j in S_u if i < j)
    sub.addConstrs(w[k,i,j] == v[k-1,i,j] + (p[i,j]/1800)*(d_val[(k,i)] - d_val[(k-1,i)]) for k in K_u for i in S_u for j in S_u if i < j if k != 1)
    sub.addConstrs(w[10,i,j] == (p[i,j]/1800)*120 for i in S_d for j in S_d if i < j)
    sub.addConstrs(w[k,i,j] == v[k-1,i,j] + (p[i,j]/1800)*(d_val[(k,i)] - d_val[(k-1,i)]) for k in K_d for i in S_d for j in S_d if i < j if k != 10)
    sub.addConstrs(w_b[k,i] == quicksum(w_b1[k,i,j] for j in S_u if i < j) for k in K_u for i in S_u)
    sub.addConstrs(w_b[k,i] == quicksum(w_b1[k,i,j] for j in S_d if i < j) for k in K_d for i in S_d)
    sub.addConstrs(n_a[k,i] == quicksum(n_b1[k,j,i] for j in S_u if j < i) for k in K_u for i in S_u)
    sub.addConstrs(n_a[k,i] == quicksum(n_b1[k,j,i] for j in S_d if j < i) for k in K_d for i in S_d)
    sub.addConstrs(n_b[k,i] - n_b1[k,i,j] == w_b[k,i] - w_b1[k,i,j] for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(n_b[k,i] - n_b1[k,i,j] == w_b[k,i] - w_b1[k,i,j] for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(w_b1[k,i,j] >= 0 for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(w_b1[k,i,j] <= w[k,i,j] for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(w_b1[k,i,j] <= M*x_val[(k,i)] for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(w_b1[k,i,j] <= M*x_val[(k,j)] for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M*(2-x_val[(k,i)]-x_val[(k,j)]) for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(w_b1[k,i,j] >= 0 for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(w_b1[k,i,j] <= w[k,i,j] for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(w_b1[k,i,j] <= M*x_val[(k,i)] for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(w_b1[k,i,j] <= M*x_val[(k,j)] for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M*(2-x_val[(k,i)]-x_val[(k,j)]) for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_u for i in S_u)
    sub.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_u for i in S_u if i!= 1)
    sub.addConstrs(n_b[k,i] >= w_b[k,i] - M*(1-sai[k,i]) for k in K_u for i in S_u)
    sub.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M*sai[k,i] for k in K_u for i in S_u if i != 1)
    sub.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_d for i in S_d)
    sub.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_d for i in S_d if i!= 9)
    sub.addConstrs(n_b[k,i] >= w_b[k,i] - M*(1-sai[k,i]) for k in K_d for i in S_d)
    sub.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M*sai[k,i] for k in K_d for i in S_d if i != 9)
    sub.addConstrs(n1[k,6] <= C * (1 - z_val[k,1,6]) for k in K_u)
    sub.addConstrs(n1[k,6] <= C * (1 - z_val[k,3,6]) for k in K_u)
    sub.addConstrs(n1[k,8] <= C * (1 - z_val[k,1,8]) for k in K_u)
    sub.addConstrs(n1[k,8] <= C * (1 - z_val[k,3,8]) for k in K_u)
    sub.addConstrs(n1[k,14] <= C * (1 - z_val[k,9,14]) for k in K_d)
    sub.addConstrs(n1[k,14] <= C * (1 - z_val[k,11,14]) for k in K_d)
    sub.addConstrs(n1[k,16] <= C * (1 - z_val[k,9,16]) for k in K_d)
    sub.addConstrs(n1[k,16] <= C * (1 - z_val[k,11,16]) for k in K_d)
    sub.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_u for i in S_u if i != 1)
    sub.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_u for i in S_u if i != 1)
    sub.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_d for i in S_d if i != 9)
    sub.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M for k in K_d for i in S_d if i != 9)
    sub.addConstrs(n1[k,1] == n_b[k,1] for k in K_u)
    sub.addConstrs(n1[k,9] == n_b[k,9] for k in K_d)
    sub.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_u for i in S_u if i != 1)
    sub.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M*(2 - x_val[(k,i)] - x_val[(k,j)]) for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M*(2 - x_val[(k,i)] - x_val[(k,j)]) for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(v[k,i,j] <= M*x_val[(k,i)] for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(v[k,i,j] <= M*x_val[(k,j)] for k in K_u for i in S_u for j in S_u if i < j)
    sub.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_d for i in S_d if i != 9)
    sub.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M*(2 - x_val[(k,i)] - x_val[(k,j)]) for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M*(2 - x_val[(k,i)] - x_val[(k,j)]) for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(v[k,i,j] <= M*x_val[(k,i)] for k in K_d for i in S_d for j in S_d if i < j)
    sub.addConstrs(v[k,i,j] <= M*x_val[(k,j)] for k in K_d for i in S_d for j in S_d if i < j)
    sub.setObjective(0)
    sub.optimize()
    if sub.Status != GRB.OPTIMAL:
        return None
    return {
        "w":   {(k,i,j): w[k,i,j].X   for (k,i,j) in B},
        "w_b": {(k,i):   w_b[k,i].X   for (k,i)   in C1},
        "n1":  {(k,i):   n1[k,i].X    for (k,i)   in C1},
        "v":   {(k,i,j): v[k,i,j].X   for (k,i,j) in B},
        "n_b1":{(k,i,j): n_b1[k,i,j].X for (k,i,j) in B},
    }

#Benders Callback
feasibility_cache = {}
def benders_callback(model, where):
    if where == GRB.Callback.MIPSOL:
        x_val = {(k,i): model.cbGetSolution(x[k,i]) for (k,i) in C1}
        d_val = {(k,i): model.cbGetSolution(d[k,i]) for (k,i) in C1}
        z_val = {(k,m,n): model.cbGetSolution(z[k,m,n]) for (k,m,n) in B}
        pattern_key = tuple(
            (k,i) for (k,i) in C1 if x_val[(k,i)] > 0.5
        )
        if pattern_key in feasibility_cache:
            if not feasibility_cache[pattern_key]:
                violated_keys = list(pattern_key)
                model.cbLazy(
                    quicksum(x[k,i] for (k,i) in violated_keys)
                    <= len(violated_keys) - 1
                )
            return
        sub_vars = solve_passenger_subproblem(x_val, d_val, z_val)
        feasibility_cache[pattern_key] = (sub_vars is not None)
        if sub_vars is not None:
            return
        violated_keys = [(k,i) for (k,i) in C1 if x_val[(k,i)] > 0.5]
        model.cbLazy(
            quicksum(x[k,i] for (k,i) in violated_keys)
            <= len(violated_keys) - 1
        )
def analyze_passenger_solution(
    solution_id,
    tau_val,
    x_val,
    d_val,
    z_val,
    h_val,
    sub_vars,
    output_prefix="Final"
):
    w     = sub_vars["w"]
    w_b   = sub_vars["w_b"]
    n1    = sub_vars["n1"]
    v     = sub_vars["v"]
    n_b1  = sub_vars["n_b1"]

    eps = 1e-6
    cap = C
    #Performance Metrics
    plf = max( (n1[(k, i)] / cap for (k, i) in n1 if tau_val.get(k, 0) > 0.5), default=0.0 )
    overloaded = [(k, i, n1[(k, i)] / cap) for (k, i) in n1 if tau_val.get(k, 0) > 0.5 and n1[(k, i)] / cap > 0.85]
    cong_seconds = sum(h_val.get((k - 1, k), h_min) for k, _, _ in overloaded)
    sum_v = sum(v[key] for key in v if x_val.get((key[0], key[1]), 0) > 0.5 and x_val.get((key[0], key[2]), 0) > 0.5)
    sum_w = sum(w[key] for key in w if x_val.get((key[0], key[1]), 0) > 0.5 and x_val.get((key[0], key[2]), 0) > 0.5)

    left_behind = sum_v / sum_w if sum_w > eps else 0.0
    dcr = {}
    for (k, i), arr in w_b.items():
        if tau_val.get(k, 0) < 0.5:
            continue
        hw = h_val.get((k - 1, k), h_min)
        dcr[(k, i)] = arr * hw / cap
    peak_dcr = max(dcr.values(), default=0.0)
MAX_Y = 12  # safe upper bound
y_level.lb = 0
#ε-Constraint controller 
mdl.optimize(benders_callback)
mdl.update()
y_vals = mdl.getAttr("X", y)
y_star = int(round(sum(y_vals.values())))
level = y_star

while level <= MAX_Y:
    y_level.lb = level
    mdl.update()
    mdl.setParam(GRB.Param.TimeLimit, 14400)
    mdl.optimize(benders_callback)
    if mdl.Status == GRB.TIME_LIMIT:
        if mdl.SolCount > 0:
            print(f"[TIME LIMIT] Extracting incumbent for ε = {level} (with gap)")
            mdl.update()
            x_val   = dict(mdl.getAttr("X", x))
            d_val   = dict(mdl.getAttr("X", d))
            z_val   = dict(mdl.getAttr("X", z))
            h_val   = dict(mdl.getAttr("X", h))
            tau_val = dict(mdl.getAttr("X", tau))
            sub_vars = solve_passenger_subproblem(x_val, d_val, z_val)
            if sub_vars is not None:
                analyze_passenger_solution(
                    solution_id=level,
                    tau_val=tau_val,
                    x_val=x_val,
                    d_val=d_val,
                    z_val=z_val,
                    h_val=h_val,
                    sub_vars=sub_vars,
                    output_prefix="TIME_LIMIT")
        else:
            print(f"[TIME LIMIT] No incumbent available for ε = {level}")
        level += 1
        continue
    if mdl.Status == GRB.INFEASIBLE or mdl.SolCount == 0:
        print(f"No feasible solution for ε = {level}")
    print(f"Hybrid ε+Benders solution found for y_level = {level}")
    #Extract master solution
    if mdl.SolCount == 0:
        print(f"[Skip] No incumbent available at ε = {level}")
        level += 1
        continue
    mdl.update()
    x_val   = dict(mdl.getAttr("X", x))
    d_val   = dict(mdl.getAttr("X", d))
    z_val   = dict(mdl.getAttr("X", z))
    h_val   = dict(mdl.getAttr("X", h))
    tau_val = dict(mdl.getAttr("X", tau))
    #Solve passenger subproblem once
    sub_vars = solve_passenger_subproblem(x_val, d_val, z_val)
    if sub_vars is not None:
        analyze_passenger_solution(
            solution_id=level,
            tau_val=tau_val,
            x_val=x_val,
            d_val=d_val,
            z_val=z_val,
            h_val=h_val,
            sub_vars=sub_vars)
    level += 1