import pandas as pd
# Model Parameters
S_up = 8
S_dn = 16
K_up = 6
K_dn = 12
h_max = 360
h_min = 90
delta_min = 135
RS = 5
C = 250
M_time = 100000.0
M_flow = 5000.0

#Reading data from excel
df1= pd.read_excel('Input_Parameters.xlsx', 'Running_time')
my_tuple = [tuple(x) for x in df1.values]
r =dict((tuple((a, b)), c) for a,b,c in df1.values)
df1 = pd.read_excel('Input_Parameters.xlsx', 'Dwelling_time')
my_tuple = [tuple(x) for x in df1.values]
e = dict((a, c) for a, c in my_tuple)
df1= pd.read_excel('Input_Parameters.xlsx', 'Demand_30(7.30-8.00am)mint')
my_tuple = [tuple(x) for x in df1.values]
p =dict((tuple((a, b)), c) for a,b,c in df1.values)
OD_p = dict(p)

#Index Sets
S_d = [i for i in range(9,S_dn+1)]
S_dd = [i for i in range(1,S_dn+1)]
K_d = [i for i in range(7,K_dn+1)]
K_dd = [i for i in range(1,K_dn+1)]
K_u = [i for i in range(1,K_up+1)]
S_u = [i for i in range(1,S_up+1)]
B = [(k,m,n) for k in K_dd for m in S_dd for n in S_dd]
A = [(k,l,m) for k in K_dd for l in K_dd for m in S_dd]
C1 = [(k,i) for k in K_dd for i in S_dd]
D = [(k,l) for k in K_dd for l in K_dd]

from gurobipy import Model, GRB, quicksum
mdl =Model('Model_1_Minimizing_Opearation_Cost')

#Decision Variables
tau = mdl.addVars(K_dd, vtype=GRB.BINARY)
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
w = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='w')
w_b = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='w_b')
w_b1 = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='w_b1')
n_b = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='n_b')
n_b1 = mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='n_b1')
n_a = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0.0, name='n_a')
n1 = mdl.addVars(C1, vtype=GRB.CONTINUOUS, lb=0, ub=C, name="n1")
v =  mdl.addVars(B, vtype=GRB.CONTINUOUS, lb=0.0, name='v')
sai = mdl.addVars(C1, vtype=GRB.BINARY)

#Operation Zone Selection Constraints
mdl.addConstrs(z[k,1,6] + z[k,1,8] + z[k,3,6]+ z[k,3,8]  == tau[k] for k in K_u)
mdl.addConstrs(z[k,9,14] + z[k,9,16] + z[k,11,14] + z[k,11,16]  == tau[k] for k in K_d)

#Served Station Constraints
mdl.addConstrs(x[k,i] <= tau[k] for i in S_u for k in K_u)
mdl.addConstrs(x[k,i] <= tau[k] for i in S_d for k in K_d)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M_flow*(1-z[k,1,6])for k in K_u)
mdl.addConstrs(x[k,i] >= z[k,1,6] for k in K_u for i in S_u if i >= 1 if i <= 6)
mdl.addConstrs(x[k,i] >= z[k,1,8] for k in K_u for i in S_u if i >= 1 if i <= 8)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M_flow*(1-z[k,3,6])for k in K_u)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i>6) <= M_flow*(1-z[k,3,6])for k in K_u)
mdl.addConstrs(x[k,i] >= z[k,3,6] for k in K_u for i in S_u if i >= 3 if i <= 6)
mdl.addConstrs(quicksum(x[k,i]for i in S_u if i<3) <= M_flow*(1-z[k,3,8])for k in K_u)
mdl.addConstrs(x[k,i] >= z[k,3,8] for k in K_u for i in S_u if i >= 3 if i <= 8)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M_flow*(1-z[k,9,14])for k in K_d)
mdl.addConstrs(x[k,i] >= z[k,9,14] for k in K_d for i in S_d if i >= 9 if i <= 14)
mdl.addConstrs(x[k,i] >= z[k,9,16] for k in K_d for i in S_d if i >= 9 if i <= 16)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M_flow*(1-z[k,11,14])for k in K_d)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i>14) <= M_flow*(1-z[k,11,14])for k in K_d)
mdl.addConstrs(x[k,i] >= z[k,11,14] for k in K_d for i in S_d if i >= 11 if i <= 14)
mdl.addConstrs(quicksum(x[k,i]for i in S_d if i<11) <= M_flow*(1-z[k,11,16])for k in K_d)
mdl.addConstrs(x[k,i] >= z[k,11,16] for k in K_d for i in S_d if i >= 11 if i <= 16)
mdl.addConstrs(x[k-1,i] + x[k,i] >= 1 for k in K_u for i in S_u if k != 1)
mdl.addConstrs(x[k-1,i] + x[k,i] >= 1 for k in K_d for i in S_d if k != 7)

#Arrival and Departure Time Constraints
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

#Headway Constraints
mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_u if k!= 1)
mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_u if k!= 1)
mdl.addConstrs(h_min*tau[k] <= h[k-1,k] for k in K_d if k!= 7)
mdl.addConstrs(h[k-1,k] <= h_max*tau[k] for k in K_d if k!= 7)
mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_u for i in S_u if k != 1)
mdl.addConstrs(d[k,i] == d[k-1,i] + h[k-1,k] for k in K_d for i in S_d if k != 7)

#Turnaround Constraints
mdl.addConstrs(2*y[k,l,6] <= z[k,1,6] +z[l,11,16] + z[k,3,6] + z[l,11,14] for k in K_u for l in K_d)
mdl.addConstrs(2*y[k,l,8] <= z[k,1,8] + z[l,9,16] +z[k,3,8] + z[l,9,14] for k in K_u for l in K_d)
mdl.addConstrs(2*y[k,l,14] <= z[k,9,14] + z[l,3,8]+z[k,11,14] + z[l,3,6] for k in K_d for l in K_u)
mdl.addConstrs(2*y[k,l,16] <= z[k,9,16] + z[l,1,8]+z[k,11,16] + z[l,1,6] for k in K_d for l in K_u)
mdl.addConstrs(a1[l,9] - d[k,8] >= delta_min + M_time*(y[k,l,8]-1) for k in K_u for l in K_d)
mdl.addConstrs(a1[l,11] - d[k,6] >= delta_min + M_time*(y[k,l,6]-1) for k in K_u for l in K_d)
mdl.addConstrs(a1[l,1] - d[k,16] >= delta_min + M_time*(y[k,l,16]-1) for k in K_d for l in K_u)
mdl.addConstrs(a1[l,3] - d[k,14] >= delta_min + M_time*(y[k,l,14]-1) for k in K_d for l in K_u)

#Rolling Stock Assignment Constraints
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
mdl.addConstrs(w_b1[k,i,j] <= M_flow*x[k,i] for k in K_u for i in S_u for j in S_u if i < j)
mdl.addConstrs(w_b1[k,i,j] <= M_flow*x[k,j] for k in K_u for i in S_u for j in S_u if i < j)
mdl.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M_flow*(2-x[k,i]-x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
mdl.addConstrs(w_b1[k,i,j] >= 0 for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(w_b1[k,i,j] <= w[k,i,j] for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(w_b1[k,i,j] <= M_flow*x[k,i] for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(w_b1[k,i,j] <= M_flow*x[k,j] for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(w_b1[k,i,j] >= w[k,i,j] - M_flow*(2-x[k,i]-x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_u for i in S_u)
mdl.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_u for i in S_u if i!= 1)
mdl.addConstrs(n_b[k,i] >= w_b[k,i] - M_flow*(1-sai[k,i]) for k in K_u for i in S_u)
mdl.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M_flow*sai[k,i] for k in K_u for i in S_u if i != 1)
mdl.addConstrs(n_b[k,i] <= w_b[k,i] for k in K_d for i in S_d)
mdl.addConstrs(n_b[k,i] <= C- n1[k,i-1] + n_a[k,i] for k in K_d for i in S_d if i!= 9)
mdl.addConstrs(n_b[k,i] >= w_b[k,i] - M_flow*(1-sai[k,i]) for k in K_d for i in S_d)
mdl.addConstrs(n_b[k,i] >= C - n1[k,i-1] + n_a[k,i]- M_flow*sai[k,i] for k in K_d for i in S_d if i != 9)
mdl.addConstrs(n1[k,6] <= C * (1 - z[k,1,6]) for k in K_u)
mdl.addConstrs(n1[k,6] <= C * (1 - z[k,3,6]) for k in K_u)
mdl.addConstrs(n1[k,8] <= C * (1 - z[k,1,8]) for k in K_u)
mdl.addConstrs(n1[k,8] <= C * (1 - z[k,3,8]) for k in K_u)
mdl.addConstrs(n1[k,14] <= C * (1 - z[k,9,14]) for k in K_d)
mdl.addConstrs(n1[k,14] <= C * (1 - z[k,11,14]) for k in K_d)
mdl.addConstrs(n1[k,16] <= C * (1 - z[k,9,16]) for k in K_d)
mdl.addConstrs(n1[k,16] <= C * (1 - z[k,11,16]) for k in K_d)
mdl.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M_flow for k in K_u for i in S_u if i != 1)
mdl.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M_flow for k in K_u for i in S_u if i != 1)
mdl.addConstrs(sai[k,i] >= ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M_flow for k in K_d for i in S_d if i != 9)
mdl.addConstrs(sai[k,i] <= 1 + ((C - n1[k,i-1] + n_a[k,i])- w_b[k,i])/M_flow for k in K_d for i in S_d if i != 9)
mdl.addConstrs(n1[k,1] == n_b[k,1] for k in K_u)
mdl.addConstrs(n1[k,9] == n_b[k,9] for k in K_d)
mdl.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_u for i in S_u if i != 1)
mdl.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M_flow*(2 - x[k,i] - x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
mdl.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M_flow*(2 - x[k,i] - x[k,j]) for k in K_u for i in S_u for j in S_u if i < j)
mdl.addConstrs(v[k,i,j] <= M_flow*x[k,i] for k in K_u for i in S_u for j in S_u if i < j)
mdl.addConstrs(v[k,i,j] <= M_flow*x[k,j] for k in K_u for i in S_u for j in S_u if i < j)
mdl.addConstrs(n1[k,i] == n1[k,i-1] - n_a[k,i] + n_b[k,i] for k in K_d for i in S_d if i != 9)
mdl.addConstrs(v[k,i,j] <= w[k,i,j] - n_b1[k,i,j] + M_flow*(2 - x[k,i] - x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(v[k,i,j] >= w[k,i,j] - n_b1[k,i,j] - M_flow*(2 - x[k,i] - x[k,j]) for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(v[k,i,j] <= M_flow*x[k,i] for k in K_d for i in S_d for j in S_d if i < j)
mdl.addConstrs(v[k,i,j] <= M_flow*x[k,j] for k in K_d for i in S_d for j in S_d if i < j)

#Objective Function 
mdl.setObjective(quicksum(y[l,k,6] for l in K_u for k in K_d) + quicksum(y[l,k,8] for l in K_u for k in K_d) +
                  quicksum(y[l,k,14] for l in K_d for k in K_u) + quicksum(y[l,k,16] for l in K_d for k in K_u))
mdl.modelSense = GRB.MAXIMIZE
mdl.Params.TimeLimit = 14400
mdl.optimize()

runtime = mdl.Runtime
print("The run time is %f" % runtime)
print("\nObjective Function Value:")
print("ObjVal: {}".format(mdl.ObjVal))
for k in K_u:
     if tau[k].x > 0:
         print("tau[{}]: {}".format(k, tau[k].x))
for k in K_d:
     if tau[k].x > 0:
         print("tau[{}]: {}".format(k, tau[k].x))
if RS1.x > 0:
     print("RS1: {}".format(RS1.x))
if RS2.x > 0:
     print("RS2: {}".format(RS2.x))
if RS3.x > 0:
     print("RS3: {}".format(RS3.x))
if RS4.x > 0:
     print("RS4: {}".format(RS4.x))
for d_val in D:
     if h[d_val].x > 0:
         print("h[{}]: {}".format(d_val, h[d_val].x))
for b in B:
     if z[b].x > 0:
         print("z[{}]: {}".format(b, z[b].x))
for a in A:
     if y[a].x > 0:
         print("y[{}]: {}".format(a, y[a].x))
for c in C1:
     if x[c].x > 0:
         print("x[{}]: {}".format(c, x[c].x))
for c in C1:
     if a1[c].x > 0:
         print("a1[{}]: {}".format(c, a1[c].x))
for c in C1:
     if d[c].x > 0:
         print("d[{}]: {}".format(c, d[c].x))
for k in K_u:
     if alpha_up1[k].x > 0:
         print("alpha_u1[{}]: {}".format(k, alpha_up1[k].x))
for k in K_d:
     if alpha_dn2[k].x > 0:
         print("alpha_d2[{}]: {}".format(k, alpha_dn2[k].x))
for k in K_u:
     if beta_up2[k].x > 0:
         print("beta_u2[{}]: {}".format(k, beta_up2[k].x))
for k in K_d:
     if beta_dn1[k].x > 0:
         print("beta_d1[{}]: {}".format(k, beta_dn1[k].x))
for k in K_u:
     if alpha_up3[k].x > 0:
         print("alpha_u3[{}]: {}".format(k, alpha_up3[k].x))
for k in K_d:
     if alpha_dn4[k].x > 0:
         print("alpha_d4[{}]: {}".format(k, alpha_dn4[k].x))
for k in K_u:
     if beta_up4[k].x > 0:
         print("beta_u4[{}]: {}".format(k, beta_up4[k].x))
for k in K_d:
     if beta_dn3[k].x > 0:
         print("beta_d3[{}]: {}".format(k, beta_dn3[k].x))
for b in B:
     if w[b].x > 0:
         print("w[{}]: {}".format(b, w[b].x))
for c in C1:
     if w_b[c].x > 0:
         print("w_b[{}]: {}".format(c, w_b[c].x))
for b in B:
     if w_b1[b].x > 0:
         print("w_b1[{}]: {}".format(b, w_b1[b].x))
for c in C1:
     if n_b[c].x > 0:
         print("n_b[{}]: {}".format(c, n_b[c].x))
for b in B:
     if n_b1[b].x > 0:
         print("n_b1[{}]: {}".format(b, n_b1[b].x))
for c in C1:
     if n_a[c].x > 0:
         print("n_a[{}]: {}".format(c, n_a[c].x))
for c in C1:
     if n1[c].x > 0:
         print("n1[{}]: {}".format(c, n1[c].x))
for b in B:
     if v[b].x > 0:
         print("v[{}]: {}".format(b, v[b].x))

# Performance metrics: PLF, Congestion Duration, Left-Behind Rate, D-to-C
eps = 1e-6
cap = C
# 1) Peak Load Factor (PLF) = max_{k,i} n1[k,i] / C
plf_by_ki = {}
for (k, i) in n1.keys():
    try:
        if (k in tau.keys()) and (tau[k].x > 0.5):
            val = n1[(k, i)].x
            plf = val / float(cap) if cap > eps else float('inf')
            plf_by_ki[(k, i)] = plf
    except Exception:
        continue
max_plf = max(plf_by_ki.values()) if plf_by_ki else 0.0

# list overloaded (PLF>0.85)
overloaded = [(k, i, p) for (k, i), p in plf_by_ki.items() if p > 0.85]

# 2) Congestion Duration: sum of headways (seconds) for services that have PLF>0.85
cong_seconds = 0.0
services_with_cong = set([k for (k, i, p) in overloaded])
for k in services_with_cong:
    idx = (k - 1, k)
    if idx in h.keys():
        cong_seconds += h[idx].x
    else:
        cong_seconds += 120.0
cong_minutes = cong_seconds / 60.0

# 3) Left-Behind Passenger Rate (only count ODs where train serves both origin and destination)
sum_v = 0.0
sum_w = 0.0
for (k, i, j) in B:
    try:
        if ((k, i) in x.keys() and (k, j) in x.keys() and x[(k, i)].x > 0.5 and x[(k, j)].x > 0.5):
            sum_v += float(v[(k, i, j)].x)
            sum_w += float(w[(k, i, j)].x)
    except Exception:
        continue

left_behind_rate = (sum_v / sum_w) if sum_w > eps else 0.0

# 4) Demand-to-Capacity Ratio (DCR) = arrivals / (capacity/headway)
dcr_by_ki = {}
hw_by_k = {}   

for (k, i) in w_b.keys():
    try:
        if (k in tau.keys()) and (tau[k].x > 0.5):
            arrivals = float(w_b[(k, i)].x)
            idx = (k - 1, k)
            if idx in h.keys():
                hw = max(h[idx].x, eps)
            else:
                hw = h_min
            hw_by_k[k] = hw
            cap_per_time = float(cap) / hw if hw > eps else float('inf')
            dcr_val = arrivals / cap_per_time if cap_per_time > eps else float('inf')
            dcr_by_ki[(k, i)] = dcr_val

    except Exception:
        continue
peak_dcr = max(dcr_by_ki.values()) if dcr_by_ki else 0.0

print('\nPerformance metrics:')
print('Peak Load Factor (PLF) = max_k,i n1[k,i]/C = {:.4f}'.format(max_plf))
if overloaded:
    print('Overloaded (PLF>0.85) (k,i,PLF):')
    for (k, i, p) in overloaded:
        print('  {}: {:.4f}'.format((k, i), p))
else:
    print('No overloaded (PLF>0.85) stations/services.')
print('Congestion duration (seconds): {:.1f}  | minutes: {:.2f}'.format(cong_seconds, cong_minutes))
print('Left-Behind Passenger Rate (sum v / sum w): {:.4f}'.format(left_behind_rate))
print('Peak Demand-to-Capacity Ratio (DCR_peak): {:.4f}'.format(peak_dcr))