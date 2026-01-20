import pandas as pd
# Model Parameters
S_U = 8
S_D = 16
K_U = 6
K_D = 12
h_max = 360
h_min = 90
delta_min = 135
RS = 14
C = 250
M_time = 100000.0
M_flow = 5000.0
M = M_flow

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

#Index Sets
S_d = [i for i in range(9,S_D+1)]
S_dd = [i for i in range(1,S_D+1)]
K_d = [i for i in range(7,K_D+1)]
K_dd = [i for i in range(1,K_D+1)]
K_u = [i for i in range(1,K_U+1)]
S_u = [i for i in range(1,S_U+1)]
B = [(k,m,n) for k in K_dd for m in S_dd for n in S_dd]
A = [(k,l,m) for k in K_dd for l in K_dd for m in S_dd]
C1 = [(k,i) for k in K_dd for i in S_dd]
D = [(k,l) for k in K_dd for l in K_dd]

from gurobipy import Model, GRB, quicksum
mdl =Model('TRB_Santiago_16station')

#Decision Variables
tau = mdl.addVars(K_dd, vtype=GRB.BINARY)
z = mdl.addVars(B, vtype=GRB.BINARY)
y = mdl.addVars(A, vtype=GRB.BINARY)
x = mdl.addVars(C1, vtype=GRB.BINARY)
a1 = mdl.addVars(C1, vtype=GRB.CONTINUOUS)
d = mdl.addVars(C1, vtype=GRB.CONTINUOUS)
h = mdl.addVars(D, vtype=GRB.CONTINUOUS)
alpha_u1 = mdl.addVars(K_u, vtype=GRB.BINARY)
alpha_d2 = mdl.addVars(K_d, vtype=GRB.BINARY)
alpha_u3 = mdl.addVars(K_u, vtype=GRB.BINARY)
alpha_d4 = mdl.addVars(K_d, vtype=GRB.BINARY)
beta_u2 = mdl.addVars(K_u, vtype=GRB.BINARY)
beta_d1 = mdl.addVars(K_d, vtype=GRB.BINARY)
beta_u4 = mdl.addVars(K_u, vtype=GRB.BINARY)
beta_d3 = mdl.addVars(K_d, vtype=GRB.BINARY)
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
zz = mdl.addVars(B, vtype=GRB.CONTINUOUS)

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
               
mdl.setObjective(quicksum(zz[k,1,6]for k in K_u) +quicksum(zz[k,1,8]for k in K_u) +quicksum(zz[k,3,6]for k in K_u)
                 +quicksum(zz[k,3,8]for k in K_u) +quicksum(zz[k,9,14]for k in K_d)+quicksum(zz[k,9,16]for k in K_d)
                 +quicksum(zz[k,11,14]for k in K_d)+quicksum(zz[k,11,16]for k in K_d)
                 +quicksum(d[k,i]-d[k-1,i] for k in K_u for i in S_u if k!=1) +quicksum(d[k,i]-d[k-1,i] for k in K_d for i in S_d if k!=7))

mdl.modelSense = GRB.MINIMIZE
mdl.setParam('OutputFlag', 0)
mdl.optimize()
       
y_val = {} 
for a in A:
    y_val[a] = y[a].x
    if y[a].x > 0:
        print("y[{}]: {}".format(a, y[a].x))
print("\nObjective Function Value:")
print("ObjVal: {}".format(mdl.ObjVal))
    
kk_expr = (quicksum(y_val[l,k,6] for l in K_u for k in K_d) + quicksum(y_val[l,k,8] for l in K_u for k in K_d) +
                 quicksum(y_val[l,k,14] for l in K_d for k in K_u) + quicksum(y_val[l,k,16] for l in K_d for k in K_u))
print("Value of kk:", kk_expr.getValue())  

# Performance metrics: PLF, Congestion Duration, Left-Behind Rate, D-to-C
eps = 1e-6
cap = C
# 1) Peak Load Factor (PLF) = max_{k,i} n1[k,i] / C1
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

solution_count = 1  # Initialize solution count            
while mdl.Status == GRB.OPTIMAL:
    mdl.addConstr(quicksum(y[l,k,6] for l in K_u for k in K_d) + quicksum(y[l,k,8] for l in K_u for k in K_d) +
                 quicksum(y[l,k,14] for l in K_d for k in K_u) + quicksum(y[l,k,16] for l in K_d for k in K_u) >= kk_expr+1)
    mdl.setParam('OutputFlag', 0)
    mdl.optimize()
    y_val = {} 
    for a in A:
        y_val[a] = y[a].x
        if y[a].x > 0:
            print("y[{}]: {}".format(a, y[a].x))
    
    kk_expr = (quicksum(y_val[l,k,6] for l in K_u for k in K_d) + quicksum(y_val[l,k,8] for l in K_u for k in K_d) +
                 quicksum(y_val[l,k,14] for l in K_d for k in K_u) + quicksum(y_val[l,k,16] for l in K_d for k in K_u))
    # Print the value of kk_expr
    print("\nObjective Function Value:")
    print("ObjVal: {}".format(mdl.ObjVal))
    print("Value of kk:", kk_expr.getValue()) 
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