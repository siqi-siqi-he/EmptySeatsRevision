import cvxpy as cp
import numpy as np
import Branch_Bound as BB
import ADP_NL_cases as cases
import math
import time

c=4
T=c*2
choice=2
eps = 1e-5

if choice == 0:
    a1, a2, a3, b, tau = cases.homo_seats(c)
elif choice == 1:
    a1, a2, a3, b, tau = cases.incre_seats(c)
elif choice == 2:
    a1, a2, a3, b, tau = cases.aw_seats(c)

p0 = cp.Variable(1, name="p0")
p1 = cp.Variable(1, name="p1")
p2 = cp.Variable(1, name="p2")
p3 = cp.Variable(1, name="p3")
p2j = cp.Variable(c, name="p2j")
p3j = cp.Variable(c, name="p3j")
z = cp.Variable(c, name="z")

objective_cp = T*(1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1  \
                   + cp.sum(
        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j]
         for j in range(c)]) \
                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                   + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j]  for j in range(c)]) \
                   + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0))

objective = cp.Maximize(objective_cp)

constraints=[]
for j in range(c):
    constraints = constraints + [T*(p2j[j]+p3j[j]+p3j[j - (-1) ** (j + 1)])+z[j]<=1]

constraints = constraints+[p2j>=eps,
               p3j>=eps,
               z>=eps,
               p1>=eps,
               T*p1<=cp.sum(z),
               p0+p1+p2+p3==1,
               p2 == cp.sum(p2j),
               p3 == cp.sum(p3j)]

subp = cp.Problem(objective, constraints)

subp.solve(solver=cp.MOSEK)
dual=[0]*c
for j in range(c):
    print(constraints[j].dual_value)
    dual[j]=constraints[j].dual_value

print(constraints[c].dual_value)
print(subp.value)

folder_path = "SBD_NL"
file_name = "DLPflex" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, dual)

