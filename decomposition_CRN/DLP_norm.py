import cvxpy as cp
import numpy as np
import Branch_Bound as BB
import ADP_NL_cases as cases
import math
import time

def UP_t(c, a1, a2, a3, b, tau):
    result = 1 - 1 / (1 + math.exp(a1 * tau[0]) + sum(math.exp(a2[j]) for j in range(c)) ** tau[1]
                      + sum(math.exp(a3[j]) for j in range(c)) ** tau[2])
    return result

def DLP_UB(c):
    start=time.time()
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
        constraints = constraints + [T*(p2j[j]+p3j[j]+p3j[j - (-1) ** (j + 1)])<=1]

    constraints = constraints+[T*(p1+p2+p3*2)<=c]
    constraints=constraints+[p2j>=0,
                   p3j>=0,
                   p1>=0,
                   p0+p1+p2+p3==1,
                   p2 == cp.sum(p2j),
                   p3 == cp.sum(p3j),
                   p0 >= 1 - UP_t(c, a1, a2, a3, b, tau)]

    subp = cp.Problem(objective, constraints)

    subp.solve(solver=cp.MOSEK)
    print(subp.value)
    end=time.time()
    dual=[0]*(c+1)
    for j in range(c):
        dual[j]=constraints[j].dual_value
    d_temp=constraints[c].dual_value
    dual[c]=d_temp[0]

    folder_path = "SBD_NL"
    file_name = "DLPnorm" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, dual)

    return end-start


eps = 1e-5
results=[0]*11
for choice in range(3):
    for i in range(1,11):
        c=i*8
        T = c * 2
        results[i]=DLP_UB(c)
    print(results)