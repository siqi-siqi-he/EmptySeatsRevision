import cvxpy as cp
import numpy as np
#import Branch_Bound as BB
import Fig7_cases as cases
import math
import time
import os

def UP_t(c, a1, a2, a3, b, tau):
    result = 1 - 1 / (1 + math.exp(a1 * tau[0]) + sum(math.exp(a2[j]) for j in range(c)) ** tau[1]
                      + sum(math.exp(a3[j]) for j in range(c)) ** tau[2])
    return result

def DLP_UB(c):
    #for the problem formulation, see DPP in the paper
    start=time.time()

    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p2j = cp.Variable(c, name="p2j")

    objective_cp = T*(1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + (a1 -a0/tau[0])/ b[0] * p1\
                       + cp.sum(
            [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + (a2[j]-a0/tau[1]) / b[1] * p2j[j]
             for j in range(c)]) \
                       + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) )

    objective = cp.Maximize(objective_cp)

    constraints=[]
    for j in range(c):
        constraints = constraints + [T*(p2j[j])<=1]

    constraints = constraints+[T*(p1+p2)<=c]
    constraints=constraints+[p2j>=0,
                   p1>=0,
                   p2 == cp.sum(p2j),
                   p0+p2+p1==1,
                   p0 >= 0]

    subp = cp.Problem(objective, constraints)

    subp.solve(solver=cp.MOSEK)
    print(subp.value)
    end=time.time()
    dual=[0]*(c+1)
    for j in range(c):
        dual[j]=constraints[j].dual_value
    d_temp=constraints[c].dual_value
    dual[c]=d_temp[0]
    folder_path = "revision/Fig7 Extension/wo3_DLP_NL"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "DLPnorm" + str(preference) + "a0" + str(a0) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, dual)

    return subp.value


eps = 1e-5
results=np.zeros((3,5))
c=8
T=16
a1, a2, a3, b, tau = cases.homo_seats(c)
folder_path="revision/Fig7 Extension/AveragePrices"
for preference in range(3):
    #high, medium, low
    #+-0.2
    for strength in range(5):
        a0=strength*0.5-1
        if preference==0:
            file_name="a1_a0_"+str(a0)+"_a3decr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a1=np.loadtxt(file_path)
            file_name="a2_a0_"+str(a0)+"_a3decr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a2=np.loadtxt(file_path)
            a3=[0.4]*8
        if preference==1:
            a1=0.2
            a3=[0.6]*8
            a2=[0.4]*8
        if preference==2:
            file_name="a1_a0_"+str(a0)+"_a3incr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a1=np.loadtxt(file_path)
            file_name="a2_a0_"+str(a0)+"_a3incr0.2.txt"
            file_path = f"{folder_path}/{file_name}"
            a2=np.loadtxt(file_path)
            a3=[0.8]*8
        print(a0)
        results[preference][strength]=DLP_UB(c)
    print(results)