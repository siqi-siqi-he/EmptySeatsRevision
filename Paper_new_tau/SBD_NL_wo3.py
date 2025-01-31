import cvxpy as cp
import numpy as np
import new_cases as cases
import math
import time
import os

def UP_t(c, a1, a2, a3, b, tau):
    #calculate upper bound for probability variable
    result = 1 - 1 / (1 + math.exp(a1 * tau[0]) + sum(math.exp(a2[j]) for j in range(c)) ** tau[1]
                      + sum(math.exp(a3[j]) for j in range(c)) ** tau[2])
    return result

def read(c,choice):
    folder_path = "wo3_SBD_NL"
    file_name = "DLPnorm" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    w = np.loadtxt(file_path)
    return w

def r1(p1, p2, p3):
    #price of product 1 given purchase probability p1 p2 p3
    if p1 <= 0:
        result = 0
    else:

        result = a1 / b[0] + 1 / b[0] * 1 / tau[0] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p1))

    return result

def r2(j, p1, p2, p3):
    #price of product 2 given purchase probability p1 p2 p3
    if p2[j] <= 0:
        result = 0
    elif sum(p2) <= 0:
        for i in range(len(p2)):
            if p2[i] <= 0:
                p2[i] = 0
        result = a2[j] / b[1] + 1 / b[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p2[j])) + 1 / b[1] * (
                1 - tau[1]) / tau[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p2)))
    else:
        result = a2[j] / b[1] + 1 / b[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p2[j])) + 1 / b[1] * (
                    1 - tau[1]) / tau[1] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p2)))
    return result

def r3(j, p1, p2, p3):
    # price of product 3 given purchase probability p1 p2 p3
    if p3[j] <= 0:
        result = 0
    elif sum(p3) <= 0:
        for i in range(len(p3)):
            if p3[i] <= 0:
                p3[i] = 0
        result = a3[j] / b[2] + 1 / b[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p3[j])) + 1 / b[2] * (
                    1 - tau[2]) / tau[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p3)))
    else:
        result = a3[j] / b[2] + 1 / b[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p3[j])) + 1 / b[2] * (
                    1 - tau[2]) / tau[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p3)))

    return result

def computeDP():

    #Compute DP table for SBD using CVXPY.
    # Initialize value table V with zeros
    V = np.full((T+1, c+1), 0.0)

    eps = 1e-5

    # Define optimization variables
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    w=read(c,choice)

    start=time.time()

    # Dynamic programming loop over time horizon and capacity
    for t in range(1,T+1):
        V[t, 0] = V[t - 1, 0] # Base case initialization
        for y in range(1, c + 1):
            if y>=1:
            # if y==1, there will be no selling of product 3
                objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 +p1*V[t-1,y-1]\
                               + cp.sum(
                    [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] +(V[t-1,y-1]- w[j]) * p2j[j]
                     for j in range(c)]) \
                               + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0)+p0*V[t-1,y]

            
                constraints = [p1 <= 1,
                            p1 >= eps,
                            p0 >= eps,
                            p2 >= eps,
                            p2 <= 1,
                            p1 + p0 + p2 + p3 == 1,
                            p2 == cp.sum(p2j),
                            p3 == cp.sum(p3j),
                            p0 >= 1 - UP_t(c, a1, a2, a3, b, tau),
                            p2j <= 1,
                            p2j >= eps / c,
                            p3j == eps / c]
                objective = cp.Maximize(objective_cp)

                subp = cp.Problem(objective, constraints)

                try:
                    subp.solve(solver=cp.MOSEK)
                except cp.error.SolverError as e:
                    print(e)
                    try:
                        subp.solve(solver=cp.SCS)
                    except cp.error.SolverError as e:
                        print(e)
                        subp.solve(solver=cp.ECOS)

                V[t, y] = subp.value

    #after the loop for DP table ends, save the results
    end=time.time()
    #print(V)
    print('total time',end-start)
    UB_SBD_y = V[T, c] + sum(w) - w[c]
    print('Upper bounds:', UB_SBD_y)
    folder_path = "wo3_SBD_NL"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "SBD_NL_DPtable" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, V)
    return end-start, UB_SBD_y

results=[0]*11
bounds=[0]*11
for choice in range(1,2):
    print("choice:", choice)
    for i in range(1,2):
        #this means that we take capacity of bus from 8 to 80
        c=i*8
        print("c:", c)
        T=c*2
        if choice == 1:
            a1, a2, a3, b, tau = cases.homo_seats(c)
        elif choice == 3:
            a1, a2, a3, b, tau = cases.incre_seats(c)
        elif choice == 2:
            a1, a2, a3, b, tau = cases.aw_seats(c)
        results[i], bounds[i]=computeDP()
    folder_path = "wo3_results"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "runtime_SBD" + str(choice) + "capacity.txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, results)
    file_name = "UB_SBD" + str(choice) + "capacity.txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, bounds)
    #print(results)
