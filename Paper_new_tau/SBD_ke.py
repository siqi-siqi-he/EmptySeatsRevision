import cvxpy as cp
import numpy as np
import new_cases as cases
import math
import time
import os

def UP_t(c, a1, a2, a3, b, tau):
    result = 1 - 1 / (1 + math.exp(a1 * tau[0]) + sum(math.exp(a2[j]) for j in range(c)) ** tau[1]
                      + sum(math.exp(a3[j]) for j in range(c)) ** tau[2])
    return result

def read(c,choice):
    folder_path = "SBD_NL"
    file_name = "DLPnorm" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    w = np.loadtxt(file_path)
    return w

def r1(p1, p2, p3):
    if p1 <= 0:
        result = 0
    else:
        try:
            result = a1 / b[0] + 1 / b[0] * 1 / tau[0] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p1))
        except ValueError as e:
            # Print the variable causing the ValueError
            print(f"ValueError: {e}")
            print("Printing variable values:")
            print(f"p1: {p1}")
            print(f"p2: {p2}")
            print(f"p3: {p3}")
            print(f"tau[2]: {tau[2]}")
            result = e
    return result

def r2(j, p1, p2, p3):
    try:
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
    except ValueError as e:
        print(f"ValueError: {e}")
        print("Printing variable values:")
        print(f"a2[j]: {a2[j]}")
        print(f"b[2]: {b[2]}")
        print(f"p1: {p1}")
        print(f"p2: {p2}")
        print(f"p3: {p3}")
        print(f"tau[1]: {tau[1]}")
        print(f"j: {j}")
    return result

def r3(j, p1, p2, p3):
    if p3[j] <= 0:
        result = 0
    elif sum(p3) <= 0:
        for i in range(len(p3)):
            if p3[i] <= 0:
                p3[i] = 0
        result = a3[j] / b[2] + 1 / b[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p3[j])) + 1 / b[2] * (
                    1 - tau[2]) / tau[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p3)))
    else:
        try:
            result = a3[j] / b[2] + 1 / b[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p3[j])) + 1 / b[2] * (
                        1 - tau[2]) / tau[2] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(sum(p3)))

        except ValueError as e:
            print(f"ValueError: {e}")
            print("Printing variable values:")
            print(f"a3[j]: {a3[j]}")
            print(f"b[2]: {b[2]}")
            print(f"p1: {p1}")
            print(f"p2: {p2}")
            print(f"p3: {p3}")
            print(f"tau[2]: {tau[2]}")
            print(f"j: {j}")

    return result


def computeDP(choice):
    eps = 1e-5

    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")


    w=read(c,choice)
    start=time.time()
    UB_SBD_ke = 100000
    for i in range(c):
        start_i = time.time()
        #print("component:", i)
        V = np.full((T + 1, 2), 0.0)
        for t in range(1,T+1):
            #print("time period:", t)
            for j in range(2):
                if j ==0:
                    objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 + (
                                V[t - 1, 0] - w[c]) * p1 \
                                   + cp.sum(
                        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (w[j] + w[c]) *
                         p2j[j] + V[t - 1,0] * p2j[j]
                         for j in range(c)]) \
                                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                                   + cp.sum(
                        [1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (
                                    w[j] + w[j - (-1) ** (j + 1)] + 2 * w[c]) *
                         p3j[j] + V[t - 1, 0] * p3j[j]
                         for j in range(c)]) \
                                   + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) \
                                   + p0 * V[t - 1, 0]
                    constraints = [p1 <= 1,
                                   p1 >= eps,
                                   p0 >= eps,
                                   p2 >= eps,
                                   p2 <= 1,
                                   p3 >= eps,
                                   p3 <= 1,
                                   p1 + p0 + p2 + p3 == 1,
                                   p2 == cp.sum(p2j),
                                   p3 == cp.sum(p3j),
                                   p0 >= 1 - UP_t(c, a1, a2, a3, b, tau),
                                   p2j <= 1,
                                   p3j <= 1,
                                   p2j >= eps / c,
                                   p3j >= eps / c,
                                   p2j[i]==eps/c,
                                   p3j[i]==eps/c,
                                   p3j[i - (-1) ** (i + 1)]==eps/c]
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

                    V[t, 0] = subp.value

                else:
                    objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 +(V[t-1,1]-w[c])*p1\
                                   + cp.sum(
                        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (w[j]+w[c]) * p2j[j] +V[t-1,1]*p2j[j]
                         for j in range(c)]) \
                                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0)-(V[t-1,1]-V[t-1,0])*p2j[i]+w[i]*p2j[i] \
                    + cp.sum(
                        [1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (w[j] +w[j - (-1) ** (j + 1)]+ 2*w[c]) *
                         p3j[j] + V[t - 1, 1] * p3j[j]
                         for j in range(c)]) \
                    + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) - (V[t - 1, 1]-V[t-1,0]-w[i]) * (p3j[i]+p3j[i - (-1) ** (i + 1)]) \
                    +p0*V[t-1,1]

                    constraints = [p1 <= 1,
                                   p1 >= eps,
                                   p0 >= eps,
                                   p2 >= eps,
                                   p2 <= 1,
                                   p3 >= eps,
                                   p3 <= 1,
                                   p1 + p0 + p2 + p3 == 1,
                                   p2 == cp.sum(p2j),
                                   p3 == cp.sum(p3j),
                                   p0 >= 1 - UP_t(c, a1, a2, a3, b, tau),
                                   p2j <= 1,
                                   p3j <= 1,
                                   p2j >= eps / c,
                                   p3j >= eps / c]


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

                    V[t, 1] = subp.value

        end_i=time.time()
        #print(V)
        #print('total time at i',end_i-start_i)

        folder_path = "SBD_NL_ke/seat" + str(c) + "choice" + str(choice)
        file_name = "SBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo"+str(i)+".txt"
        file_path = f"{folder_path}/{file_name}"
        os.makedirs(folder_path, exist_ok=True)
        np.savetxt(file_path, V)
        temp_UB = V[T, 1] + sum(w) - w[i] + w[c] * (c - 1)
        if temp_UB <= UB_SBD_ke:
            UB_SBD_ke = temp_UB


    end=time.time()
    print("total time:",end-start)
    print("UB for x (y to be tested):", UB_SBD_ke)
    return end-start, UB_SBD_ke

results=[0]*11
bounds=[0]*11
for choice in range(3,4):
    print("choice:", choice)
    for i in range(9,11):
        c=i*8
        T=c*2
        print("c:", c)
        if choice == 1:
            a1, a2, a3, b, tau = cases.homo_seats(c)
        elif choice == 3:
            a1, a2, a3, b, tau = cases.incre_seats(c)
        elif choice == 2:
            a1, a2, a3, b, tau = cases.aw_seats(c)
        results[i], bounds[i]=computeDP(choice)
    folder_path = "results"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "runtime_SBD_ke" + str(choice) + "capacity.txt"
    file_path = f"{folder_path}/{file_name}"
    #np.savetxt(file_path, results)
    file_name = "UB_SBD_ke" + str(choice) + "capacity.txt"
    file_path = f"{folder_path}/{file_name}"
    #np.savetxt(file_path, bounds)
    print(results)