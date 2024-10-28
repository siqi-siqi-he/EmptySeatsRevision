import cvxpy as cp
import numpy as np
import Branch_Bound as BB
import ADP_NL_cases as cases
import math
import time
import os

def UP_t(c, a1, a2, a3, b, tau):
    result = 1 - 1 / (1 + math.exp(a1 * tau[0]) + sum(math.exp(a2[j]) for j in range(c)) ** tau[1]
                      + sum(math.exp(a3[j]) for j in range(c)) ** tau[2])
    return result

def read(c,choice):
    folder_path = "BB_test"
    file_name = "v_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    v = np.loadtxt(file_path)
    file_name = "w_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    w = np.loadtxt(file_path)
    file_name = "theta_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    theta = np.loadtxt(file_path)

    return v,w,theta

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
    #price of product 3 given purchase probability p1 p2 p3
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

def obj_func1(v, w, theta, t, p0, p10, p20, p30, x, y, V):
    #obj when product 3 is not selling
    obj_value = p10 * r1(p10, p20, p30) +V[t-1,y-1] * p10 + sum(
        [p20[j] * (r2(j, p10, p20, p30) - v[t - 1, j] +V[t-1,y-1]) for j in range(c)]) \
                +p0*V[t-1,y] - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
    return obj_value

def obj_func2(v, w, theta, t, p0, p10, p20, p30, x, y, V):
    #obj when product 1, 2 and 3 are all selling
    obj_value = p10 * r1(p10, p20, p30) +V[t-1,y-1] * p10 + sum(
        [p20[j] * (r2(j, p10, p20, p30) - v[t - 1, j] +V[t-1,y-1]) for j in range(c)]) \
                + sum(
        [p30[j] * (r3(j, p10, p20, p30) - v[t - 1, j] - v[t - 1, j - (-1) ** (j + 1)] +V[t-1,y-2]) for j in
         range(c)]) \
                +p0*V[t-1,y] - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
    return obj_value

def computeDP(choice):
    #compute the DP table
    V = np.full((T+1, c+1), 0.0)
    eps = 1e-5
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    x = cp.Variable(c, name="x")
    v,w,theta=read(c,choice)
    start=time.time()

    for t in range(1,T+1):
        #at different time horizon t
        V[t, 0] = V[t - 1, 0] + sum((v[t, j] - v[t - 1, j]) for j in range(c))
        for y in range(1,c+1):
            if y==1:
                #product 3 cannot be sold when y==1
                objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1+p1*V[t-1,y-1] \
                               + cp.sum(
                    [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] +(V[t-1,y-1]-v[t - 1, j]) * p2j[j]
                     for j in range(c)]) \
                               + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                               +p0*V[t-1,y]- sum((v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
                constraints = [p1 <= 1,
                               p1 >= eps,
                               p0 >= eps,
                               p2 >= eps,
                               p2 <= 1,
                               p1 + p0 + p2 == 1,
                               p2 == cp.sum(p2j),
                               p0 >= 1 - UP_t(c, a1, a2, a3, b, tau),
                               x >= 0,
                               x <= 1,
                               p2j <= 1,
                               p2j >= eps / c,
                               p2j <= x + eps / c,
                               p3j == 0]

                bool_vars = [x[i] for i in range(c)]
                vars = {p0[0], p2[0], p3[0]}
                func = lambda p0, p10, p20, p30, x: obj_func1(v, w, theta, t, p0, p10, p20, p30, x, y, V)
                root = BB.BBNode(vars=vars, constraints=constraints, objective=objective_cp, bool_vars=bool_vars,
                                     p0_vars={p0}, p1_vars={p1},
                                     p2j_vars={p2j[j] for j in range(c)}, p3j_vars={p3j[j] for j in range(c)}, func=func)
                res, sol = root.bbsolve()
                V[t, 1] = res

            if y>1:
                #product 3 can be considered
                objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 + p1 * V[t - 1, y - 1] \
                               + cp.sum(
                    [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (v[t - 1, j]) * p2j[j] + V[
                        t - 1, y - 1] * p2j[j] for j in range(c)]) \
                               + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                               + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (
                        v[t - 1, j] + v[t - 1, j - (-1) ** (j + 1)] - V[t - 1, y - 2]) * p3j[j] for j in range(c)]) \
                               + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) \
                               + p0 * V[t - 1, y] - sum((v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
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
                               x >= 0,
                               x <= 1,
                               p2j <= 1,
                               p3j <= 1,
                               p2j >= eps / c,
                               p3j >= eps / c,
                               p2j <= x + eps / c,
                               p3j <= x + eps / c]

                for j in range(c):
                    constraints = constraints + [p3j[j] <= x[j - (-1) ** (j + 1)] + eps / c]
                bool_vars = [x[i] for i in range(c)]
                vars = {p0[0], p2[0], p3[0]}
                func = lambda p0, p10, p20, p30, x: obj_func2(v, w, theta, t, p0, p10, p20, p30, x, y, V)
                root = BB.BBNode(vars=vars, constraints=constraints, objective=objective_cp, bool_vars=bool_vars,
                                     p0_vars={p0}, p1_vars={p1},
                                     p2j_vars={p2j[j] for j in range(c)}, p3j_vars={p3j[j] for j in range(c)}, func=func)
                res, sol = root.bbsolve()
                V[t, y] = res

    end=time.time()
    #print(V)
    print('total time',end-start)
    UB_DBD_y = V[T, c] + sum(v[T, :])
    print("Upper Bound:", UB_DBD_y)
    folder_path = "DBD_NL"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_name = "DBD_NL_DPtable" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    np.savetxt(file_path, V)

for i in range(6):
    #do DBD from c=8 to c=48
    c=(i+1)*8
    T = c * 2
    print("c",c)
    choice = 0
    a1, a2, a3, b, tau = cases.homo_seats(c)

    print("choice",choice)
    computeDP(choice)