import cvxpy as cp
import numpy as np
import Branch_Bound_y as BB
import ADP_NL_cases as cases
import math
import time
import os
import multiprocessing as multip
c=16
T=c*2

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
    if p1 <= 0:
        result = 0
    else:
        try:
            result = a1 / b[0] + 1 / b[0] * 1 / tau[0] * (math.log(1 - p1 - sum(p2) - sum(p3)) - math.log(p1))
        except ValueError as e:
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

def obj_func0(v, w, theta, t, p0, p10, p20, p30, x, y,V, i):
    obj_value = p10 * r1(p10, p20, p30) + (V[t - 1, 0]-w[t-1]) * p10 + sum(
        [p20[j] * (r2(j, p10, p20, p30) - v[t - 1, j] -w[t-1]+ V[t - 1, 0]) for j in range(c)]) \
                + sum(
        [p30[j] * (r3(j, p10, p20, p30) - v[t - 1, j] - v[t - 1, j - (-1) ** (j + 1)]-2*w[t-1] + V[t - 1, 0]) for j in
         range(c)]) \
                + p0 * V[t - 1, 0] - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))+(v[t, i] - v[t - 1, i]) * x[i]+(w[t-1]-w[t])*(y+1)
    return obj_value

def obj_func1(v, w, theta, t, p0, p10, p20, p30, x, y,V, i):
    obj_value = p10 * r1(p10, p20, p30) + (V[t - 1, 1]-w[t-1]) * p10 + sum(
        [p20[j] * (r2(j, p10, p20, p30) - v[t - 1, j] -w[t-1]+ V[t - 1, 1]) for j in range(c)]) \
    -p20[i]*(V[t-1,1]-V[t-1,0])+v[t-1,i]*p20[i]\
                + sum(
        [p30[j] * (r3(j, p10, p20, p30) - v[t - 1, j] - v[t - 1, j - (-1) ** (j + 1)]-2*w[t-1] + V[t - 1, 1]) for j in
         range(c)]) \
                - (p30[i]+p30[i - (-1) ** (i + 1)]) * (V[t - 1, 1]-V[t - 1, 0] - v[t - 1, i] ) \
                + p0 * V[t - 1, 1] - sum(
        (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))+(v[t, i] - v[t - 1, i]) * x[i]+(w[t-1]-w[t])*(y+1)
    return obj_value

def computeDP(choice):
    UB_DBD_ke=100000
    eps = 1e-5
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    x = cp.Variable(c, name="x")
    y = cp.Variable(1, name="y")

    v,w,theta=read(c,choice)

    start=time.time()
    for i in range(c):
        start_i = time.time()
        V = np.full((T + 1, 2), 0.0)
        for t in range(1,T+1):
            for situation in range(2):
                if situation == 0:
                    objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 + (
                                V[t - 1, 0] - w[t - 1]) * p1 \
                                   + cp.sum(
                        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (
                                    v[t - 1, j] + w[t - 1]) * p2j[j] + V[t - 1, 0] * p2j[j]
                         for j in range(c)]) \
                                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                                   + cp.sum(
                        [1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (
                                    v[t - 1, j] + v[t - 1, j - (-1) ** (j + 1)] + 2 * w[t - 1]) *
                         p3j[j] + V[t - 1, 0] * p3j[j]
                         for j in range(c)]) \
                                   + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) \
                                   + p0 * V[t - 1, 1] - sum((v[t, j] - v[t - 1, j]) * x[j] for j in range(c))\
                                   +(v[t,i]-v[t-1,i])*x[i]+(w[t-1]-w[t])*(y+1)

                    constraints = [p1 <= 1,
                                   p1 >= eps,
                                   p0 >= eps,
                                   p2 >= eps,
                                   p2 <= 1,
                                   p3 >= eps,
                                   p3 <= 1,
                                   x[i]==0,
                                   p1 + p0 + p2 + p3 == 1,
                                   p2 == cp.sum(p2j),
                                   p3 == cp.sum(p3j),
                                   p0 >= 1 - UP_t(c, a1, a2, a3, b, tau),
                                   x >= 0,
                                   x <= 1,
                                   y<=1,
                                   y>=0,
                                   p1<=y+1+eps/c,
                                   p2j <= 1,
                                   p3j <= 1,
                                   p2j >= eps / c,
                                   p3j >= eps / c,
                                   p2j <= x + eps / c,
                                   p3j <= x + eps / c,
                                   p2j <= y+1 + eps / c,
                                   p3j <= y + eps / c,
                                   cp.sum(x)>=1+y]

                    for j in range(c):
                        constraints = constraints + [p3j[j] <= x[j - (-1) ** (j + 1)] + eps / c]

                    bool_vars = [x[i] for i in range(c)]
                    vars = {p0[0], p2[0], p3[0]}
                    func = lambda p0, p10, p20, p30, x, y: obj_func0(v, w, theta, t, p0, p10, p20, p30, x, y, V, i)
                    root = BB.BBNode(vars=vars, constraints=constraints, objective=objective_cp,
                                         bool_vars=bool_vars,
                                         p0_vars={p0}, p1_vars={p1},
                                         p2j_vars={p2j[j] for j in range(c)}, p3j_vars={p3j[j] for j in range(c)},
                                         func=func,y_vars=y)
                    res, sol = root.bbsolve()
                    res_temp=sum(v[t-1,:]-v[t,:])-(v[t-1,i]-v[t,i])
                    if res_temp>res:
                        res=res_temp
                    V[t, 0] = res
                else:
                    objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 +(V[t-1,1]-w[t-1])*p1\
                                   + cp.sum(
                        [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (v[t - 1, j]+w[t-1]) * p2j[j] +V[t-1,1]*p2j[j]
                         for j in range(c)]) \
                                   + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0)-(V[t-1,1]-V[t-1,0])*p2j[i]+v[t-1,i]*p2j[i] \
                    + cp.sum(
                        [1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (v[t - 1, j] +v[t - 1, j - (-1) ** (j + 1)]+ 2*w[t - 1]) *
                         p3j[j] + V[t - 1, 1] * p3j[j]
                         for j in range(c)]) \
                    + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) - (V[t - 1, 1]-V[t-1,0]-v[t - 1, i]) * (p3j[i]+p3j[i - (-1) ** (i + 1)]) \
                    +p0*V[t-1,1]- sum((v[t, j] - v[t - 1, j]) * x[j] for j in range(c))+(v[t, i] - v[t - 1, i]) * x[i]+(w[t-1]-w[t])*(y+1)
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
                                   y <= 1,
                                   y >= 0,
                                   x[i]==0,
                                   p2j <= 1,
                                   p3j <= 1,
                                   p2j >= eps / c,
                                   p3j >= eps / c,
                                   p2j <= x + eps / c,
                                   p3j <= x + eps / c,
                                   p3j <= y + eps / c,
                                   cp.sum(x)>=y+1]

                    for j in range(c):
                        constraints = constraints + [p3j[j] <= x[j - (-1) ** (j + 1)] + eps / c]

                    bool_vars = [x[i] for i in range(c)]
                    vars = {p0[0], p2[0], p3[0]}
                    func = lambda p0, p10, p20, p30, x, y: obj_func1(v, w, theta, t, p0, p10, p20, p30, x, y, V, i)
                    root = BB.BBNode(vars=vars, constraints=constraints, objective=objective_cp, bool_vars=bool_vars,
                                         p0_vars={p0}, p1_vars={p1},
                                         p2j_vars={p2j[j] for j in range(c)}, p3j_vars={p3j[j] for j in range(c)}, func=func,y_vars=y)
                    res, sol = root.bbsolve()
                    res_temp = sum(v[t - 1, :] - v[t, :]) - (v[t - 1, i] - v[t, i])
                    if res_temp > res:
                        res = res_temp
                    V[t, 1] = res

        end_i=time.time()

        folder_path = "DBD_NL_ke/seat" + str(c) + "choice" + str(choice)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_name = "DBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo"+str(i)+".txt"
        file_path = f"{folder_path}/{file_name}"
        os.makedirs(folder_path, exist_ok=True)
        np.savetxt(file_path, V)
        temp_UB = V[T, 1] + sum(v[T, :]) - v[T, i] + w[T] * c
        if temp_UB <= UB_DBD_ke:
            UB_DBD_ke = temp_UB
    end=time.time()
    total_time=end-start
    print("total time", total_time)
    print("UB for x (y to be tested):", UB_DBD_ke)
    return total_time

for i in range(1,4):
    c=i*8
    T = c * 2
    print("c",c)
    choice=1
    if choice == 0:
        a1, a2, a3, b, tau = cases.homo_seats(c)
    elif choice == 1:
        a1, a2, a3, b, tau = cases.incre_seats(c)
    elif choice == 2:
        a1, a2, a3, b, tau = cases.aw_seats(c)
    print("choice",choice)
    computeDP(choice)