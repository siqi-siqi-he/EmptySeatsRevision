import pandas as pd
import numpy as np
import sys
import cplex
from cplex.exceptions import CplexSolverError
from cplex import SparsePair
import docplex
from docplex.mp.model import Model
import cases as cases
import cvxpy as cp
import math
import time
from cvxpy.error import SolverError
import mosek
choice=0
c=8
eps=1e-5

def UP_t(c,a1, a2, a3, b, tau):
    result=1-1/(1+math.exp(a1*tau[0])+sum(math.exp(a2[j]) for j in range(c))**tau[1]
          +sum(math.exp(a3[j]) for j in range(c))**tau[2])
    return result

def UP_1(c,a1, a2, a3, b, tau):
    result=math.exp(tau[0]*(a1))/(math.exp(tau[0]*(a1))+1)
    return result

def UP_2(j, c, a1, a2, a3, b, tau):
    result = math.exp(a2[j]) * sum(math.exp(a2[j]) for j in range(c)) ** (tau[1] - 1) / (
                1 + sum(math.exp(a2[j]) for j in range(c)) ** (tau[1]))
    return result

def UP_3(j, c, a1, a2, a3, b, tau):
    result = math.exp(a3[j]) * sum(math.exp(a3[j]) for j in range(c)) ** (tau[2] - 1) / (
            1 + sum(math.exp(a3[j]) for j in range(c)) ** (tau[2]))
    return result

def V(x, y, t,c,choice):
    v, w, theta = read(c, choice)
    V = sum(x[j] * v[t][j] for j in range(c)) + w[t] * y + theta[t]
    return V

def E(j,c):
    result=[0 for j in range(c)]
    result[j]=1
    return result

def r1(p1,p2,p3,a1, a2, a3, b, tau):
    if p1<=0:
        result=0
    else:
        result=a1/b[0]+1/b[0]*1/tau[0]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p1))
    return result

def r2(j,p1,p2,p3,a1, a2, a3, b, tau):
    if p2[j]<=0:
        result=0
    else:
        result=a2[j]/b[1]+1/b[1]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p2[j]))+1/b[1]*(1-tau[1])/tau[1]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(sum(p2)))
    return result

def r3(j,p1,p2,p3,a1, a2, a3, b, tau):
    if p3[j]<=0:
        result=0
    else:
        result=a3[j]/b[2]+1/b[2]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p3[j]))+1/b[2]*(1-tau[2])/tau[2]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(sum(p3)))

    return result

def read(c,choice):
    folder_path = "BB"
    file_name = "v_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    v = np.loadtxt(file_path)
    file_name = "w_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    w = np.loadtxt(file_path)
    file_name = "theta_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    theta = np.loadtxt(file_path)

    return v,w,theta


def Price_Calculate():

    T=c*2

    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    #maybe other variables
    folder_path = "BB"
    file_name = "v_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    v = np.loadtxt(file_path)
    file_name = "w_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    w = np.loadtxt(file_path)
    file_name = "theta_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
    file_path = f"{folder_path}/{file_name}"
    theta = np.loadtxt(file_path)

    x = np.array([0,0,1,1,0,0,0,0])
    y = 2
    price1=[0]*(T+1)
    price2=[0]*(T+1)
    price3=[0]*(T+1)

    for t in range(T,0,-1):
        print('t',t)
        print(V(x, y, t,c,choice))
        objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 - w[t - 1] * p1 \
                    + cp.sum(
            [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (v[t - 1, j] + w[t - 1]) * p2j[j]
            for j in range(c)]) \
                    + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                    + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (
                    v[t - 1, j] + v[t - 1, j - (-1) ** (j + 1)] + 2 * w[t - 1]) * p3j[j] for j in range(c)]) \
                    + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) \
                    - sum(
            (v[t, j] - v[t - 1, j]) * x[j] for j in range(c))
        # objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 + V(x, y - 1, t - 1, c,
        #                                                                                          choice) * p1 \
        #                + cp.sum(
        #     [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] + V(x - E(j, c),
        #                                                                                    y - 1, t - 1,
        #                                                                                    c, choice) *
        #      p2j[j] for j in range(c)]) \
        #                + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
        #                + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] + V(
        #     x - E(j, c) - E(j - (-1) ** (j + 1), c), y - 2, t - 1, c, choice) * p3j[j] for j in range(c)]) \
        #                + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0) + p0 * V(x, y, t - 1, c,
        #                                                                                             choice)
        
        objective = cp.Maximize(objective_cp)

        constraints = [p1 <= 1,
                       p1 >= eps,
                       p0 >= eps,
                       p0 <= 1,
                       p2 >= eps,
                       p2 <= 1,
                       p3 >= eps,
                       p3 <= 1,
                       p1 + p0 + p2 + p3 == 1,
                       p2 == cp.sum(p2j),
                       p3 == cp.sum(p3j),
                       p0>= 1-UP_t(c, a1, a2, a3, b, tau)]
        for j in range(c):
            constraints = constraints + [p2j[j] <= x[j]+eps/c,
                                         p3j[j] <= x[j]+eps/c,
                                         p3j[j] <= x[j - (-1) ** (j + 1)]+eps/c,
                                         p1 <= y+eps,
                                         p2j[j] <= y+eps/c,
                                         p3j[j] <= math.floor(y / 2)+eps/c,
                                         p2j[j] <= 1,
                                         p3j[j] <= 1,
                                         p2j[j] >= eps/c,
                                         p3j[j] >= eps/c]
        subp = cp.Problem(objective, constraints)
        try:
            subp.solve(solver=cp.MOSEK)
            if math.isnan(subp.value):
                subp.solve(solver=cp.SCS)
            pass
        except cp.error.SolverError as e:
            # print(f"SolverError: {e}")
            subp.solve(solver=cp.SCS)
        if p0[0].value is None or p1[0].value is None or p2[0].value is None or p3[0].value is None:
            subp.solve(solver=cp.SCS)


        if p0[0].value is None or p1[0].value is None or p2[0].value is None or p3[0].value is None:
            raise ValueError("None again, tried=", p0[0].value ,p1[0].value,p2[0].value,p3[0].value)
        p00 = p0[0].value
        p10 = p1[0].value
        p20 = p2[0].value
        p30 = p3[0].value
        p2j0 = p2j.value
        p3j0 = p3j.value
        obj_value = subp.value
        print('r1', r1(p10,p2j0,p3j0,a1, a2, a3, b, tau))
        price1[t]=r1(p10,p2j0,p3j0,a1, a2, a3, b, tau)
        price2[t]=r2(3,p10,p2j0,p3j0,a1, a2, a3, b, tau)
        print('r2j',[r2(j, p10,p2j0,p3j0,a1, a2, a3, b, tau) for j in range(c)])
        print('r2',price2[t])
        price3[t] = [r3(j, p10,p2j0,p3j0,a1, a2, a3, b, tau) for j in range(c)][3]
        print('r3',price3[t])

    return price1, price2, price3



c = 8
choice=1
a1, a2, a3, b, tau = cases.homo_seats(c)
out_p1=[]
out_p2=[]
out_p3=[]
for step in range(10,11):
    a3=[0.1*step for i in range(len(a3))]
    p1,p2,p3=Price_Calculate()
    out_p1.append(p1[1:])
    out_p2.append(p2[1:])
    out_p3.append(p3[1:])
    
df1 = pd.DataFrame(out_p1)
df2 = pd.DataFrame(out_p2)
df3 = pd.DataFrame(out_p3)
output_file = "ADP_prices.xlsx"  # Replace with your desired file path
sheet1 = "price1"  # Replace with your desired sheet name
sheet2 = "price2"
sheet3 = "price3"
#Save to a specific sheet
with pd.ExcelWriter(output_file, mode='w', engine='openpyxl') as writer:
    df1.to_excel(writer, sheet_name=sheet1, index=False, header=True)
    df2.to_excel(writer, sheet_name=sheet2, index=False, header=True)
    df3.to_excel(writer, sheet_name=sheet3, index=False, header=True)