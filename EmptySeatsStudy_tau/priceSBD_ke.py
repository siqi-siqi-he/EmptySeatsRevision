import pandas as pd
import numpy as np

import cases as cases
import cvxpy as cp
import math
from cvxpy.error import SolverError
import mosek
import time


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

def read(c, choice):
    V = np.full((c, c * 2 + 1, 2), 0.0)

    folder_path = "SBD_ke"
    for i in range(c):
        file_name = "SBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo" + str(i) + "step" + str(
            step) + ".txt"
        file_path = f"{folder_path}/{file_name}"
        V[i, :, :] = np.loadtxt(file_path)

    folder_path = "SBD"
    file_name = "SBD_NL_DPtable" + str(choice) + "capacity" + str(c) + "step" + str(step) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    Vy = np.loadtxt(file_path)

    return V, Vy


def Price_Calculate():

    T=c*2
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")

    x = np.array([1,1,0,0,0,0,0,0])
    y = 2

    V, Vy = read(c, choice)
    price1=[0]*(T+1)
    price2=[0]*(T+1)
    price3=[0]*(T+1)
    for t in range(T,0,-1):

        Vyd = 0
        Vyd2 = 0
        Vxd = [0] * c
        if y > 0:
            Vyd = Vy[t - 1, y] - Vy[t - 1, y - 1]
            if y > 1:
                Vyd2 = Vy[t - 1, y] - Vy[t - 1, y - 2]
        for i in range(c):
            if x[i] > 0:
                Vxd[i] = V[i, t - 1, 1] - V[i, t - 1, 0]

        objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 - p1 * Vyd \
                       + cp.sum(
            [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (Vxd[j] + Vyd) *
             p2j[j] for j in range(c)]) \
                       + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                       + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] - (
                    Vxd[j] + Vxd[j - (-1) ** (j + 1)] + Vyd2) * p3j[j] for j in range(c)]) \
                       + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0)
        objective = cp.Maximize(objective_cp)

        if y == 0:
            break

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
                       p0 >= 1 - UP_t(c, a1, a2, a3, b, tau)]
        for j in range(c):
            constraints = constraints + [p2j[j] <= x[j] + eps / c,
                                         p3j[j] <= x[j] + eps / c,
                                         p3j[j] <= x[j - (-1) ** (j + 1)] + eps / c,
                                         p1 <= y + eps,
                                         p2j[j] <= y + eps / c,
                                         p3j[j] <= math.floor(y / 2) + eps / c,
                                         p2j[j] <= 1,
                                         p3j[j] <= 1,
                                         p2j[j] >= eps / c,
                                         p3j[j] >= eps / c]
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
            raise ValueError("None again, tried=", p0[0].value, p1[0].value, p2[0].value, p3[0].value)
        p00 = p0[0].value
        p10 = p1[0].value
        p20 = p2[0].value
        p30 = p3[0].value
        p2j0 = p2j.value
        p3j0 = p3j.value
        obj_value = subp.value

        price1[t]=r1(p10,p2j0,p3j0,a1, a2, a3, b, tau)
        price2[t]=[r2(j,p10,p2j0,p3j0,a1, a2, a3, b, tau) for j in range(c)][0]
        price3[t] = [r3(j, p10,p2j0,p3j0,a1, a2, a3, b, tau) for j in range(c)][0]

    return price1, price2, price3



c = 8
choice=0
a1, a2, a3, b, tau = cases.homo_seats(c)
out_p1=[]
out_p2=[]
out_p3=[]
for step in range(6,7):
    a3=[0.1*step for i in range(len(a3))]
    p1,p2,p3=Price_Calculate()
    out_p1.append(p1[1:])
    out_p2.append(p2[1:])
    out_p3.append(p3[1:])
    
df1 = pd.DataFrame(out_p1)
df2 = pd.DataFrame(out_p2)
df3 = pd.DataFrame(out_p3)
output_file = "SBD_ke_prices.xlsx"  # Replace with your desired file path
sheet1 = "price1"  # Replace with your desired sheet name
sheet2 = "price2"
sheet3 = "price3"
# Save to a specific sheet
with pd.ExcelWriter(output_file, mode='w', engine='openpyxl') as writer:
    df1.to_excel(writer, sheet_name=sheet1, index=False, header=True)
    df2.to_excel(writer, sheet_name=sheet2, index=False, header=True)
    df3.to_excel(writer, sheet_name=sheet3, index=False, header=True)


# output_arrays = []
# for i in range(10):  # Replace with your actual loop
#     array = [i, i * 2, i * 3]  # Replace this with your actual array
#     output_arrays.append(array)

# # Convert the list of arrays to a DataFrame
# df = pd.DataFrame(output_arrays)

# # Add a column to indicate the iteration (optional)
# df.insert(0, "Iteration", range(1, len(output_arrays) + 1))

# # Save the DataFrame to Excel
# output_file = "output.xlsx"  # Specify your desired file path
# df.to_excel(output_file, index=False, header=True)

# print(f"Data saved to {output_file}")