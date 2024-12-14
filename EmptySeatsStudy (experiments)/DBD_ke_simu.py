import pandas as pd
import numpy as np

import cases as cases
import cvxpy as cp
import math
from cvxpy.error import SolverError
import mosek
import time
eps=1e-9
c=8
choice=0

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


def findseats(rand,p1,pj):
    sum=p1
    for i in range(len(pj)):
        sum=sum+pj[i]
        if sum>=rand:
            return i
    raise TypeError("Math Wrong")

def read(c,choice):

    V = np.full((c,c*2 + 1, 2), 0.0)
    folder_path = "EmptySeatsStudy (experiments)/DBD"
    for i in range(c):
        file_name = "DBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo" + str(i) + "step" + str(step) + ".txt"
        file_path = f"{folder_path}/{file_name}"
        V[i,:,:] = np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/DBD"
    file_name = "DBD_NL_DPtable" + str(choice) + "capacity" + str(c) + "step" + str(step) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    Vy = np.loadtxt(file_path)

    return V, Vy

def Simulation(ii):
    p3num=0
    T=c*2
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    #maybe other variables

    x = np.full(c,1)
    y = c

    Revenue=0
    V, Vy=read(c,choice)

    for t in range(T,0,-1):

        folder_path = "EmptySeatsStudy (experiments)/RandomNumbers"
        file_name = "randomdecisions_capacity8_choice0_sim100_a3_4_t" + str(t) + ".txt"
        file_path = f"{folder_path}/{file_name}"
        random = np.loadtxt(file_path)

        if y==0:
            break
        Vyd = 0
        Vyd2=0
        Vxd = [0] * c
        if y>0:
            Vyd=Vy[t-1,y]-Vy[t-1,y-1]
            if y>1:
                Vyd2 = Vy[t - 1, y] - Vy[t - 1, y - 2]
        for i in range(c):
            if x[i]>0:
                Vxd[i]=V[i,t-1,1]-V[i,t-1,0]
        objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 - p1*Vyd \
                       + cp.sum(
            [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (Vxd[j]+Vyd) *
             p2j[j] for j in range(c)]) \
                       + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                       + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] -(Vxd[j]+Vxd[j - (-1) ** (j + 1)]+Vyd2) * p3j[j] for j in range(c)]) \
                       + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0)
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

        price1=r1(p10,p2j0,p3j0,a1, a2, a3, b, tau)
        price2=[r2(j,p10,p2j0,p3j0,a1, a2, a3, b, tau) for j in range(c)]
        price3 = [r3(j, p10,p2j0,p3j0,a1, a2, a3, b, tau) for j in range(c)]

        if random[ii] > p10 + p20 + p30:
            continue
        elif random[ii] < p10:
            Revenue = Revenue + price1
            y = y - 1
            continue
        elif random[ii] < p10 + p20:
            seat = findseats(random[ii], p10, p2j0)
            if x[seat] == 0:
                continue
                # raise ValueError("Algorithm is incorrect: x is already sold")
            else:
                Revenue = Revenue + price2[seat]
                y = y - 1
                x[seat] = 0
                continue
        else:
            seat = findseats(random[ii], p10 + p20, p3j0)
            if x[seat] == 0:
                continue
                # raise ValueError("Algorithm is incorrect: x is already sold")
            else:

                if x[seat] == 0 or x[seat - (-1) ** (seat + 1)] == 0:
                    continue
                    # raise ValueError("Algorithm is incorrect: x is already sold")
                else:
                    y = y - 2
                    Revenue = Revenue + price3[seat]
                    x[seat] = 0
                    x[seat - (-1) ** (seat + 1)] = 0
                    p3num=p3num+1
                    continue

    return Revenue, p3num

def run_Simulator():
    start_time = time.time()
    results = [0] * 100
    p3 = [0] * 100
    for i in range(100):
        results[i],p3[i] = Simulation(i)

    end_time = time.time()

    mean = np.mean(results)
    var = np.var(results)
    p3n=np.mean(p3)
    print(results)
    print('average:', mean)
    print('variance:', var)
    print('p3 average', p3n)
    print('time:', end_time - start_time)
    return mean, var, p3n

c = 8
choice=0
a1, a2, a3, b, tau = cases.homo_seats(c)
means=[0]*51
vars=[0]*51
p3=[0]*51
for step in range(51):
    a3=[0.1*step for i in range(len(a3))]
    means[step], vars[step], p3[step]=run_Simulator()

folder_path = "EmptySeatsStudy (experiments)/DBD"
file_name = "ke_simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, means)
file_name = "ke_simu_vars_choice" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, vars)
file_name = "ke_simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, p3)