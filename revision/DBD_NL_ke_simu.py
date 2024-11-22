import pandas as pd
import numpy as np
import ADP_NL_cases as cases
import cvxpy as cp
import math
from cvxpy.error import SolverError
import mosek
import os

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

def read(c,choice):

    V = np.full((c,c*2 + 1, 2), 0.0)
    folder_path = "DBD_NL_ke/seat" + str(c) + "choice" + str(choice)
    for i in range(c):
        file_name = "DBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo" + str(i) + ".txt"
        file_path = f"{folder_path}/{file_name}"
        V[i,:,:] = np.loadtxt(file_path)
    folder_path = "DBD_NL/"
    file_name = "DBD_NL_DPtable" + str(choice) + "capacity" + str(c)  + ".txt"
    file_path = f"{folder_path}/{file_name}"
    Vy = np.loadtxt(file_path)

    return V, Vy

def choose_choice(choice,c):

    if choice == 4:
        a1, a2, a3, b, tau = cases.hom_de_in_seats(c)
    elif choice == 5:
        a1, a2, a3, b, tau = cases.hom_hom_in_seats(c)
    return a1, a2, a3, b, tau


def findseats(rand,p1,pj):
    sum=p1
    for i in range(len(pj)):
        sum=sum+pj[i]
        if sum>=rand:
            return i
    raise TypeError("Math Wrong")

def Simulation(random):

    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")
    x = np.full(c,1)
    y = c

    Revenue=0
    V, Vy=read(c,choice)

    for t in range(T,0,-1):
        if y==0:
            break
        Vyd = 1000
        Vyd2=1000
        Vxd = [1000] * c
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
        if random[t - 1] > p10 + p20 + p30:
            continue
        elif random[t - 1] < p10:
            Revenue = Revenue + price1
            y = y - 1
        elif random[t - 1] < p10 + p20:
            seat = findseats(random[t - 1], p10, p2j0)
            if x[seat] == 0:
                continue
            else:
                Revenue = Revenue + price2[seat]
                y = y - 1
                x[seat] = 0
        else:
            seat = findseats(random[t - 1], p10 + p20, p3j0)
            if x[seat] == 0:
                continue
            else:
                if x[seat] == 0 or x[seat - (-1) ** (seat + 1)] == 0:
                    continue
                else:
                    y = y - 2
                    Revenue = Revenue + price3[seat]
                    x[seat] = 0
                    x[seat - (-1) ** (seat + 1)] = 0
    return Revenue


num_sim=100
mean_size=np.zeros((10,6))
var_size=np.zeros((10,6))
for j in range(1,6):
    c=j*8
    T=c*2
    folder_path = "RandomNumbers/"+str(c)
    file_name = "rand.txt"
    file_path = f"{folder_path}/{file_name}"
    random = np.loadtxt(file_path)
    for choice in range(4,6):
        a1, a2, a3, b, tau=choose_choice(choice,c)
        results=[0]*num_sim
        for i in range(num_sim):
            rand=random[:,i]
            results[i-num_sim]=Simulation(rand)

        directory = "results/DBD_ke"
        os.makedirs(directory, exist_ok=True)
        full_path = directory + "/" + str(choice) + str(c) + "CRN.txt"
        np.savetxt(full_path, results)

        mean_size[j-1,choice]=np.mean(results)
        var_size[j-1,choice]=np.var(results)
        print("c,choice:",c,choice)
        print("mean,var:",mean_size[j-1,choice],var_size[j-1,choice])