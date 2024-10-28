import pandas as pd
import numpy as np
import ADP_NL_cases as cases
import cvxpy as cp
import math
from cvxpy.error import SolverError
import mosek
import os

eps=1e-5

def UP_t(c,a1, a2, a3, b, tau):
    #The NL upper bound for p0
    result=1-1/(1+math.exp(a1*tau[0])+sum(math.exp(a2[j]) for j in range(c))**tau[1]
          +sum(math.exp(a3[j]) for j in range(c))**tau[2])
    return result

def UP_1(c,a1, a2, a3, b, tau):
    # The NL upper bound for p1
    result=math.exp(tau[0]*(a1))/(math.exp(tau[0]*(a1))+1)
    return result

def UP_2(j, c, a1, a2, a3, b, tau):
    # The NL upper bound for p2
    result = math.exp(a2[j]) * sum(math.exp(a2[j]) for j in range(c)) ** (tau[1] - 1) / (
                1 + sum(math.exp(a2[j]) for j in range(c)) ** (tau[1]))
    return result

def UP_3(j, c, a1, a2, a3, b, tau):
    # The NL upper bound for p3
    result = math.exp(a3[j]) * sum(math.exp(a3[j]) for j in range(c)) ** (tau[2] - 1) / (
            1 + sum(math.exp(a3[j]) for j in range(c)) ** (tau[2]))
    return result

def E(j,c):
    result=[0 for j in range(c)]
    result[j]=1
    return result

def r1(p1,p2,p3,a1, a2, a3, b, tau):
    #the price of product 1 given purchase probability p1 p2 p3
    if p1<=0:
        result=0
    else:
        result=a1/b[0]+1/b[0]*1/tau[0]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p1))
    return result

def r2(j,p1,p2,p3,a1, a2, a3, b, tau):
    # the price of product 2 given purchase probability p1 p2 p3
    if p2[j]<=0:
        result=0
    else:
        result=a2[j]/b[1]+1/b[1]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p2[j]))+1/b[1]*(1-tau[1])/tau[1]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(sum(p2)))
    return result

def r3(j,p1,p2,p3,a1, a2, a3, b, tau):
    # the price of product 3 given purchase probability p1 p2 p3
    if p3[j]<=0:
        result=0
    else:
        result=a3[j]/b[2]+1/b[2]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(p3[j]))+1/b[2]*(1-tau[2])/tau[2]*(math.log(1-p1-sum(p2)-sum(p3))-math.log(sum(p3)))

    return result

def read(c,choice):
    folder_path = "SBD_NL"
    file_name = "DLPnorm" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    w = np.loadtxt(file_path)

    folder_path = "SBD_NL/"
    file_name = "SBD_NL_DPtable" + str(choice) + "capacity" + str(c)  + ".txt"
    file_path = f"{folder_path}/{file_name}"
    Vy = np.loadtxt(file_path)

    return Vy,w

def choose_choice(choice,c):
    if choice==0:
        a1, a2, a3, b, tau = cases.homo_seats(c)
    elif choice==1:
        a1, a2, a3, b, tau = cases.incre_seats(c)
    elif choice == 2:
        a1, a2, a3, b, tau = cases.aw_seats(c)
    return a1, a2, a3, b, tau

def findseats(rand,p1,pj):
    #given a random number, find which seat is chosen for product 2 or 3
    sum=p1
    for i in range(len(pj)):
        sum=sum+pj[i]
        if sum>=rand:
            return i
    raise TypeError("Math Wrong")

def Simulation(random):
    # Initialize CVXPY variables for probabilities
    p0 = cp.Variable(1, name="p0")
    p1 = cp.Variable(1, name="p1")
    p2 = cp.Variable(1, name="p2")
    p3 = cp.Variable(1, name="p3")
    p2j = cp.Variable(c, name="p2j")
    p3j = cp.Variable(c, name="p3j")

    #initialize the availability parameters
    x = np.full(c,1)
    y = c
    Revenue=0
    Vy,w=read(c,choice)

    for t in range(T,0,-1):
        if y==0:
        #no sale, sale over
            break

        if y==1:
        #no sale of product 3
            objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 + p1 * Vy[t - 1, y - 1] \
                           + cp.sum(
                [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (
                            w[j] - Vy[t - 1, y - 1]) *
                 p2j[j] for j in range(c)]) \
                           + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                           +p0*Vy[t-1,y]+ sum(w[j] * x[j] for j in range(c))
        else:
        #sale for all 3 products
            objective_cp = 1 / (b[0] * tau[0]) * (-cp.kl_div(p1, p0) + p0 - p1) + a1 / b[0] * p1 + p1*Vy[t-1,y-1] \
                           + cp.sum(
                [1 / b[1] * (-cp.kl_div(p2j[j], p0) + p0 - p2j[j]) + a2[j] / b[1] * p2j[j] - (w[j]-Vy[t-1,y-1]) *
                 p2j[j] for j in range(c)]) \
                           + 1 / b[1] * (1 - tau[1]) / tau[1] * (-cp.kl_div(p2, p0) - p2 + p0) \
                           + cp.sum([1 / b[2] * (-cp.kl_div(p3j[j], p0) + p0 - p3j[j]) + a3[j] / b[2] * p3j[j] -(w[j]+w[j - (-1) ** (j + 1)]-Vy[t-1,y-2]) * p3j[j] for j in range(c)]) \
                           + 1 / b[2] * (1 - tau[2]) / tau[2] * (-cp.kl_div(p3, p0) - p3 + p0)+p0*Vy[t-1,y]+sum(w[j]*x[j] for j in range(c))
        objective = cp.Maximize(objective_cp)
        #adding constraints
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
        #solving the problem
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

        '''
        #for debugging
        if p0[0].value is None or p1[0].value is None or p2[0].value is None or p3[0].value is None:
            raise ValueError("None again, tried=", p0[0].value ,p1[0].value,p2[0].value,p3[0].value)
        '''
        #obtain the results
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

        #decide on the seat based on random number
        if random[t - 1] > p10 + p20 + p30:
            #no purchase
            continue
        elif random[t - 1] < p10:
            #buy product 1
            Revenue = Revenue + price1
            y = y - 1
        elif random[t - 1] < p10 + p20:
            #buy product 2
            seat = findseats(random[t - 1], p10, p2j0)
            if x[seat] == 0:
                raise ValueError("the seat for product 2 is unavailable for purchase")
            else:
                Revenue = Revenue + price2[seat]
                y = y - 1
                x[seat] = 0
        else:
            #buy product 3
            seat = findseats(random[t - 1], p10 + p20, p3j0)
            if x[seat] == 0:
                raise ValueError("the seat for product 3 is unavailable for purchase")
            else:
                if x[seat - (-1) ** (seat + 1)] == 0:
                    raise ValueError("the adj seat for product 3 is unavailable for purchase")
                else:
                    y = y - 2
                    Revenue = Revenue + price3[seat]
                    x[seat] = 0
                    x[seat - (-1) ** (seat + 1)] = 0
    return Revenue

num_sim=100
mean_size=np.zeros((10,3))
var_size=np.zeros((10,3))
for j in range(1,11):
    #capacity from 8 to 80
    c=j*8
    T=c*2
    folder_path = "RandomNumbers/"+str(c)
    file_name = "rand.txt"
    file_path = f"{folder_path}/{file_name}"
    random = np.loadtxt(file_path)
    for choice in range(3):
        a1, a2, a3, b, tau=choose_choice(choice,c)
        results=[0]*num_sim
        #simulation process
        for i in range(num_sim,2*num_sim):
            rand=random[:,i]
            results[i-num_sim]=Simulation(rand)

        directory = "results/SBD"
        os.makedirs(directory, exist_ok=True)
        full_path = directory + "/" + str(choice) + str(c) + "CRN.txt"
        np.savetxt(full_path, results)

        mean_size[j-1,choice]=np.mean(results)
        var_size[j-1,choice]=np.var(results)
        print("c,choice:",c,choice)
        print("mean,var:",mean_size[j-1,choice],var_size[j-1,choice])