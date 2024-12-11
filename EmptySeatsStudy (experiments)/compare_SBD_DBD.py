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

def read_SBD(c,choice):
    folder_path = "EmptySeatsStudy (experiments)/DP"
    file_name = "DLPnorm" + str(choice) + "capacity" + str(c) + "step" + str(step) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    w = np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/SBD/"
    file_name = "SBD_NL_DPtable" + str(choice) + "capacity" + str(c) + "step"+str(step)+".txt"
    file_path = f"{folder_path}/{file_name}"
    Vy = np.loadtxt(file_path)

    return Vy,w

def read_SBD_ke(c, choice):
    V = np.full((c, c * 2 + 1, 2), 0.0)

    folder_path = "EmptySeatsStudy (experiments)/SBD"
    for i in range(c):
        file_name = "SBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo" + str(i) + "step" + str(
            step) + ".txt"
        file_path = f"{folder_path}/{file_name}"
        V[i, :, :] = np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/SBD"
    file_name = "SBD_NL_DPtable" + str(choice) + "capacity" + str(c) + "step" + str(step) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    Vy = np.loadtxt(file_path)

    return V, Vy

def read_DBD_ke(c,choice):

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

def read_DBD(c,choice):
    folder_path = "EmptySeatsStudy (experiments)/ADP"
    file_name = "v_value_choice" + str(choice) + "capacity" + str(c) + "step" + str(step) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    v = np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/DBD"
    file_name = "DBD_NL_DPtable" + str(choice) + "capacity" + str(c) + "step" + str(step) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    Vy = np.loadtxt(file_path)

    return Vy,v

def findseats(rand,p1,pj):
    sum=p1
    for i in range(len(pj)):
        sum=sum+pj[i]
        if sum>=rand:
            return i
    raise TypeError("Math Wrong")

def Simulation(i):

    V_D, Vy_D=read_DBD_ke(c,choice)
    Vy_S, w_S=read_SBD(c,choice)
    Vy_D,v_D=read_DBD(c,choice)
    V_S, Vy_S=read_SBD_ke(c, choice)

    for t in range(T,0,-1):
        print("SBD: ", w_S, Vy_S[t,:])
        print("SBD_ke when full: ", V_S[:,t,1], Vy_S[t,:])
        print("SBD_ke when full: ", V_S[:,t,0], Vy_S[t,:])
        print("DBD_ke when full: ", V_D[:,t,1], Vy_D[t,:])
        print("DBD_ke when full: ", V_D[:,t,0], Vy_D[t,:])
        print("DBD: ", v_D[t,:],Vy_D[t,:])



def run_Simulator():

    for ii in range(100):
        Simulation(ii)




c = 8
choice=0
a1, a2, a3, b, tau = cases.homo_seats(c)
means=[0]*51
vars=[0]*51
p3=[0]*51
T=c*2
for step in range(51):
    a3=[0.1*step for i in range(len(a3))]
    run_Simulator()

