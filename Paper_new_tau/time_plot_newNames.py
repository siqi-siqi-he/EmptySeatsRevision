import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
x=[8,16,24,32,40,48,56,64,72,80]

def regular(a):
    a = [a[i] / (i + 1) / 8 for i in range(len(a))]
    return a
def fill(a):
    n=len(x)-len(a)
    if n <0:
        raise Exception("x too short")
    a=np.append(a,np.full(n,np.nan))
    return a

def check(a):
    a=regular(a)
    a=fill(a)
    return a

def read(choice):

    folder_path = "results"
    file_name = "runtime_ADP"+str(choice)+"capacity"+".csv"
    file_path = f"{folder_path}/{file_name}"
    ADP=pd.read_csv(file_path)
    ADP=ADP.sort_values('i')
    ADP=fill(ADP['Result'].to_numpy())
    '''
    folder_path = "results/DBD"
    file_name = "ke_simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_ke_simu=np.loadtxt(file_path)
    '''
    folder_path = "results"
    file_name = "runtime_DBD"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD=np.loadtxt(file_path)
        DBD=DBD+ADP[:len(DBD)]
        DBD=fill(DBD[1:])
    except:
        DBD=fill([])

    folder_path = "results"
    file_name = "runtime_DLPnorm"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=fill(DLP[1:])

    folder_path = "results"
    file_name = "runtime_SBD"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD=np.loadtxt(file_path)
    SBD=fill(SBD[1:])
    SBD=SBD+DLP
    #checked

    folder_path = "results"
    file_name = "runtime_SBD_ke"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke=np.loadtxt(file_path)
    SBD_ke=fill(SBD_ke[1:])
    SBD_ke=SBD_ke+SBD

    '''
    folder_path = "EmptySeatsStudy (experiments)/sbADP"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    sbADP_simu=np.loadtxt(file_path)
    '''
    return ADP, SBD, SBD_ke, DLP, DBD

width=2
marksize=10
plt.rcParams['font.size'] = 14
#homo
plt.figure(figsize=(10, 6))
ADP, SBD, SBD_ke, DLP, DBD=read(1)

plt.plot(x, SBD, label='DPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD, label='TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke, label='DPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
#plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x,DLP,label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)

plt.xticks(x)
plt.ylim(-5000,51000)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Runtime',fontsize=14)
plt.legend(ncol=2)


#AW
plt.figure(figsize=(10, 6))
ADP, SBD, SBD_ke, DLP, DBD=read(2)

plt.plot(x, SBD, label='DPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD, label='TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke, label='DPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
#plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x,DLP,label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)

plt.xticks(x)
plt.ylim(-5000,51000)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Runtime',fontsize=14)
plt.legend(ncol=2)


#HETE
plt.figure(figsize=(10, 6))
ADP, SBD, SBD_ke, DLP, DBD=read(3)

plt.plot(x, SBD, label='DPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD, label='TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke, label='DPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
#plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x,DLP,label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)

plt.xticks(x)
plt.ylim(-5000,51000)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Runtime',fontsize=14)
plt.legend(ncol=2)

plt.show()
