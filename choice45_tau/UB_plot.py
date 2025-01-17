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
    file_name = "UB_ADP"+str(choice)+"capacity"+".csv"
    file_path = f"{folder_path}/{file_name}"
    ADP=pd.read_csv(file_path)
    ADP=ADP.sort_values('i')
    ADP=fill(ADP['Result'].to_numpy())

    folder_path = "results"
    file_name = "UB_DBD"+str(choice)+"capacity"+".csv"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD=pd.read_csv(file_path)
        DBD=DBD.sort_values('i')
        DBD=fill(DBD['Result'].to_numpy())
    except:
        DBD=fill([])

    folder_path = "results"
    file_name = "UB_DBD_ke"+str(choice)+"capacity"+".csv"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD_ke=pd.read_csv(file_path)
        DBD_ke=DBD_ke.sort_values('i')
        DBD_ke=fill(DBD_ke['Result'].to_numpy())
        DBD_ke=np.minimum(DBD,DBD_ke)
    except:
        DBD_ke=fill([])

    folder_path = "results"
    file_name = "UB_SBD"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD=np.loadtxt(file_path)
    SBD=fill(SBD[1:])
    #checked

    folder_path = "results"
    file_name = "UB_SBD_ke"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke=np.loadtxt(file_path)
    SBD_ke=fill(SBD_ke[1:])
    SBD_ke=np.minimum(SBD_ke,SBD)
    
    folder_path = "results"
    file_name = "UB_DLPnorm"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=fill(DLP[1:])

    '''
    folder_path = "EmptySeatsStudy (experiments)/sbADP"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    sbADP_simu=np.loadtxt(file_path)
    '''

    folder_path = "simu_results"
    file_name = "mean_ADP_simu_choice_"+str(choice)+"_CRN.txt"
    file_path = f"{folder_path}/{file_name}"
    ADP_simu=np.loadtxt(file_path)
    ADP_simu=fill(ADP_simu)

    folder_path = "simu_results"
    file_name = "mean_DBD_NL_ke_simu_choice_"+str(choice)+"_CRN.txt"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD_ke_simu=np.loadtxt(file_path)
        DBD_ke_simu=fill(DBD_ke_simu)
    except:
        DBD_ke_simu=fill([])
    
    folder_path = "simu_results"
    file_name = "mean_DBD_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD_simu=np.loadtxt(file_path)
        DBD_simu=fill(DBD_simu)
    except:
        DBD_simu=fill([])

    folder_path = "simu_results"
    file_name = "mean_SBD_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_simu=np.loadtxt(file_path)
    SBD_simu=fill(SBD_simu)
    #checked

    folder_path = "simu_results"
    file_name = "mean_SBD_ke_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke_simu=np.loadtxt(file_path)
    SBD_ke_simu=fill(SBD_ke_simu)
    
    folder_path = "simu_results"
    file_name = "mean_DLP_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    DLP_simu=np.loadtxt(file_path)
    DLP_simu=fill(DLP_simu)

    return ADP, SBD, SBD_ke, DLP, DBD, DBD_ke, ADP_simu, DBD_ke_simu, DBD_simu, SBD_simu, SBD_ke_simu, DLP_simu

width=2
marksize=10
plt.rcParams['font.size'] = 14
#homo
plt.figure(figsize=(10, 6))
ADP, SBD, SBD_ke, DLP, DBD, DBD_ke, ADP_simu, DBD_ke_simu, DBD_simu, SBD_simu, SBD_ke_simu, DLP_simu=read(5)

plt.plot(x, SBD/SBD, label='UB DPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD/SBD, label='UB TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/SBD, label='UB DPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, DBD_ke/SBD, label='UB TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP/SBD, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x,DLP/SBD,label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
plt.plot(x, SBD_simu/SBD, label='Policy DPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD, label='Policy TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD, label='Policy DPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD, label='Policy TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP_simu/SBD, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)


plt.xticks(x)
plt.ylim(0.85, 1.3)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Expected Revenue Upper Bound',fontsize=14)
plt.legend(ncol=2)

plt.show()