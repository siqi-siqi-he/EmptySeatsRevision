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

    folder_path = "simu_results"
    file_name = "UB_ADP"+str(choice)+"capacity"+".csv"
    file_path = f"{folder_path}/{file_name}"
    # ADP=pd.read_csv(file_path)
    # ADP=ADP.sort_values('i')
    # ADP=fill(ADP['Result'].to_numpy())
    '''
    folder_path = "results/DBD"
    file_name = "ke_simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_ke_simu=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/DBD"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_simu=np.loadtxt(file_path)
    '''

    folder_path = "simu_results"
    file_name = "mean_SBD_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD=np.loadtxt(file_path)
    SBD=fill(SBD[1:])
    #checked

    folder_path = "simu_results"
    file_name = "mean_SBD_NL_ke_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke=np.loadtxt(file_path)
    SBD_ke=fill(SBD_ke[1:])
    
    folder_path = "simu_results"
    file_name = "mean_DLP_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=fill(DLP[1:])

    '''
    folder_path = "EmptySeatsStudy (experiments)/sbADP"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    sbADP_simu=np.loadtxt(file_path)
    '''
    return SBD, SBD_ke, DLP

width=2
marksize=10
plt.rcParams['font.size'] = 14
#homo

#pricing policies
plt.figure(figsize=(10, 6))
SBD, SBD_ke, DLP=read(1)

plt.plot(x, SBD/DLP, label='Policy DPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, DBD_simu/DLP, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/DLP, label='Policy DPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP/DLP,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=100)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)

plt.xticks(x)
plt.ylim(0.85, 1.2)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Homogeneous Seats (scaled by Policy DPP)')
plt.legend(ncol=2)

#aw2

plt.figure(figsize=(10, 6))

SBD, SBD_ke, DLP=read(2)
plt.plot(x, SBD/DLP, label='Policy DPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/DLP, label='Policy DPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP/DLP,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)

plt.xticks(x)
plt.ylim(0.85, 1.2)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Window Aisle Seats')
plt.legend(ncol=2)

#hete1
plt.figure(figsize=(10, 6))
SBD, SBD_ke, DLP=read(3)

plt.plot(x, SBD/DLP, label='Policy DPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/DLP, label='Policy DPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP/DLP,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)

plt.xticks(x)
plt.ylim(0.85, 1.2)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Heterogeneous Seats')
plt.legend(ncol=2)

plt.show()