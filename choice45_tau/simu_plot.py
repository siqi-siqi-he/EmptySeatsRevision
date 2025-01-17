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
    file_name = "UB_SBD"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    UB=np.loadtxt(file_path)
    UB=fill(UB[1:])

    folder_path = "simu_results"
    file_name = "mean_ADP_simu_choice_"+str(choice)+"_CRN.txt"
    file_path = f"{folder_path}/{file_name}"
    ADP=np.loadtxt(file_path)
    ADP=fill(ADP)

    folder_path = "simu_results"
    file_name = "mean_DBD_NL_ke_simu_choice_"+str(choice)+"_CRN.txt"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD_ke=np.loadtxt(file_path)
        DBD_ke=fill(DBD_ke)
    except:
        DBD_ke=fill([])
    
    folder_path = "simu_results"
    file_name = "mean_DBD_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD=np.loadtxt(file_path)
        DBD=fill(DBD)
    except:
        DBD=fill([])


    folder_path = "simu_results"
    file_name = "mean_SBD_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD=np.loadtxt(file_path)
    SBD=fill(SBD)
    #checked

    folder_path = "simu_results"
    file_name = "mean_SBD_NL_ke_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke=np.loadtxt(file_path)
    SBD_ke=fill(SBD_ke)
    
    folder_path = "simu_results"
    file_name = "mean_DLP_simu_choice_"+str(choice)+"_CRN"+".txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=fill(DLP)

    '''
    folder_path = "EmptySeatsStudy (experiments)/sbADP"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    sbADP_simu=np.loadtxt(file_path)
    '''
    return UB, SBD, SBD_ke, DLP, ADP, DBD, DBD_ke

width=2
marksize=10
plt.rcParams['font.size'] = 14
#homo1

#pricing policies
plt.figure(figsize=(10, 6))
UB, SBD, SBD_ke, DLP, ADP, DBD, DBD_ke=read(1)

print(UB)
print(SBD)
plt.plot(x, SBD/UB, label='Policy DPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD/UB, label='Policy TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/UB, label='Policy DPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke/UB, label='Policy TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP/UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP/UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=100)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)

plt.xticks(x)
plt.ylim(0.85, 1)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Homogeneous Seats (scaled by UB DPD)')
plt.legend(ncol=2)

#aw2

plt.figure(figsize=(10, 6))

UB, SBD, SBD_ke, DLP, ADP, DBD, DBD_ke=read(2)
plt.plot(x, SBD/UB, label='Policy DPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD/UB, label='Policy TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/UB, label='Policy DPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke/UB, label='Policy TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP/UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP/UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)

plt.xticks(x)
plt.ylim(0.85, 1)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Window Aisle Seats (scaled by UB DPD)')
plt.legend(ncol=2)

#hete3
plt.figure(figsize=(10, 6))
UB, SBD, SBD_ke, DLP, ADP, DBD, DBD_ke=read(3)

plt.plot(x, SBD/UB, label='Policy DPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD/UB, label='Policy TDPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/UB, label='Policy DPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke/UB, label='Policy TDPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP/UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP/UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)

plt.xticks(x)
plt.ylim(0.85, 1)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Heterogeneous Seats (scaled by UB DPD)')
plt.legend(ncol=2)

plt.show()