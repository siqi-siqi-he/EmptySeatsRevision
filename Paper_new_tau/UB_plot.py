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
    file_name = "UB_DBD"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD=np.loadtxt(file_path)
    DBD=fill(DBD[1:])

    folder_path = "results"
    file_name = "UB_DBD_ke"+str(choice)+"capacity"+".txt"
    file_path = f"{folder_path}/{file_name}"
    try:
        DBD_ke=np.loadtxt(file_path)
        DBD_ke=fill(DBD_ke[1:])
        DBD_ke=np.minimum(DBD_ke,DBD)
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
    return ADP, SBD, SBD_ke, DLP, DBD, DBD_ke

width=4
marksize=10
plt.rcParams['font.size'] = 20
fs = 20
figs = (15,10)
#homo
plt.figure(figsize=figs)
ADP, SBD, SBD_ke, DLP, DBD, DBD_ke=read(1)

plt.plot(x,DLP/SBD,label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, SBD/SBD, label='UB DPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, DBD/SBD, label='UB TD-DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP/SBD, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/SBD, label='UB DPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke/SBD, label='UB TD-DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)




plt.xticks(x)
plt.ylim(0.98, 1.12)
plt.xlabel('Bus size',fontsize=fs)
plt.ylabel('UB relative to UB DPD',fontsize=fs)
plt.legend(ncol=2)


#aw
plt.figure(figsize=figs)
ADP, SBD, SBD_ke, DLP, DBD, DBD_ke=read(2)

plt.plot(x,DLP/SBD,label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, SBD/SBD, label='UB DPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, DBD/SBD, label='UB TD-DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP/SBD, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/SBD, label='UB DPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke/SBD, label='UB TD-DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)




plt.xticks(x)
plt.ylim(0.98, 1.12)
plt.xlabel('Bus size',fontsize=fs)
plt.ylabel('UB relative to UB DPD',fontsize=fs)
plt.legend(ncol=2)


#het
plt.figure(figsize=figs)
ADP, SBD, SBD_ke, DLP, DBD, DBD_ke=read(3)

plt.plot(x,DLP/SBD,label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, SBD/SBD, label='UB DPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, DBD/SBD, label='UB TD-DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, ADP/SBD, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke/SBD, label='UB DPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke/SBD, label='UB TD-DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)



plt.xticks(x)
plt.ylim(0.98, 1.12)
plt.xlabel('Bus size',fontsize=fs)
plt.ylabel('UB relative to UB DPD',fontsize=fs)
plt.legend(ncol=2)

plt.show()
