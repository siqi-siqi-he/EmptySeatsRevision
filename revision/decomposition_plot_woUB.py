import matplotlib.pyplot as plt
import numpy as np
import os
x=[8,16,24,32,40,48,56,64,72,80]

#This plot is used to generate graphs for the submitted paper
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

def read_DLP(choice):
    DLP=[]
    directory = "revision/results/DLP"
    for c in range(8,81,8):
        full_path = directory + "/"+str(choice)+str(c)+"CRN.txt"
        if not os.path.exists(full_path):
            DLP_simu=fill(DLP)
            return DLP_simu
        DLP_c=np.mean(np.loadtxt(full_path))
        DLP=DLP+[DLP_c]
    return DLP

def read_ALP(choice):
    ALP=[]
    directory = "revision/results/ADP"
    for c in range(8,81,8):
        full_path = directory + "/"+str(choice)+str(c)+"CRN.txt"
        if not os.path.exists(full_path):
            ALP_simu=fill(ALP)
            return ALP_simu
        ALP_c=np.mean(np.loadtxt(full_path))
        ALP=ALP+[ALP_c]
    return ALP

def read_SBD(choice):
    SBD=[]
    directory = "revision/results/SBD"
    for c in range(8,81,8):
        full_path = directory + "/"+str(choice)+str(c)+"CRN.txt"
        if not os.path.exists(full_path):
            SBD_simu=fill(SBD)
            return SBD_simu
        SBD_c=np.mean(np.loadtxt(full_path))
        SBD=SBD+[SBD_c]
    return SBD

def read_SBD_ke(choice):
    SBD_ke=[]
    directory = "revision/results/SBD_ke"
    for c in range(8,81,8):
        full_path = directory + "/"+str(choice)+str(c)+"CRN.txt"
        if not os.path.exists(full_path):
            SBD_ke_simu=fill(SBD_ke)
            return SBD_ke_simu
        SBD_ke_c=np.mean(np.loadtxt(full_path))
        SBD_ke=SBD_ke+[SBD_ke_c]
    return SBD_ke

def read_DBD(choice):
    DBD=[]
    directory = "revision/results/DBD"
    for c in range(8,81,8):
        full_path = directory + "/"+str(choice)+str(c)+"CRN.txt"
        if not os.path.exists(full_path):
            DBD_simu=fill(DBD)
            return DBD_simu
        DBD_c=np.mean(np.loadtxt(full_path))
        DBD=DBD+[DBD_c]
    return DBD

def read_DBD_ke(choice):
    DBD_ke=[]
    directory = "revision/results/DBD_ke"
    for c in range(8,81,8):
        full_path = directory + "/"+str(choice)+str(c)+"CRN.txt"
        if not os.path.exists(full_path):
            DBD_ke_simu=fill(DBD_ke)
            return DBD_ke_simu
        DBD_ke_c=np.mean(np.loadtxt(full_path))
        DBD_ke=DBD_ke+[DBD_ke_c]
    return DBD_ke


width=2
marksize=10
plt.rcParams['font.size'] = 14

#choice 4
SBD_UB=np.array([68.6715819035891,158.999831998647,258.879837610151,365.629534999981,477.975412909706,595.177922039927,716.74938252263,842.334912906324,971.668948389566,1104.54486440026])#done,10 #done
DBD_UB=fill(np.array([68.6795620072274,159.140833508944,258.997945936926,365.739462060416]))
SBDb_UB=np.array([68.6715819035891,158.999831998647,258.879837610151,365.629534999981,477.975412909706,595.177922039927,716.74938252263,842.334912906324,971.668948389566,1104.54486440026])#done,10 #done
#SBDb_UB_woy=np.array([73.4808988822258,165.935911900997,267.284163308995,375.17681565103,488.470401867839,606.490207670472,728.779126008739,855.004026835559,984.916921375032,1118.3232773819])#done,10 #done
DBDb_UB=fill(np.array([]))
aALP=fill(np.array([73.1648187824978,165.481227121692,266.64690363421,374.280216708899,487.151061348171]))#done,7
DLPflex=fill(np.array([73.9884542687187,166.498083154143,267.895439478212,375.836449240972,489.179553435683,607.248130545888,729.586667194739,855.865528932043,985.835661228318,1119.30166244716]))#done,10 #done

SBD_simu=read_SBD(4)
SBD_ke_simu=read_SBD_ke(4)
DBD_simu=read_DBD(4)
DBD_ke_simu=read_DBD_ke(4)
aALP_simu=read_ALP(4)
#sbADP_simu=fill(np.array([64.24429505,143.5083651,207.1127125,303.1144468,488.4390255,702.6241617,598.6003484,805.3462068]))#done,10 #done
DLP_simu=read_DLP(4)

#pricing policies
plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=100)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
'''
plt.plot(x, SBD_UB/SBD_UB, label='UB SDPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD_UB/SBD_UB, label='UB DPD',color='lightsteelblue', marker='s', markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, SBDb_UB/SBD_UB, label='UB SDPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, DBDb_UB/SBD_UB, label='UB DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, aALP/SBD_UB, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DLPflex/SBD_UB, label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
'''
plt.xticks(x)
plt.ylim(0.88, 0.96)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for hom_de_in_seats')
plt.legend(ncol=2)


#choice 5
SBD_UB=fill(np.array([68.2417792891811,157.125405898257,254.399615042252,357.286937677139,464.437725259773,575.044921767565,688.569207247698,804.623485419348,922.916059754153,1043.21923788135]))
DBD_UB=fill(np.array([68.3276274717675,157.179376820253,254.6065222414,357.666471763378]))
SBDb_UB=fill(np.array([68.2417792891811,157.125405898257,254.399615042252,357.286937677139,464.437725259773,575.044921767565,688.569207247698,804.623485419348,922.916059754153,1043.21923788135]))
#SBDb_UB_woy=fill(np.array([73.0423366698426,164.038409933635,262.758985742266,366.75922800273,474.819750250712,586.198952932876,700.395153193808,817.045030978347,935.872814484985,1056.66345147933]))
DBDb_UB=fill(np.array([]))
aALP=fill(np.array([72.7291883494042,163.607569896348,262.224464801718,366.138918220246,474.135437089136]))
DLPflex=fill(np.array([73.5463378022181,164.584524119174,263.334627070838,367.358251532822,475.438554594031,586.835307110368,701.047417623319,817.712069915188,936.553783880678,1057.3559761638]))

SBD_simu=read_SBD(5)
SBD_ke_simu=read_SBD_ke(5)
DBD_simu=read_DBD(5)
DBD_ke_simu=read_DBD_ke(5)
aALP_simu=read_ALP(5)
#sbADP_simu=fill(np.array([49.78824925,145.3247069,213.2878257,313.8770599,484.9166528,620.5563582,693.5898398,578.9596915]))
DLP_simu=read_DLP(5)


plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
'''
plt.plot(x, SBD_UB/SBD_UB, label='UB SDPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD_UB/SBD_UB, label='UB DPD',color='lightsteelblue', marker='s', markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, SBDb_UB/SBD_UB, label='UB SDPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, DBDb_UB/SBD_UB, label='UB DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, aALP/SBD_UB, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DLPflex/SBD_UB, label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
'''
plt.xticks(x)
plt.ylim(0.88, 0.96)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for hom_hom_in_seats')
plt.legend(ncol=2)

plt.show()