import matplotlib.pyplot as plt
import numpy as np
x0=[8,16,24,32,40,48,56,64,72,80]

#this is the run time on the server, after the submission
def regular(a):
    a = [a[i] / (i + 1) / 8 for i in range(len(a))]
    return a
def fill(a):
    n=len(x0)-len(a)
    if n <0:
        raise Exception("x too short")
    a=np.append(a,np.full(n,np.nan))
    return a

def check(a):
    a=regular(a)
    a=fill(a)
    return a

plt.rcParams['font.size'] = 14
width=2
marksize=10
#time plot
#homo 0
y0=[0,10000,20000,30000,40000,50000]

#DBD_time=fill(np.array([22, 1509,5603,3202,8210])) this is the time to run DBD itself
DBD_time=fill(np.array([22+408, 1509+704,5603+1500,3202+4522,8210+13155]))
#DBD_ke_time=fill(np.array([194,2011,15476]))  this is the time to run DBD_ke itself
DBD_ke_time=fill(np.array([194+22+408,2011+1509+704,15476+5603+1500]))
SBD_time=fill(np.array([11,70,215,485,920,1543,2418,3581,4892,6708]))
#SBD_ke_time=fill(np.array([22,143,442,993,1916,3155,4958,7345,10393,14582])) this is the time to run SBD_ke itself
SBD_ke_time=fill(np.array([33,213,657,1478,2836,3155+1543,4958+2418,7345+3581,10393+4892,14582+6708]))
aALP_time=fill(np.array([408,704,1500,4522,13155,35077]))
sbADP_time=fill(np.array([52,105,219,462,1106,2196,4163,7478,13430,21077]))
DLPflex=fill(np.array([0.119138002,0.154306412,0.206588268,0.271484137,0.324224234,0.36230588,0.444319725,0.48730135,0.482412815,0.582036018]))


plt.figure(figsize=(10, 6))

plt.plot(x0, SBD_time, label='SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_time, label='DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, SBD_ke_time, label='SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_ke_time, label='DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, aALP_time, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DLPflex, label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
plt.plot(x0, sbADP_time, label='sbADP (M=100)',color='salmon', marker='<', linestyle='-', linewidth=width, markersize=marksize)

plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Run Time',fontsize=14)
#plt.title('Run Time (seconds) for Homogeneous Seats')
plt.legend()
plt.xticks(x0)
plt.yticks(y0)
#aw2
#DBD_time=fill(np.array([219, 2492,12349,13735,24191])) to compute DBD itself
DBD_time=fill(np.array([219+146, 2492+911,12349+1643,13735+4114,24191+13381]))
#DBD_ke_time=fill(np.array([381,3664,16503])) to compute DBD_ke itself
DBD_ke_time=fill(np.array([381+219+146,3664+2492+911,16503+12349+1643]))
SBD_time=fill(np.array([11,70,213,483,914,1540,2412,3464,4928,6660]))
#SBD_ke_time=fill(np.array([22,144,443,989,1891,3156,4971,7345,10524,14083])) SBD_ke itself
SBD_ke_time=fill(np.array([22+11,144+70,443+213,989+483,1891+914,3156+1540,4971+2412,7345+3464,10524+4928,14083+6660]))
aALP_time=fill(np.array([146,911,1643,4114,13381,34267]))
sbADP_time_a=fill(np.array([264,528, 1127, 2500,6023,12373, 23225, 42083]))
sbADP_time=fill(np.array([0,264,528, 1127, 2500,6023,12373, 23225, 42083,73170]))
DLPflex=fill(np.array([0.081054211,0.127912998,0.201182604,0.28028059,0.282226801,0.370129824,0.384749174,0.462885141,0.486329556,0.589525461]))

plt.figure(figsize=(10, 6))


plt.plot(x0, SBD_time, label='SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_time, label='DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, SBD_ke_time, label='SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_ke_time, label='DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, aALP_time, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DLPflex, label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
plt.plot(x0, sbADP_time, label='sbADP (M=500)',color='salmon', marker='<', linestyle='-', linewidth=width, markersize=marksize)

plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Run Time',fontsize=14)
#plt.title('Run Time (seconds) for Window Aisle Seats')
plt.legend()
plt.xticks(x0)
plt.yticks(y0)
#hete1
#DBD_time=fill(np.array([203, 2321,12356,15104,36293])) DBD itself
DBD_time=fill(np.array([203+190, 2321+484,12356+1326,15104+5235,36293+14810]))
#DBD_ke_time=fill(np.array([311,3385])) DBD_ke itself
DBD_ke_time=fill(np.array([311+203+190,3385+2321+484]))
SBD_time=fill(np.array([11,70,214,485,915,1535,2418,3573,4911,6676]))
#SBD_ke_time=fill(np.array([23,144,442,986,1914,3160,4977,7330,10579,14328])) SBD_ke itself
SBD_ke_time=fill(np.array([23+11,144+70,442+214,986+485,1914+915,3160+1535,4977+2418,7330+3575,10579+4911,14328+6676]))
aALP_time=fill(np.array([190,484,1326,5235,14810,39616]))
sbADP_time_h=fill(np.array([198,582,1236,2674,6583,13063,23957]))
sbADP_time=fill(np.array([0,198,582,1236,2674,6583,13063,23957,26668,72921]))
DLPflex=fill(np.array([0.102527142,0.146480799,0.18261385,0.247063637,0.326700687,0.334962845,0.403319359,0.459959984,0.486333847,0.57419467]))

plt.figure(figsize=(10, 6))

plt.plot(x0, SBD_time, label='SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_time, label='DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, SBD_ke_time, label='SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_ke_time, label='DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, aALP_time, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DLPflex, label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
plt.plot(x0, sbADP_time, label='sbADP (M=500)',color='salmon', marker='<', linestyle='-', linewidth=width, markersize=marksize)


plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Run Time',fontsize=14)
#plt.title('Run Time (seconds) for Heterogeneous Seats')
plt.legend()
plt.xticks(x0)
plt.yticks(y0)


plt.show()