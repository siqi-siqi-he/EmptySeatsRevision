import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

x=np.arange(-5,5.5,0.5)



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

def read_DLP():
    folder_path = "DLP_NL"
    file_name = "simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

def read_DLP_wo3():
    folder_path = "wo3_DLP_NL"
    file_name = "simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

def read_SBD_ke():
    folder_path = "SBD_NL_ke"
    file_name = "ke_simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "ke_simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

def read_SBD_ke_wo3():
    folder_path = "wo3_SBD_NL_ke"
    file_name = "ke_simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "ke_simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

def read_SBD():
    folder_path = "SBD_NL"
    file_name = "simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

def read_SBD_wo3():
    folder_path = "wo3_SBD_NL"
    file_name = "simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

def read_ADP():
    folder_path = "ADP"
    file_name = "simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

def read_ADP_wo3():
    folder_path = "wo3_ADP"
    file_name = "simu_means_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    DLP=np.loadtxt(file_path)
    DLP=DLP[1][:]
    print(DLP)

    file_name = "simu_SLF_choice.txt"
    file_path = f"{folder_path}/{file_name}"
    SL=np.loadtxt(file_path)
    SL=SL[1][:]
    return DLP[0::2]*100, SL[0::2]*100

width=2
marksize=10
#homo
plt.figure(figsize=(15, 9))
DLP, DLP_SL = read_DLP()
DLP_wo3, DLP_SL_wo3 = read_DLP_wo3()
SBD_ke, SBD_ke_SL = read_SBD_ke()
SBD_ke_wo3, SBD_ke_SL_wo3 = read_SBD_ke_wo3()
SBD, SBD_SL = read_SBD()
SBD_wo3, SBD_SL_wo3 = read_SBD_wo3()
ADP, ADP_SL = read_ADP()
ADP_wo3, ADP_SL_wo3 = read_ADP_wo3()

x = np.arange(-5,6)
xticksvector = [r'Strong, $a_0=-5$', '', '', '', '', '', '', '', '', '', r'Weak, $a_0=5$']
plt.rcParams['font.size'] = 14

plt.plot(x, DLP/DLP_wo3*100, linestyle='-', marker='o', label="Relative Expected Revenue (Policy DPP)", color='violet')
plt.plot(x, SBD_ke/SBD_ke_wo3*100, linestyle='-', marker='o', label="Relative Expected Revenue (Policy DPD-Benchmark)", color='red')
plt.plot(x, SBD/SBD_wo3*100, linestyle='-', marker='o', label="Relative Expected Revenue (Policy DPD)", color='green')
plt.plot(x, ADP/ADP_wo3*100, linestyle='-', marker='o', label="Relative Expected Revenue (Policy AFF)", color='blue')
plt.plot(x, DLP_SL/8, linestyle='--', marker='o', label="SLF (Policy DPP) with extra seats", color='violet')
plt.plot(x, DLP_SL_wo3/8, linestyle='--', marker='s', label="SLF (Policy DPP) without extra seats", color='violet')
# plt.plot(x, SBD_ke_SL/8, linestyle='--', marker='o', label="SLF (Policy DPD-Benchmark) with extra seats", color='red')
# plt.plot(x, SBD_ke_SL_wo3/8, linestyle='--', marker='s', label="SLF (Policy DPD-Benchmark) without extra seats", color='red')
plt.plot(x, SBD_SL/8, linestyle='--', marker='o', label="SLF (Policy DPD) with extra seats", color='green')
plt.plot(x, SBD_SL_wo3/8, linestyle='--', marker='s', label="SLF (Policy DPD) without extra seats", color='green')
# plt.plot(x, ADP_SL/8, linestyle='--', marker='o', label="SLF (Policy AFF) with extra seats", color='blue')
# plt.plot(x, ADP_SL_wo3/8, linestyle='--', marker='s', label="SLF (Policy AFF) without extra seats", color='blue')

plt.xlabel('Demand')
plt.ylabel('Percentage')
# plt.title(f'Relative Expected Revenue by selling Extra Seats and SLF over Demand Ratio (HOMOG) for a3={a3:.1f}')
plt.xticks(x, xticksvector)
plt.legend(fontsize=14,loc='best',ncol=2)

plt.xlim(min(x) - 0.3, max(x) + 0.3)
plt.ylim(0, 240)


plt.show()
