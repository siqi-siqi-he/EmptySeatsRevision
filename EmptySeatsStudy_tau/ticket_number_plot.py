import numpy as np
import matplotlib.pyplot as plt

choice=0
c=8

def read():
    # folder_path = "ADP"
    # file_name = "ADP_NL_results" + str(choice) + "capacity" + str(c) + ".txt"
    # file_path = f"{folder_path}/{file_name}"
    # ADP = np.loadtxt(file_path)

    folder_path = "ADP"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    ADP_simu=np.loadtxt(file_path)

    # folder_path = "DBD"
    # file_name = "DBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    # file_path = f"{folder_path}/{file_name}"
    # DBD_ke=np.loadtxt(file_path)

    folder_path = "DBD_ke"
    file_name = "ke_simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_ke_simu=np.loadtxt(file_path)

    folder_path = "DBD"
    file_name = "DBD_NL_results0capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD=np.loadtxt(file_path)

    # for i in range(len(DBD_ke)):
    #     if DBD_ke[i]>DBD[i]:
    #         DBD_ke[i] = DBD[i]

    folder_path = "DBD"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_simu=np.loadtxt(file_path)

    # folder_path = "SBD"
    # file_name = "SBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    # file_path = f"{folder_path}/{file_name}"
    # SBD_ke=np.loadtxt(file_path)

    folder_path = "SBD_ke"
    file_name = "ke_simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke_simu=np.loadtxt(file_path)

    folder_path = "SBD"
    file_name = "SBD_NL_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD=np.loadtxt(file_path)

    # for i in range(len(SBD_ke)):
    #     if SBD_ke[i]>SBD[i]:
    #         SBD_ke[i] = SBD[i]

    folder_path = "SBD"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_simu=np.loadtxt(file_path)

    folder_path = "DP"
    file_name = "simu_p3_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DP_simu=np.loadtxt(file_path)

    return SBD_simu, SBD, DBD, SBD_ke_simu, DBD_simu, DBD_ke_simu, ADP_simu, DP_simu

width=4
marksize=10
plt.rcParams['font.size'] = 20
fs = 20
figs = (15,10)
SBD_simu, SBD, DBD, SBD_ke_simu, DBD_simu, DBD_ke_simu, ADP_simu, DP_simu=read()
back_simu=[1.2,1.24,1.26,1.28,1.3,1.31,1.33,1.35,1.36,1.34,1.34,1.33,1.35,1.35,1.36,1.35,1.41,1.42,1.45,1.45,1.47,1.47,1.52,1.54,1.54,1.57,1.6,1.61,1.63,1.66,1.69,1.7,1.74,1.76,1.77,1.79,1.8,1.84,1.85,1.87,1.88,1.87,1.91,1.91,1.93,1.99,2,2.01,2.03,2.06,2.09]
sbADP_simu=[1.21,1.19,1.22,1.21,1.23,1.26,1.28,1.3,1.3,1.3,1.31,1.32,1.33,1.33,1.36,1.36,1.38,1.37,1.37,1.37,1.4,1.4,1.42,1.47,1.51,1.55,1.56,1.57,1.58,1.61,1.65,1.68,1.7,1.7,1.7,1.71,1.72,1.74,1.75,1.8,1.83,1.83,1.88,1.91,1.92,1.95,1.99,2.05,2.07,2.08,2.1]

x0=[0,1,2,3,4,5]
plt.figure(figsize=figs)
x=[i*0.1 for i in range(51)]
#plt.plot(x, SBD, label='UB SDPD (SDPD-Benchmark)',color='forestgreen', marker='x',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
#plt.plot(x, DBD, label='UB DPD (DPD-Benchmark)',color='lightsteelblue', marker='s', markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, SBD_ke, label='UB ',color='red', marker='o', markerfacecolor='none', linestyle='-.', linewidth=width, markersize=marksize)
#plt.plot(x, DBD, label='UB DPD-Benchmark',color='yellow', marker='s', markerfacecolor='none', linestyle=' ', linewidth=width, markersize=marksize)
#plt.plot(x, ADP, label='UB AFF',color='black', marker='x', linestyle='-', linewidth=1.5, markersize=5)
plt.plot(x, DP_simu, label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, SBD_simu, label='Policy DPD',color='forestgreen', marker='o', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, DBD_simu, label='Policy TD-DPD',color='lightsteelblue', marker='s', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, back_simu, label='Policy BE',color='cyan', marker='^', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, ADP_simu, label='Policy AFF',color='black', marker='x', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu, label='Policy DPD-Benchmark',color='red', marker='o', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu, label='Policy TD-DPD-Benchmark',color='yellow', marker='s', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, sbADP_simu, label='Policy sbADP (M=100)',color='salmon', marker='+', linestyle=':', linewidth=width, markersize=marksize)


x0=[i*0.5 for i in range(11)]

plt.ylim(0.8, 2.25)
plt.xticks(x0)
plt.xlabel('Quality index of product type 3',fontsize=fs)
plt.ylabel('Number of tickets sold for product type 3 ',fontsize=fs)
#plt.title('Homogeneous Seats')
plt.legend(ncol=2)

plt.show()