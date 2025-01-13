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


width=1.5
marksize=5
SBD_simu, SBD, DBD, SBD_ke_simu, DBD_simu, DBD_ke_simu, ADP_simu, DP_simu=read()
#back_simu=[0.44,0.44,0.49,0.53,0.55,0.61,0.66,0.72,0.72,0.72,0.76,0.78,0.78,0.8,0.84,0.85,0.83,0.84,0.9,0.92,0.95,0.97,0.99,1,1.03,1.08,1.09,1.17,1.16,1.2,1.22,1.25,1.3,1.32,1.35,1.36,1.35,1.39,1.41,1.44,1.5,1.51,1.54,1.56,1.59,1.64,1.66,1.69,1.72,1.78,1.82]
#sbADP_simu=[0.48,0.5,0.55,0.54,0.63,0.59,0.71,0.74,0.76,0.76,0.74,0.78,0.81,0.8,0.82,0.85,0.87,0.88,0.93,0.94,0.94,1,1.02,1.03,1.03,1.08,1.13,1.14,1.16,1.19,1.19,1.19,1.23,1.29,1.31,1.36,1.43,1.43,1.43,1.45,1.49,1.54,1.5,1.54,1.57,1.64,1.64,1.69,1.7,1.73,1.77]

x0=[0,1,2,3,4,5]
plt.figure(figsize=(10, 6))
plt.rcParams['font.size'] = 14
x=[i*0.1 for i in range(51)]
#plt.plot(x, SBD, label='UB SDPD (SDPD-Benchmark)',color='forestgreen', marker='x',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
#plt.plot(x, DBD, label='UB DPD (DPD-Benchmark)',color='lightsteelblue', marker='s', markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, SBD_ke, label='UB ',color='red', marker='o', markerfacecolor='none', linestyle='-.', linewidth=width, markersize=marksize)
#plt.plot(x, DBD, label='UB DPD-Benchmark',color='yellow', marker='s', markerfacecolor='none', linestyle=' ', linewidth=width, markersize=marksize)
#plt.plot(x, ADP, label='UB AFF',color='black', marker='x', linestyle='-', linewidth=1.5, markersize=5)

plt.plot(x, SBD_simu, label='Policy SDPD',color='forestgreen', marker='o', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, SBD_ke_simu, label='Policy SDPD-Benchmark',color='red', marker='o', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, DBD_simu, label='Policy DPD',color='lightsteelblue', marker='s', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, DBD_ke_simu, label='Policy DPD-Benchmark',color='yellow', marker='s', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, ADP_simu, label='Policy AFF',color='black', marker='x', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, DP_simu, label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu, label='Policy sbADP (M=100)',color='salmon', marker='+', linestyle=':', linewidth=1.5, markersize=5)
#plt.plot(x, back_simu, label='Policy BE',color='cyan', marker='^', linestyle=':', linewidth=1.5, markersize=5)

x0=[i*0.5 for i in range(11)]

plt.ylim(0.75, 2.25)
plt.xticks(x0)
plt.xlabel('Quality Index of Product Type 3')
plt.ylabel('Number of Tickets Sold for Product Type 3 ')
#plt.title('Homogeneous Seats')
plt.legend(ncol=2)

plt.show()