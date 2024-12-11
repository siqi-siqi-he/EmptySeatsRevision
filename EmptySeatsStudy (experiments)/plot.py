import numpy as np
import matplotlib.pyplot as plt

choice=0
c=8

def read():
    folder_path = "EmptySeatsStudy (experiments)/ADP"
    file_name = "ADP_NL_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    ADP = np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/ADP"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    ADP_simu=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/DBD"
    file_name = "DBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_ke=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/DBD"
    file_name = "ke_simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_ke_simu=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/DBD"
    file_name = "DBD_NL_results0capacity" + str(c) + "2.txt"
    file_path = f"{folder_path}/{file_name}"
    DBD=np.loadtxt(file_path)

    for i in range(len(DBD_ke)):
        if DBD_ke[i]>DBD[i]:
            DBD_ke[i] = DBD[i]

    folder_path = "EmptySeatsStudy (experiments)/DBD"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_simu=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/SBD"
    file_name = "SBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/SBD"
    file_name = "ke_simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke_simu=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/SBD"
    file_name = "SBD_NL_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD=np.loadtxt(file_path)

    for i in range(len(SBD_ke)):
        if SBD_ke[i]>SBD[i]:
            SBD_ke[i] = SBD[i]

    folder_path = "EmptySeatsStudy (experiments)/SBD"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_simu=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/DP"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DP_simu=np.loadtxt(file_path)

    folder_path = "EmptySeatsStudy (experiments)/sbADP"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    sbADP_simu=np.loadtxt(file_path)

    return SBD_simu, SBD, SBD_ke, SBD_ke_simu, DBD, DBD_simu, DBD_ke, DBD_ke_simu, ADP, ADP_simu, DP_simu, sbADP_simu


width=1.5
marksize=5
SBD_simu, SBD, SBD_ke, SBD_ke_simu, DBD, DBD_simu, DBD_ke, DBD_ke_simu, ADP, ADP_simu, DP_simu, sbADP_simu=read()

plt.figure(figsize=(10, 6))
plt.rcParams['font.size'] = 14
x=[i*0.1 for i in range(51)]
x00=[0,1,2,3,4,5]
back=[63.48090692,64.47043224,65.87998328,67.75406393,70.0999021,72.89937247]
back_simu=[65.65510235,65.67131896,65.75212986,65.91768584,66.11203923,66.11339464,65.97867942,65.95458266,66.06448423,66.26503826,66.49277512,66.64467551,66.63482759,66.69249495,66.83278497,67.05617255,67.37263299,67.56042949,67.25393189,67.62014215,67.98684891,68.24697849,68.62041555,68.88793616,68.80859758,68.74825876,68.95693421,69.24236443,69.61660982,69.56070769,69.92478534,70.36681123,70.31699272,70.66209698,70.90929098,71.10459976,71.07708152,71.07935235,71.48605962,71.22271161,71.30461644,71.5195365,71.94949418,72.4127505,72.62690389,72.78730301,73.0923973,73.61213484,73.7621953,73.90671906,74.36345422]
plt.plot(x, SBD, label='UB SDPD (SDPD-Benchmark)',color='forestgreen', marker='x',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, DBD, label='UB DPD (DPD-Benchmark)',color='lightsteelblue', marker='s', markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
#plt.plot(x00, back, label='UB exact DP',color='cyan', marker='^',markerfacecolor='none', linestyle='-', linewidth=1.5, markersize=5)
#plt.plot(x, SBD_ke, label='UB ',color='red', marker='o', markerfacecolor='none', linestyle='-.', linewidth=width, markersize=marksize)
#plt.plot(x, DBD, label='UB DPD-Benchmark',color='yellow', marker='s', markerfacecolor='none', linestyle=' ', linewidth=width, markersize=marksize)
#plt.plot(x, ADP, label='UB AFF',color='black', marker='x', linestyle='-', linewidth=1.5, markersize=5)


plt.plot(x, SBD_simu, label='Policy SDPD',color='forestgreen', marker='o', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, SBD_ke_simu, label='Policy SDPD-Benchmark',color='red', marker='o', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, DBD_simu, label='Policy DPD',color='lightsteelblue', marker='s', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, DBD_ke_simu, label='Policy DPD-Benchmark',color='yellow', marker='s', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, ADP_simu, label='Policy AFF',color='black', marker='x', linestyle=':', linewidth=1.5, markersize=5)
plt.plot(x, DP_simu, label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, sbADP_simu, label='Policy sbADP (M=100)',color='salmon', marker='+', linestyle=':', linewidth=1.5, markersize=5)
#plt.plot(x, back_simu, label='Policy BE',color='cyan', marker='^', linestyle=':', linewidth=1.5, markersize=5)

x0=[i*0.5 for i in range(11)]
plt.xticks(x0)

plt.ylim(60, 85)
plt.xlabel('Quality Index of Product Type 3')
plt.ylabel('Revenue')
#plt.title('Homogeneous Seats')
plt.legend(ncol=2)

plt.show()