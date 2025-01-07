import numpy as np
import matplotlib.pyplot as plt

choice=0
c=8

def read():
    folder_path = "ADP"
    file_name = "ADP_NL_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    ADP = np.loadtxt(file_path)

    folder_path = "ADP"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    ADP_simu=np.loadtxt(file_path)

    folder_path = "DBD"
    file_name = "DBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_ke=np.loadtxt(file_path)

    folder_path = "DBD"
    file_name = "ke_simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_ke_simu=np.loadtxt(file_path)

    folder_path = "DBD"
    file_name = "DBD_NL_results0capacity" + str(c) + "2.txt"
    file_path = f"{folder_path}/{file_name}"
    DBD=np.loadtxt(file_path)

    for i in range(len(DBD_ke)):
        if DBD_ke[i]>DBD[i]:
            DBD_ke[i] = DBD[i]

    folder_path = "DBD"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_simu=np.loadtxt(file_path)

    folder_path = "SBD"
    file_name = "SBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke=np.loadtxt(file_path)

    folder_path = "SBD"
    file_name = "ke_simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_ke_simu=np.loadtxt(file_path)

    folder_path = "SBD"
    file_name = "SBD_NL_results" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD=np.loadtxt(file_path)

    for i in range(len(SBD_ke)):
        if SBD_ke[i]>SBD[i]:
            SBD_ke[i] = SBD[i]

    folder_path = "SBD"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_simu=np.loadtxt(file_path)

    folder_path = "DP"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DP_simu=np.loadtxt(file_path)

    folder_path = "sbADP"
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
back_simu=[66.3037905360154,66.3076630366673,66.1844208875921,66.311977846484,66.5763890332364,66.7379725519168,66.8653144351569,66.6704155504446,66.7683456432169,66.9209150894345,67.0803409137416,66.9628714153628,67.1320846319956,67.2944794706529,67.475326625114,67.84276651227,68.1488877284806,68.3602618329351,68.30898212694,68.2073673332018,68.3135544520792,68.4767497856754,68.7173473336762,68.8180650562869,69.2026855249813,69.0801467477101,69.0670333163683,69.2502319856361,69.4866431459213,69.472798832494,69.7629778485827,69.8026852555247,69.8671252544415,69.985259084133,70.221783393474,70.3512158622362,70.6199152137167,70.8089093663304,70.7687785495535,71.0750031460329,71.2894660072374,71.2251121567624,71.4285557977877,71.9016865435222,72.3474874576786,72.6576816949845,72.8464895013263,72.9148958620734,73.4471833684394,73.861025216043,74.0012135587128]
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
plt.plot(x, back_simu, label='Policy BE',color='cyan', marker='^', linestyle=':', linewidth=1.5, markersize=5)

x0=[i*0.5 for i in range(11)]
plt.xticks(x0)

plt.ylim(60, 85)
plt.xlabel('Quality Index of Product Type 3')
plt.ylabel('Revenue')
#plt.title('Homogeneous Seats')
plt.legend(ncol=2)

plt.show()