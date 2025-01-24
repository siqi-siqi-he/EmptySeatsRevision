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
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    ADP_simu=np.loadtxt(file_path)

    # folder_path = "DBD"
    # file_name = "DBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    # file_path = f"{folder_path}/{file_name}"
    # DBD_ke=np.loadtxt(file_path)

    folder_path = "DBD_ke"
    file_name = "ke_simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
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
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DBD_simu=np.loadtxt(file_path)

    # folder_path = "SBD"
    # file_name = "SBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
    # file_path = f"{folder_path}/{file_name}"
    # SBD_ke=np.loadtxt(file_path)

    folder_path = "SBD_ke"
    file_name = "ke_simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
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
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    SBD_simu=np.loadtxt(file_path)

    folder_path = "DP"
    file_name = "simu_means_choice" + str(choice) + "capacity" + str(c) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    DP_simu=np.loadtxt(file_path)

    return SBD_simu, SBD, DBD, SBD_ke_simu, DBD_simu, DBD_ke_simu, ADP_simu, DP_simu


width=4
marksize=10
plt.rcParams['font.size'] = 20
fs = 20
figs = (15,10)
SBD_simu, SBD, DBD, SBD_ke_simu, DBD_simu, DBD_ke_simu, ADP_simu, DP_simu=read()

sbADP_simu=[74.1537312452605,74.2367933642168,74.5383989598361,74.8442545126494,75.1570725033662,75.0830170108328,75.4341720860732,75.7011555065234,75.8352737386407,76.0772662529862,76.3562362776216,76.5878620677655,76.773552771101,76.8959064574097,77.4118605736404,77.6179995406188,77.9706865477905,78.1153102087869,78.3487376249387,78.3602530879663,78.5946946587639,78.8238807219534,79.311856385674,79.8433645054052,80.3518253965741,80.5649827361087,80.9398113541407,81.3036453899208,81.2963397236947,81.7532539740003,82.1007782613648,82.378735140952,82.6552972302016,82.7920814151021,83.0856500524482,83.257937056669,83.7008661224365,84.0528972721104,84.186873301251,84.5902126738983,85.3316785492079,85.6095681605576,86.1303694761636,86.7267217753216,87.0292587358855,87.3363777445834,87.696473390251,87.9757251988349,88.4821620378651,88.8310713210367,89.2554821634013]

plt.figure(figsize=figs)
x=[i*0.1 for i in range(51)]
x00=[0,1,2,3,4,5]
back_simu=[74.7732123655279,74.9417963898215,75.1651295328667,75.5118061282089,75.823079051862,75.8010153138733,76.1002091519976,76.5395916047871,76.7416319625379,76.6738623946228,76.7161933253265,76.6028456538968,76.7852466666737,76.9622828838417,77.1779316738666,77.2709083396488,77.9763961663297,78.2738208894418,78.4227487775482,78.7544922971352,79.1240400090436,79.3658189873389,79.8233668497888,80.1578301200962,80.1218637092999,80.4053481001591,80.8717347220173,81.0088847421423,81.3902944334303,81.6919857033516,82.1277061666842,82.4595197422291,82.6563744912798,83.0935900039287,83.3072481357243,83.4773802139676,83.8249262911472,84.2383760113857,84.4978984660582,84.9817901052413,85.3892634184886,85.5600974669735,86.0492520099222,86.4387934166676,86.9762318168815,87.3745380503054,87.7439868313437,88.0673089695013,88.1692686495114,88.5305820903084,88.9679680811535]

plt.plot(x, SBD, label='UB DPD (DPD-Benchmark)',color='forestgreen', marker='x',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x, DP_simu, label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_simu, label='Policy DPD',color='forestgreen', marker='o', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, DBD_simu, label='Policy TD-DPD',color='lightsteelblue', marker='s', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, back_simu, label='Policy BE',color='cyan', marker='^', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, DBD, label='UB TD-DPD (TD-DPD-Benchmark)',color='lightsteelblue', marker='s', markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
#plt.plot(x00, back, label='UB exact DP',color='cyan', marker='^',markerfacecolor='none', linestyle='-', linewidth=1.5, markersize=5)
#plt.plot(x, SBD_ke, label='UB ',color='red', marker='o', markerfacecolor='none', linestyle='-.', linewidth=width, markersize=marksize)
#plt.plot(x, DBD, label='UB DPD-Benchmark',color='yellow', marker='s', markerfacecolor='none', linestyle=' ', linewidth=width, markersize=marksize)
#plt.plot(x, ADP, label='UB AFF',color='black', marker='x', linestyle='-', linewidth=1.5, markersize=5)
plt.plot(x, ADP_simu, label='Policy AFF',color='black', marker='x', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, SBD_ke_simu, label='Policy DPD-Benchmark',color='red', marker='o', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, DBD_ke_simu, label='Policy TD-DPD-Benchmark',color='yellow', marker='s', linestyle=':', linewidth=width, markersize=marksize)

plt.plot(x, sbADP_simu, label='Policy sbADP (M=100)',color='salmon', marker='+', linestyle=':', linewidth=width, markersize=marksize)



x0=[i*0.5 for i in range(11)]
plt.xticks(x0)

plt.ylim(73, 103)
plt.xlabel('Quality index of product type 3',fontsize=fs)
plt.ylabel('Revenue',fontsize=fs)
#plt.title('Homogeneous Seats')
plt.legend(ncol=2)

plt.show()