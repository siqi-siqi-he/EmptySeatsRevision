import numpy as np

c=56
choice=2
T=2*c

folder_path = "SBD_NL"
file_name = "DLPnorm" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
w = np.loadtxt(file_path)
folder_path = "SBD_NL/"
file_name = "SBD_NL_DPtable" + str(choice) + "capacity" + str(c)  + ".txt"
file_path = f"{folder_path}/{file_name}"
Vy = np.loadtxt(file_path)

V = np.full((c,c*2 + 1, 2), 0.0)

folder_path = "SBD_NL_ke/seat" + str(c) + "choice" + str(choice)
for i in range(c):
    file_name = "SBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo" + str(i) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    V[i,:,:] = np.loadtxt(file_path)

UB_SBD_y=Vy[T,c]+sum(w)-w[c]
print(UB_SBD_y)
UB_DBD_ke=100000
for i in range(c):
    temp_UB=V[i,T,1]+sum(w)-w[i]+w[c]*(c-1)

    if temp_UB<=UB_DBD_ke:
        UB_DBD_ke=temp_UB
print(UB_DBD_ke)