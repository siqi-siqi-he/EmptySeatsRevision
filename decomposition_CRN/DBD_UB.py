import numpy as np

c=16
choice=2
T=2*c

folder_path = "../BB"
file_name = "part_3_v_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
file_path = f"{folder_path}/{file_name}"
v = np.loadtxt(file_path)
print(v)
file_name = "part_3_w_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + "BB_ex.txt"
file_path = f"{folder_path}/{file_name}"
w = np.loadtxt(file_path)
print(w)
folder_path = "DBD_NL/"
file_name = "DBD_NL_DPtable" + str(choice) + "capacity" + str(c)  + ".txt"
file_path = f"{folder_path}/{file_name}"
Vy = np.loadtxt(file_path)

V = np.full((c,c*2 + 1, 2), 0.0)

folder_path = "DBD_NL_ke/seat" + str(c) + "choice" + str(choice)
for i in range(c):
    file_name = "DBD_NL_ke_DPtable" + str(choice) + "capacity" + str(c) + "compo" + str(i) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    V[i,:,:] = np.loadtxt(file_path)

UB_DBD_y=Vy[T,c]+sum(v[T,:])
print(UB_DBD_y)
UB_DBD_ke=100000
for i in range(c):

    temp_UB=V[i,T,1]+sum(v[T,:])-v[T,i]+w[T]*c

    if temp_UB<=UB_DBD_ke:
        UB_DBD_ke=temp_UB
print(UB_DBD_ke)