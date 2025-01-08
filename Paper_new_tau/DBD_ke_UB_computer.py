import numpy as np

c=8
choice=2
folder_path = "ADP"
file_name = "v_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) +".txt"
file_path = f"{folder_path}/{file_name}"
v=np.loadtxt(file_path)
file_name = "w_value_ADP_NL_CVXPY_choice" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
w=np.loadtxt(file_path)

UB=100000
T=c*2
for i in range(c):
    folder_path = "DBD_ke"
    file_name = "DBD_NL_ke_results" + str(choice) + "capacity" + str(c) + "compo"+str(i)+".txt"
    file_path = f"{folder_path}/{file_name}"
    V=np.loadtxt(file_path)
    temp_UB = V[T,1] + sum(v[T, :]) - v[T, i] + w[T] * c
    if UB>temp_UB:
        UB=temp_UB
print(UB)