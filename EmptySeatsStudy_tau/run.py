import cvxpy as cp
import numpy as np
import decomposition.Branch_Bound as BB
import cases as cases
import math
import time
c=8
T=c*2
choice=0


a1, a2, a3, b, tau = cases.homo_seats(c)
results=[0]*51

for step in range(51):
    folder_path = "DBD"
    file_name = "DBD_NL_DPtable" + str(choice) + "capacity" + str(c) + "step" + str(step) + ".txt"
    file_path = f"{folder_path}/{file_name}"
    V=np.loadtxt(file_path)

    results[step] = V[T, c] + sum(v[T, :])

folder_path = "DBD"
file_name = "DBD_NL_ke_results" + str(choice) + "capacity" + str(c) + ".txt"
file_path = f"{folder_path}/{file_name}"
np.savetxt(file_path, results)