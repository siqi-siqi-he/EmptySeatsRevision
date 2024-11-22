import numpy as np

def read(c,choice):
    directory = "results/DLP"
    full_path = directory + "/" + str(choice) + str(c) + ".txt"
    v=np.loadtxt(full_path)
    res=np.mean(v)
    print(res)
    return res

read(48,2)