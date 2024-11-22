import numpy as np

for j in range(1,9):
    c=j*8
    for choice in range(3):
        directory = "results/DLP"
        full_path = directory + "/"+str(choice)+str(c)+".txt"
        DLP=np.loadtxt(full_path)
        directory = "results/SBD"
        full_path = directory + "/" + str(choice) + str(c) + ".txt"
        SBD = np.loadtxt(full_path)
        '''
        directory = "results/SBD_ke"
        full_path = directory + "/" + str(choice) + str(c) + ".txt"
        SBD_ke = np.loadtxt(full_path)
        directory = "results/DBD"
        full_path = directory + "/" + str(choice) + str(c) + ".txt"
        DBD = np.loadtxt(full_path)
        directory = "results/DBD_ke"
        full_path = directory + "/" + str(choice) + str(c) + ".txt"
        DBD_ke = np.loadtxt(full_path)            
        print("new")
        print(np.var(DLP-SBD,ddof=1))
        print(np.var(DLP,ddof=1)+np.var(SBD,ddof=1)-2*np.cov(DLP,SBD)[1,0])
        '''
