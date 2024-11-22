import numpy as np
import os

np.random.seed(42)
simu_time=1000
for i in range(1,11):
    c=8*i
    rand=np.random.rand(2*c,simu_time)
    directory="RandomNumbers/"+str(c)
    os.makedirs(directory,exist_ok=True)
    full_path=directory+"/rand.txt"
    np.savetxt(full_path, rand)
