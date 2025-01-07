import numpy as np
import os

np.random.seed(42) #42 test: 50, test2: 60
simu_time=1000
i=10
c=8*i
rand=np.random.rand(2*c,simu_time)
directory="RandomNumbers_test2"
os.makedirs(directory,exist_ok=True)
full_path=directory+"/rand.txt"
np.savetxt(full_path, rand)
