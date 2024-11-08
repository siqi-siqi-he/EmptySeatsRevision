import numpy as np

def prob3(r1,r2,r3,a1,a2,a3,b,tau,c,a0=0):
    p3=[np.exp((a3[k]-b[2]*r3[k])*tau[2])/(np.exp(a0)+np.exp((a1-b[0]*r1)*tau[0])+(np.sum([np.exp(a2[j]-b[1]*r2[j]) for j in range(c)])**tau[1])
                   +(np.sum([np.exp(a3[j]-b[2]*r3[j]) for j in range(c)])**tau[2])) for k in range(c)]
    
    return p3