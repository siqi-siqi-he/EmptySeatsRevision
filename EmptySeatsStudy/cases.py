import numpy as np

def homo_seats(c):
    #low demand
    #Create a scenario that, demand is low, and it forces the policy maker to price low in the beginning to attract customers
    b=np.array([k*0.2+0.2 for k in range(3)])
    a2=np.array([0.4 for k in range(c)])
    a3=np.array([0.6 for k in range(c)])
    a1=0.2
    tau=[1,0.4,0.6]
    return a1,a2,a3,b, tau