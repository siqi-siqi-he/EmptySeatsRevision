import numpy as np

def homo_seats(c):
    #low demand
    #Create a scenario that, demand is low, and it forces the policy maker to price low in the beginning to attract customers
    b=[0.2,0.4,0.6]
    a2=np.array([0.4 for k in range(c)])
    a3=np.array([0.6 for k in range(c)])
    a1=0.2
    tau=[1,0.4,0.4]
    return a1,a2,a3,b, tau

def incre_seats(c):
    #low demand
    #Create a scenario that, demand is low, and it forces the policy maker to price low in the beginning to attract customers
    b=np.array([k*0.2+0.2 for k in range(3)])
    a2=np.array([0.4+0.01*k for k in range(c)])
    a3=np.array([0.605+0.02*((k+2)//2-1) for k in range(c)])
    a1=0.2
    tau=[1,0.4,0.6]
    return a1,a2,a3,b, tau

def aw_seats(c):
    a1=0.2
    a2=[0.4 for j in range(c)]
    for j in range(c):
        if j%4==0 or j%4==3:
            a2[j]=0.5
    a3=[0.6 for j in range(c)]
    b=np.array([k*0.2+0.2 for k in range(3)])
    tau = [1, 0.4, 0.6]

    return a1, a2, a3, b, tau

def homo_seats_wo3(c):
    # low demand
    # Create a scenario that, demand is low, and it forces the policy maker to price low in the beginning to attract customers
    b = np.array([k * 0.2 + 0.2 for k in range(3)])
    a2 = np.array([0.4 for k in range(c)])
    a3 = np.array([-500 for k in range(c)])
    a1 = 0.2
    tau = [1, 0.4, 0.6]
    return a1, a2, a3, b, tau

def incre_seats_wo3(c):
    #low demand
    #Create a scenario that, demand is low, and it forces the policy maker to price low in the beginning to attract customers
    b=np.array([k*0.2+0.2 for k in range(3)])
    a2=np.array([0.4+0.01*k for k in range(c)])
    a3=np.array([-500 for k in range(c)])
    a1=0.2
    tau=[1,0.4,0.6]
    return a1,a2,a3,b, tau

def aw_seats_wo3(c):
    a1=0.2
    a2=[0.4 for j in range(c)]
    for j in range(c):
        if j%4==0 or j%4==3:
            a2[j]=0.5
    a3=[-500 for j in range(c)]
    b=np.array([k*0.2+0.2 for k in range(3)])
    tau = [1, 0.4, 0.6]

    return a1, a2, a3, b, tau