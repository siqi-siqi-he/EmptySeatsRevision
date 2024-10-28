import matplotlib.pyplot as plt
import numpy as np
x=[8,16,24,32,40,48,56,64,72,80]

#This plot is used to generate graphs for the submitted paper
def regular(a):
    a = [a[i] / (i + 1) / 8 for i in range(len(a))]
    return a
def fill(a):
    n=len(x)-len(a)
    if n <0:
        raise Exception("x too short")
    a=np.append(a,np.full(n,np.nan))
    return a

def check(a):
    a=regular(a)
    a=fill(a)
    return a

width=2
marksize=10
plt.rcParams['font.size'] = 14
#homo
SBD_UB=np.array([68.19310778,156.886865,253.8413461,356.2919813,462.8982624,572.8590697,685.6397759,800.8566531,918.2200106,1037.503977])#done,10 #done
DBD_UB=fill(np.array([68.20755474,157.4561226,253.9768079,356.6548509,467.6632643]))
SBDb_UB=np.array([68.19310778,156.886865,253.8413461,356.2919813,462.8982624,572.8590697,685.6397759,800.8566531,918.2200106,1037.503977])#done,10 #done
DBDb_UB=fill(np.array([68.20755474,157.4561226,253.9768079]))
aALP=fill(np.array([72.8403, 163.8807,261.7514,365.4608,472.9219,583.7837]))#done,7
DLPflex=fill(np.array([73.4672,164.2791, 262.6809, 366.2465, 473.7669,584.5070,697.9593,813.7966,931.7119,1051.4982]))#done,10 #done
SBD_simu=fill(np.array([64.45809368804186,143.3558442032371,227.64853806529132,323.2557729131836,426.1535249589765,522.8804291843757,617.4411673390339,719.2647080306211,819.6401945886668,933.0951726028303]))#done,10 #done
SBD_ke_simu=fill(np.array([64.13492400371939,144.6330699263746,228.88853776814994,323.9639818161254,427.02098641286875,523.9407636260163,615.945781929311,721.1228944880239,823.4462667699713]))#done,10 #done
DBD_simu=fill(np.array([64.45809051095915,143.35673062213235,227.62414229682327,323.6187261277688,426.93141009009065,522.5696358501605]))
DBD_ke_simu=fill(np.array([64.45641323136032,143.43540471957166,225.21941510947264,323.5581634821159]))
aALP_simu=fill(np.array([61.78427238957758, 138.61018319345382,226.14938995316416, 320.5578006315373, 421.05821247178477, 516.9756549928264]))
sbADP_simu=fill(np.array([64.24429505,143.5083651,207.1127125,303.1144468,488.4390255,702.6241617,598.6003484,805.3462068]))#done,10 #done
DLP_simu=fill(np.array([61.42340582807179,138.77384748603984,225.31256181267935,319.40060639730126,419.9542325879878,515.679796399494,611.4473875845814,715.54296515527,814.386362157405,924.9763364853288])) #done

#pricing policies
plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=100)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_UB/SBD_UB, label='UB SDPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD_UB/SBD_UB, label='UB DPD',color='lightsteelblue', marker='s', markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, SBDb_UB/SBD_UB, label='UB SDPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, DBDb_UB/SBD_UB, label='UB DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, aALP/SBD_UB, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DLPflex/SBD_UB, label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)

plt.xticks(x)
plt.ylim(0.85, 1.2)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Homogeneous Seats')
plt.legend(ncol=2)

#aw2
SBD_UB=fill(np.array([68.68107321,158.0035498,255.6408439,358.8093908,466.159216,576.8840245,690.4451588,806.4563508,924.6261723,1044.727118]))#done,10
DBD_UB=fill(np.array([68.85054879,158.8803758,255.7254855,358.9721947,466.3530573]))
SBDb_UB=fill(np.array([68.68107321,158.0035498,255.6408439,358.8093908,466.159216,576.8840245,690.4451588,806.4563508,924.6261723,1044.727118]))#done,10
DBDb_UB=fill(np.array([68.85054879,158.8803758,255.7254855]))
aALP=fill(np.array([73.5945, 165.3088, 263.5421,367.8312, 476.1498,587.6852765]))
DLPflex=fill(np.array([73.9693, 165.4193, 264.5128, 368.8046, 477.0683,588.5881, 702.8386, 819.4657,938.1937,1058.8027]))#done,10
SBD_simu=fill(np.array([65.10895773073004,143.80928378987676,229.64719502380683,325.8589842082737,430.35131838156354,526.2531429955628,620.5608891545172,724.4000404771156,826.1283734223487,938.6955120764221]))#done,10
SBD_ke_simu=fill(np.array([64.33205266756296,145.18591339173017,230.60421039533148,326.8093020384833,429.5995718363687,527.6350448525543,619.912427464435,725.3625526412491,828.2820129827884,940.7906955727236]))#done,10
DBD_simu=fill(np.array([65.05083854903685,144.0057210516991,229.7121247309898,325.7704926225174,430.304856632932]))
DBD_ke_simu=fill(np.array([65.01550890761612,143.97220255866603,225.97091928164133,325.7547260396516]))
aALP_simu=fill(np.array([62.07747180347792,140.36114203554246,227.65609499276587, 322.803462225222, 423.4051142872868,519.9610140033883]))
sbADP_simu=fill(np.array([64.49672941,144.5956576,209.9636685,306.1695418,493.7417763,708.2145221,601.1257444,797.1170085]))
DLP_simu=fill(np.array([61.82903035086845,139.98436266670598,226.95453511006673,321.7007977616657,422.6028381149527,519.1756560576766,616.2971665792353,719.8185365743251,820.7219294325531,930.8937186940685]))

plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_UB/SBD_UB, label='UB SDPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD_UB/SBD_UB, label='UB DPD',color='lightsteelblue', marker='s', markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, SBDb_UB/SBD_UB, label='UB SDPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, DBDb_UB/SBD_UB, label='UB DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, aALP/SBD_UB, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DLPflex/SBD_UB, label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)

plt.xticks(x)
plt.ylim(0.85, 1.2)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Window Aisle Seats')
plt.legend(ncol=2)

#hete1
SBD_UB=fill(np.array([68.5763,158.7802,258.5238,365.1316,477.3276,594.3414,715.6614,840.9240878,969.864613,1102.264783]))
DBD_UB=fill(np.array([68.73937214,159.3815067,258.4787941,365.0597056]))
SBDb_UB=fill(np.array([68.5763,158.7802,258.5238,365.1316,477.3276,594.3414,715.6614,840.9240878,969.864613,1102.264783]))
DBDb_UB=fill(np.array([68.73937214,159.3815067,258.4787941]))
aALP=fill(np.array([73.2709,165.9500,266.4274,374.1611,487.1153]))
DLPflex=fill(np.array([73.8903,166.2737,267.5327,375.3259,488.4989,606.3513,728.4132,854.347,983.8986,1116.869]))
SBD_simu=fill(np.array([64.91859205342483,144.44060081632043,231.96201179271165,330.32877287934826,439.8849827504754,540.5428376987037,642.7059835868749,753.1165448065798,864.2548956512823,988.4190725847928]))
SBD_ke_simu=fill(np.array([64.27717819255412,145.93353891466737,232.72466824411995,331.83877206044866,439.6551335279278,542.007203892988,641.6674102498348,756.5277381765294,866.123581246961,990.1219900151912]))
DBD_simu=fill(np.array([64.90404381916471,144.70667822365897,232.03290941674229,330.4564738149191,440.15818282369617]))
DBD_ke_simu=fill(np.array([64.80777907288238,144.80001390304696,228.9944048354347,330.5568757934293]))
aALP_simu=fill(np.array([62.086445641567764,140.7567224599399,229.70741491852507,327.59495025447677,432.14190889441284]))
sbADP_simu=fill(np.array([49.78824925,145.3247069,213.2878257,313.8770599,484.9166528,620.5563582,693.5898398,578.9596915]))
DLP_simu=fill(np.array([61.82121816182555,140.3993730896393,229.28901692312252,326.9326027370243,431.84158643394056,533.1170004257158,637.013406882737,748.578756606833,857.3097408468957,978.5052961735919]))


plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_UB/SBD_UB, label='UB SDPD',color='forestgreen', marker='o',markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DBD_UB/SBD_UB, label='UB DPD',color='lightsteelblue', marker='s', markerfacecolor='none',linestyle='-', linewidth=width, markersize=10)
plt.plot(x, SBDb_UB/SBD_UB, label='UB SDPD-Benchmark',color='orange', marker='x',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, DBDb_UB/SBD_UB, label='UB DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=10)
plt.plot(x, aALP/SBD_UB, label='UB AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=10)
plt.plot(x, DLPflex/SBD_UB, label='UB DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)

plt.xticks(x)
plt.ylim(0.85, 1.2)
plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Revenue',fontsize=14)
plt.title('Pricing Policies for Heterogeneous Seats')
plt.legend(ncol=2)

plt.show()