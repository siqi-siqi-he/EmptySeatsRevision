import matplotlib.pyplot as plt
import numpy as np
x=[8,16,24,32,40,48,56,64,72,80]

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
SBD_simu=fill(np.array([64.7844594477173,143.74435856485456,230.88775249018883,320.9056885873656,417.02491756749833,505.4333590745147,607.5826057273233,719.0147250369787,816.2714141970114, 935.4162050870052]))
SBD_ke_simu=fill(np.array([65.40562751449072,144.31078163136414,232.08725249253698,321.76669078802774,418.0399106673806,505.9540173478857,608.0916117898857,721.1636032939725,816.6135209407907,937.7666834792665]))
DBD_simu=fill(np.array([64.78457014494427, 143.75756509391104, 230.85438746899072, 321.2700808318734, 417.2067457682034]))
DBD_ke_simu=fill(np.array([64.78280304296386,143.75463171295667,227.24986116515316, 321.24471024694884]))
aALP_simu=fill(np.array([64.24432143260042, 141.9893251664817, 227.87668848991788, 316.8850452978404, 411.4522774296492, 502.59316874600614]))
sbADP_simu=fill(np.array([64.39882702,143.4395477,228.3903853,322.2731476,425.1693205,524.6126715,617.0372468,725.5294735,823.7657819,929.3082735]))#done,10 #done
DLP_simu=fill(np.array([63.15657399068582, 141.26455657904478, 227.9128414992202, 317.3645549357825, 410.02728791927524, 502.0846899340256, 604.3701985421726, 713.3806687833422, 808.9917280171077,927.5455330706978])) #done

#pricing policies
plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=100)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
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
SBD_simu=fill(np.array([65.1137321662178, 144.6595425674769, 232.31596402926894, 323.94679585784166, 419.65315993160846, 509.2974695307869, 612.0381547595799, 723.9253376068085, 820.8240864827378, 941.7884587791731]))
SBD_ke_simu=fill(np.array([65.66900576019005, 145.38438661955786, 233.4417549061405, 323.89811245011276, 420.16648123788764, 509.88153761396376, 612.9192985961913, 724.1197439944955, 821.5198891591658, 943.9995086027025]))
DBD_simu=fill(np.array([65.1443783558015, 145.24327847217228, 232.51311292998574, 324.1144053237597, 420.02860121636513]))
DBD_ke_simu=fill(np.array([65.07574083194945, 144.93991413501885, 229.1035923342615, 323.66427635795276]))
aALP_simu=fill(np.array([64.43728398945966, 143.16528035127817, 229.96752996274643, 319.89408206865403, 414.37270142961574, 505.0830919668856]))
sbADP_simu=fill(np.array([64.72927981,144.6300028,224.3224936,314.4368287,415.1331849,512.8844123,609.9621112,342.0380058,728.5255798,924.6603628]))
DLP_simu=fill(np.array([63.52199223201822, 141.90701424367688, 229.7801899118867, 319.1074813684044, 413.72967681384966, 505.4918625151584, 608.883975782836, 718.781370667185, 814.4778370673899,933.2755185251626]))

plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
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
SBD_simu=fill(np.array([65.02270252655944, 145.43106222341947, 234.6481029901588, 328.9530176496257, 429.42309051401804, 524.7555072649145, 632.0244613681923, 753.6239210050913, 858.3742953836594, 992.3130694295322]))
SBD_ke_simu=fill(np.array([65.65200599109528, 146.1383998492798, 235.92195166257997, 329.4656478862897, 429.98119937651455, 524.7212250016369, 634.1047530749912, 754.9477018435795, 860.3295800924002, 993.4567156333042]))
DBD_simu=fill(np.array([65.03233434782382, 145.68712049294297, 234.63373492681737, 329.0972992027019, 429.316804837577]))
DBD_ke_simu=fill(np.array([65.09040618485795, 145.46807763746853, 231.32470999291462, 328.6036978737251]))
aALP_simu=fill(np.array([64.62209078610763, 143.4636617666422, 232.1046227135223, 324.9791340244749, 423.4695323976028, 520.0383206062221]))
sbADP_simu=fill(np.array([0.573343904,144.6010238,232.4070841,330.166788,438.2675074,542.7296438,643.2595558,632.2683895,865.8429301,985.5522434]))
DLP_simu=fill(np.array([63.52420320111835, 142.61367250594898, 231.78120724062472, 324.59630954057207, 423.4722087976518, 520.8335696965938, 627.8965273158587, 747.5564105211687, 852.5643544967196,980.775359189735]))


plt.figure(figsize=(10, 6))

plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
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