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

#choice 4
SBD_UB=np.array([68.6715819035891,158.999831998647,258.879837610151,365.629534999981,477.975412909706,595.177922039927,716.74938252263,842.334912906324,971.668948389566,1104.54486440026])#done,10 #done
DBD_UB=fill(np.array([68.6795620072274,159.140833508944,258.997945936926,365.739462060416]))
SBDb_UB=np.array([68.6715819035891,158.999831998647,258.879837610151,365.629534999981,477.975412909706,595.177922039927,716.74938252263,842.334912906324,971.668948389566,1104.54486440026])#done,10 #done
#SBDb_UB_woy=np.array([73.4808988822258,165.935911900997,267.284163308995,375.17681565103,488.470401867839,606.490207670472,728.779126008739,855.004026835559,984.916921375032,1118.3232773819])#done,10 #done
DBDb_UB=fill(np.array([]))
aALP=fill(np.array([73.1648187824978,165.481227121692,266.64690363421,374.280216708899,487.151061348171]))#done,7
DLPflex=fill(np.array([73.9884542687187,166.498083154143,267.895439478212,375.836449240972,489.179553435683,607.248130545888,729.586667194739,855.865528932043,985.835661228318,1119.30166244716]))#done,10 #done
'''
SBD_simu=fill(np.array([64.45809368804186,143.3558442032371,227.64853806529132,323.2557729131836,426.1535249589765,522.8804291843757,617.4411673390339,719.2647080306211,819.6401945886668,933.0951726028303]))#done,10 #done
SBD_ke_simu=fill(np.array([64.13492400371939,144.6330699263746,228.88853776814994,323.9639818161254,427.02098641286875,523.9407636260163,615.945781929311,721.1228944880239,823.4462667699713]))#done,10 #done
DBD_simu=fill(np.array([64.45809051095915,143.35673062213235,227.62414229682327,323.6187261277688,426.93141009009065,522.5696358501605]))
DBD_ke_simu=fill(np.array([64.45641323136032,143.43540471957166,225.21941510947264,323.5581634821159]))
aALP_simu=fill(np.array([61.78427238957758, 138.61018319345382,226.14938995316416, 320.5578006315373, 421.05821247178477, 516.9756549928264]))
sbADP_simu=fill(np.array([64.24429505,143.5083651,207.1127125,303.1144468,488.4390255,702.6241617,598.6003484,805.3462068]))#done,10 #done
DLP_simu=fill(np.array([61.42340582807179,138.77384748603984,225.31256181267935,319.40060639730126,419.9542325879878,515.679796399494,611.4473875845814,715.54296515527,814.386362157405,924.9763364853288])) #done
'''
#pricing policies
plt.figure(figsize=(10, 6))
'''
plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=100)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
'''
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
plt.title('Pricing Policies for hom_de_in_seats')
plt.legend(ncol=2)


#choice 5
SBD_UB=fill(np.array([68.2417792891811,157.125405898257,254.399615042252,357.286937677139,464.437725259773,575.044921767565,688.569207247698,804.623485419348,922.916059754153,1043.21923788135]))
DBD_UB=fill(np.array([68.3276274717675,157.179376820253,254.6065222414,357.666471763378]))
SBDb_UB=fill(np.array([68.2417792891811,157.125405898257,254.399615042252,357.286937677139,464.437725259773,575.044921767565,688.569207247698,804.623485419348,922.916059754153,1043.21923788135]))
#SBDb_UB_woy=fill(np.array([73.0423366698426,164.038409933635,262.758985742266,366.75922800273,474.819750250712,586.198952932876,700.395153193808,817.045030978347,935.872814484985,1056.66345147933]))
DBDb_UB=fill(np.array([]))
aALP=fill(np.array([72.7291883494042,163.607569896348,262.224464801718,366.138918220246,474.135437089136]))
DLPflex=fill(np.array([73.5463378022181,164.584524119174,263.334627070838,367.358251532822,475.438554594031,586.835307110368,701.047417623319,817.712069915188,936.553783880678,1057.3559761638]))
'''
SBD_simu=fill(np.array([64.91859205342483,144.44060081632043,231.96201179271165,330.32877287934826,439.8849827504754,540.5428376987037,642.7059835868749,753.1165448065798,864.2548956512823,988.4190725847928]))
SBD_ke_simu=fill(np.array([64.27717819255412,145.93353891466737,232.72466824411995,331.83877206044866,439.6551335279278,542.007203892988,641.6674102498348,756.5277381765294,866.123581246961,990.1219900151912]))
DBD_simu=fill(np.array([64.90404381916471,144.70667822365897,232.03290941674229,330.4564738149191,440.15818282369617]))
DBD_ke_simu=fill(np.array([64.80777907288238,144.80001390304696,228.9944048354347,330.5568757934293]))
aALP_simu=fill(np.array([62.086445641567764,140.7567224599399,229.70741491852507,327.59495025447677,432.14190889441284]))
sbADP_simu=fill(np.array([49.78824925,145.3247069,213.2878257,313.8770599,484.9166528,620.5563582,693.5898398,578.9596915]))
DLP_simu=fill(np.array([61.82121816182555,140.3993730896393,229.28901692312252,326.9326027370243,431.84158643394056,533.1170004257158,637.013406882737,748.578756606833,857.3097408468957,978.5052961735919]))
'''

plt.figure(figsize=(10, 6))
'''
plt.plot(x, SBD_simu/SBD_UB, label='Policy SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_simu/SBD_UB, label='Policy DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, SBD_ke_simu/SBD_UB, label='Policy SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, DBD_ke_simu/SBD_UB, label='Policy DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x, aALP_simu/SBD_UB, label='Policy AFF',color='black', marker='^', linestyle=':', linewidth=width, markersize=marksize)
plt.plot(x,DLP_simu/SBD_UB,label='Policy DPP',color='violet', marker='.', linestyle=':', linewidth=width, markersize=marksize)
#plt.plot(x, sbADP_simu/SBD_UB, label='Policy sbADP (M=500)',color='salmon', marker='<', linestyle=':', linewidth=width, markersize=marksize)
'''
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
plt.title('Pricing Policies for hom_hom_in_seats')
plt.legend(ncol=2)

plt.show()