import matplotlib.pyplot as plt
import numpy as np
x0=[8,16,24,32,40,48,56,64,72,80]

def regular(a):
    a = [a[i] / (i + 1) / 8 for i in range(len(a))]
    return a
def fill(a):
    n=len(x0)-len(a)
    if n <0:
        raise Exception("x too short")
    a=np.append(a,np.full(n,np.nan))
    return a

def check(a):
    a=regular(a)
    a=fill(a)
    return a

plt.rcParams['font.size'] = 14
width=2
marksize=10
#time plot
#homo 0
y0=[0,10000,20000,30000,40000,50000]

DBD_time=fill(np.array([345.8553593, 2878.876501,18782.09059,25362.08651,44290.1748]))
DBD_ke_time=fill(np.array([658.5907857,4339.441035,38366.94471]))
#DBD_ke_time=fill(np.array([0,528.3984311,3908.125157,28905.89837,67202.10762]))
SBD_time=fill(np.array([9.282173157,85.30460787,230.6744361,381.0603254,916.5976729,2383.854894,2955.83834,4134.055677,6772.866756,7330.182154]))
SBD_ke_time=fill(np.array([27.80594444, 200.3940618,818.3333614,1580.680403,3113.017182,5791.096684,9138.81043,14124.25367,17954.55252,23513.62532]))
aALP_time=fill(np.array([130,431,9461,13816,15744,35372]))
#aALP_time=fill(np.array([0,130,431,9461,13816,15744,35372,84648.66592]))
sbADP_time=fill(np.array([52,105,219,462,1106,2196,4163,7478,13430,21077]))
DLPflex=fill(np.array([0.07914900779724121, 0.12064266204833984, 0.18642735481262207, 0.23073077201843262, 0.2664494514465332, 0.32999467849731445, 0.3990328311920166, 0.41704225540161133, 0.48668670654296875, 0.5458629131317139]))


plt.figure(figsize=(10, 6))

plt.plot(x0, SBD_time, label='SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_time, label='DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, SBD_ke_time, label='SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_ke_time, label='DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, aALP_time, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DLPflex, label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
plt.plot(x0, sbADP_time, label='sbADP (M=100)',color='salmon', marker='<', linestyle='-', linewidth=width, markersize=marksize)

plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Run Time',fontsize=14)
#plt.title('Run Time (seconds) for Homogeneous Seats')
plt.legend()
plt.xticks(x0)
plt.yticks(y0)
#aw2
DBD_time=fill(np.array([286.1224339, 2206.501436,11554.20849,29286.40663,36663.69675]))
DBD_ke_time=fill(np.array([455.0947847,3963.751479,22157.11438]))
#DBD_ke_time=fill(np.array([0,350.4428675,3541.74771,19689.52314,75754.19555]))
SBD_time=fill(np.array([9.266380072,58.73397589,227.493367,475.9389737,1069.428683,1530.37579,2952.668684,4048.5019,6692.840921,7903.385561]))
SBD_ke_time=fill(np.array([28.9281497, 202.2409129,753.9688015,1695.522408,3250.571537,5647.914413,9329.398163,13084.49453,20125.33135,21868.89895]))
aALP_time=fill(np.array([155,1460,2468,14428,15547,38560]))
#aALP_time=fill(np.array([0,155,1460,2468,14428,15547,38560,327851.9573]))
sbADP_time_a=fill(np.array([264,528, 1127, 2500,6023,12373, 23225, 42083]))
sbADP_time=fill(np.array([0,264,528, 1127, 2500,6023,12373, 23225, 42083,73170]))
DLPflex=fill(np.array([0.09188175201416016, 0.13351035118103027, 0.1693112850189209, 0.22303390502929688, 0.29317712783813477, 0.3276815414428711, 0.3832204341888428, 0.47506117820739746, 0.5571846961975098, 0.5986065864562988]))

plt.figure(figsize=(10, 6))


plt.plot(x0, SBD_time, label='SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_time, label='DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, SBD_ke_time, label='SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_ke_time, label='DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, aALP_time, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DLPflex, label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
plt.plot(x0, sbADP_time, label='sbADP (M=500)',color='salmon', marker='<', linestyle='-', linewidth=width, markersize=marksize)

plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Run Time',fontsize=14)
#plt.title('Run Time (seconds) for Window Aisle Seats')
plt.legend()
plt.xticks(x0)
plt.yticks(y0)
#hete1
DBD_time=fill(np.array([294.8786089, 3195.408029,11719.0131,25402.84648]))
DBD_ke_time=fill(np.array([440.8157694,4902.045527,26468.85311]))
#DBD_ke_time=fill(np.array([0,312.4536974,4434.602535,24010.88961,85548.72674]))
SBD_time=fill(np.array([8.89052844,91.70611191,207.0481138,501.7657037,916.0186565,1568.886369,2924.425771,4099.07664,6735.244641,7484.137447]))
SBD_ke_time=fill(np.array([34.61270237, 205.1068294,739.0670755,1736.487378,3183.843248,5005.444528,9073.568395,14029.98374,18665.1997,21966.71746]))
aALP_time=fill(np.array([187,1113,2458,14428,43821]))
#aALP_time=fill(np.array([0,187,1113,2458,14428,43821,80617.40491,134931.6026]))
sbADP_time_h=fill(np.array([198,582,1236,2674,6583,13063,23957]))
sbADP_time=fill(np.array([0,198,582,1236,2674,6583,13063,23957,26668,72921]))
DLPflex=fill(np.array([0.07342720031738281, 0.11832618713378906, 0.181640625, 0.24383139610290527, 0.2663893699645996, 0.33180975914001465, 0.4287388324737549, 0.46768832206726074, 0.5236818790435791, 0.5246953964233398]))

plt.figure(figsize=(10, 6))

plt.plot(x0, SBD_time, label='SDPD',color='forestgreen', marker='o',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_time, label='DPD',color='lightsteelblue', marker='s',markerfacecolor='none', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, SBD_ke_time, label='SDPD-Benchmark',color='orange', marker='x', markerfacecolor='none',linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, DBD_ke_time, label='DPD-Benchmark',color='crimson', marker='+',markerfacecolor='none', linestyle='--', linewidth=width, markersize=marksize)
plt.plot(x0, aALP_time, label='AFF',color='black', marker='^', linestyle='-', linewidth=width, markersize=marksize)
plt.plot(x0, DLPflex, label='DPP',color='violet', marker='.', linestyle='-', linewidth=width, markersize=10)
plt.plot(x0, sbADP_time, label='sbADP (M=500)',color='salmon', marker='<', linestyle='-', linewidth=width, markersize=marksize)


plt.xlabel('Bus Size',fontsize=14)
plt.ylabel('Run Time',fontsize=14)
#plt.title('Run Time (seconds) for Heterogeneous Seats')
plt.legend()
plt.xticks(x0)
plt.yticks(y0)


plt.show()