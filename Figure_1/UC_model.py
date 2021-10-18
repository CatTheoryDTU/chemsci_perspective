import matplotlib.pyplot as plt
import numpy as np

xdefect = np.linspace(0,1,100)
diff = np.arange(0.02, 0.22, 0.02)
R = 8.617e-5 #eV/K
T = 298      #K
color_list = ['#BF3F3F','#F47A33','#FFE228','#7FBF3F','#3FBFBF','#3F7FBF','#3F3FBF','#7F3FBF','#BF3F7F','#BF3F3F','#333333','#000000']

for j in range(0,len(diff)):
    ratio = []
    num = []
    denom = []
    percent_xdefect = []
    percent_ratio = []
    for i in range(0,len(xdefect)):
        num.append(xdefect[i])
        denom.append(xdefect[i] + (1-xdefect[i])*(np.exp(-diff[j]/(R*T))))
        ratio.append(num[i]/denom[i])
        percent_xdefect.append(xdefect[i]*100)
        percent_ratio.append(ratio[i]*100)
    plt.plot(percent_xdefect,percent_ratio,linestyle='solid',color=color_list[j],label='{0:.2g}'.format(diff[j]),lw=2)

plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
plt.xlim(0,30)
plt.ylim(0,100)
plt.ylabel('\% of activity from under-coordinated sites', fontsize=12)
plt.xlabel('\% of under-coordinated sites', fontsize=12)

plt.plot(15, 92, marker='o', markersize=10,markeredgecolor='black',color="red")
plt.plot(3, 53, marker='o', markersize=10,markeredgecolor='black', color="red")
Pb_coverage = [3,15]
activity_reduction = [53,92]
yerr = [[9,1.5],[8,2]]

plt.errorbar(Pb_coverage,activity_reduction,yerr,color = "k",ls = "None",capsize=3)

plt.plot(Pb_coverage,activity_reduction,'ro',markersize=8,markeredgecolor='black')


plt.savefig('UC-sites-percentage.eps')
