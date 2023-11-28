#Dp_Dc_distribution
#show Dc distribution and Dp distribution and k
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import os
import csv
from scipy import optimize
import move as mv
font_manager.fontManager.addfont('/data/home/zzy/.fonts/ARIAL.TTF')
plt.rcParams['font.sans-serif']= 'ARIAL'

#f_1
def f_1(x,A,B):
    return A * x + B
#log func
def logF(va):
    if va == 0:
        return 'a'
    else:
        return math.log(va)


DataPath = '../Data/'
DataName = 'Dp_Dc_Conc.csv'
FigurePath = '../Figure/'

#Read data
Dp_Dc_Conc = []
#Read
with open(DataPath+DataName, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        float_row = [float(x) for x in row]
        Dp_Dc_Conc.append(float_row)
print(len(Dp_Dc_Conc))
Dp = Dp_Dc_Conc[:240]
Dc = Dp_Dc_Conc[240:480]
conc = Dp_Dc_Conc[480:]
print(len(Dp))
print(len(Dc))
print(len(conc))
File_number = len(Dp)
bin_gap=10
n_bin=60
#data cut
X = [i+1 for i in range(File_number)]
x = [i * bin_gap for i in range(n_bin)] # x,Dp from 0 to (Dp_max//50+1)*50 sep = 50    [0, 10, 20, ..., 300 ]   len = index +
def filter_fun(Dp,Dc,conc,bin_gap,n_bin,File_number,Dc_min,Dc_max):
    nCT = [[] for i in range(File_number)]
    for k in range(File_number):
        nCT[k] = [0 for i in range(n_bin)]
        for i in range(len(Dp[k])):
            for j in range(n_bin):
                if (j * bin_gap)<= Dp[k][i] - Dc[k][i] < (j * bin_gap + bin_gap) and Dc_min < Dc[k][i] <= Dc_max:
                    nCT[k][j] += conc[k][i]
    return nCT

nCT = [[] for i in range(File_number)]

nCT = filter_fun(Dp,Dc,conc,bin_gap,n_bin,File_number,Dc_min=0,Dc_max=10000)
#do move average

nCT_mv = list(map(list, zip(*nCT)))
for n in range(len(nCT_mv)):
    nCT_mv[n] = mv.smoothmoveavg(nCT_mv[n],24,'valid')
nCT_mv = list(map(list, zip(*nCT_mv)))
X1 = mv.smoothmoveavg(X,24,'valid')
X1 = [i-0.5 for i in X1]

print(X)
print(X1)
#  nomal fun
def nDp_to_filtered_nD(nDp):
    for i in range(len(nDp)):
        time_i_nDp_sum = sum(nDp[i])
        for j in range(len(nDp[i])):
            nDp[i][j] = nDp[i][j]/time_i_nDp_sum 
    return nDp

# do normal
norm_nCT = nDp_to_filtered_nD(nCT)
norm_nCT_mv = nDp_to_filtered_nD(nCT_mv)

#do ln
n_CT_ln = [[] for i in range(len(norm_nCT))]
for i in range(len(norm_nCT)):
    n_CT_ln[i] = list(map(logF, norm_nCT[i]))

n_CT_ln_mv = [[] for i in range(len(norm_nCT_mv))]
for i in range(len(norm_nCT_mv)):
    n_CT_ln_mv[i] = list(map(logF, norm_nCT_mv[i]))


#do filter n>0 and dp>=140
filtered_X_CT = [[] for i in range(len(n_CT_ln))]
filtered_nCT = [[] for i in range(len(n_CT_ln))]
for n in range(len(n_CT_ln)):
    filtered_X_CT[n] = [i for i, j in zip(x,n_CT_ln[n]) if j != 'a']
    filtered_nCT[n] = [j for i, j in zip(x,n_CT_ln[n]) if j != 'a']
filtered_X_CT_mv = [[] for i in range(len(n_CT_ln_mv))]
filtered_nCT_mv = [[] for i in range(len(n_CT_ln_mv))]
for n in range(len(n_CT_ln_mv)):
    filtered_X_CT_mv[n] = [i for i, j in zip(x,n_CT_ln_mv[n]) if j != 'a']
    filtered_nCT_mv[n] = [j for i, j in zip(x,n_CT_ln_mv[n]) if j != 'a']
#coculate B and K
k = [[] for i in range(len(norm_nCT))]
b = [[] for i in range(len(norm_nCT))]
for i in range(len(norm_nCT)):
#    print(len(filtered_X_CT[i]))
    if len(filtered_X_CT[i]) ==1:k[i],b[i]=0,0
    else:
        k[i],b[i] = optimize.curve_fit(f_1,filtered_X_CT[i],filtered_nCT[i])[0]

k_mv = [[] for i in range(len(norm_nCT_mv))]
b_mv = [[] for i in range(len(norm_nCT_mv))]


for i in range(len(norm_nCT_mv)):
#    print(len(filtered_X_CT_mv[i]))
    if len(filtered_X_CT_mv[i])== 1:k_mv[i],b_mv[i]=0,0
    else:
        k_mv[i],b_mv[i] = optimize.curve_fit(f_1,filtered_X_CT_mv[i],filtered_nCT_mv[i])[0]
print(max(k[-144:-1]))
print(min(k[-144:-1]))
print(len(X))
print(len(k))
plot_X = [i for i,j in zip(X,k) if j<0]
plot_k = [j for i,j in zip(X,k) if j<0]
plot_X1 = [i for i, j in zip(X1,k_mv) if j<0]
plot_k_mv = [j for i, j in zip(X1,k_mv) if j<0]


#print(X_mark)
#exit()
#plot1
FigureName = 'FigureS2_e'
fig, ax = plt.subplots(figsize=(2.3, 1.5),constrained_layout=False)
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.16)
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(right=0.95)
ax.plot(plot_X, plot_k , color='#0452d7',linewidth = 1,label = '$\mathrm{k}$')
ax.plot(plot_X1, plot_k_mv , color='#0452d7',linewidth = 1,linestyle='--',alpha=0.6,label = 'moving averaged $\mathrm{k}$')
ax.legend(loc='upper right',fontsize=5.5,frameon=True)#,bbox_to_anchor = (1,0.5))
ax.set_title('(e)', fontsize = 6.5 ,loc = 'left')
ax.set_xlabel('Time of simulation(h)', fontsize=5.5)
ax.set_ylabel(r"Slope of CT distribution", fontsize=5.5)
#ax.text(520,-3,'k = '+str(abs(round(k,3))),fontsize=10,horizontalalignment='center', verticalalignment='center')
ax.set_xlim([-5,245])
ax.set_ylim([-0.1,0])
#ax.set_yticks([-0.08,-0.06, -0.04,-0.02,0])
#ax.set_yticklabels([-10, -8,-6, -4, -2])
ax.set_xticks([24*i for i in range(11)])
ax.tick_params(axis='y', which='major', direction='in', width=0.8, length=3, labelsize=5.5, pad=1, rotation=0)
ax.tick_params(axis='x', which='major', direction='in', width=0.8, length=3, labelsize=5.5, pad=1, rotation=0)
ax.tick_params(top=True, bottom=True, left=True, right=True)
fig.savefig(FigurePath+FigureName,dpi=1000)


