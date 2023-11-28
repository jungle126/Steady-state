#Dp_Dc_distribution
#show Dc distribution and Dp distribution and k
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import os
import csv
from scipy import optimize
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
FigurePath = '../Figure/CT_dis/'

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


print(X)
#  nomal fun
def nDp_to_filtered_nD(nDp):
    for i in range(len(nDp)):
        time_i_nDp_sum = sum(nDp[i])
        for j in range(len(nDp[i])):
            nDp[i][j] = nDp[i][j]/time_i_nDp_sum 
    return nDp

# do normal
norm_nCT = nDp_to_filtered_nD(nCT)
print(sum(norm_nCT[0])) #validation


#do ln
n_CT_ln = [[] for i in range(len(norm_nCT))]
for i in range(len(norm_nCT)):
    n_CT_ln[i] = list(map(logF, norm_nCT[i]))

#do filter n>0 and dp>=140
filtered_X_CT = [[] for i in range(len(n_CT_ln))]
filtered_nCT = [[] for i in range(len(n_CT_ln))]
for n in range(len(n_CT_ln)):
    filtered_X_CT[n] = [i for i, j in zip(x,n_CT_ln[n]) if j !='a']
    filtered_nCT[n] = [j for i, j in zip(x,n_CT_ln[n]) if j !='a' ]

#k,b=optimize.curve_fit(f_1,filtered_X_CT,filtered_nCT)[0]
#X_line = np.arange(min(filtered_X_CT),max(filtered_X_CT),0.01)
#Y_line = X_line*k+b

def R2_func(filtered_X_CT,filtered_nCT):
    mean_filtered_nCT = sum(filtered_nCT)/len(filtered_nCT)
    sum_upper = 0
    sum_lower = 0
    for i in range(len(filtered_X_CT)):
        yi = filtered_nCT[i]
        yi_pred = filtered_X_CT[i]*k+b
        y_mean = mean_filtered_nCT
        upper = (yi_pred-yi)**2
        lower = (yi-y_mean)**2
        sum_upper+= upper
        sum_lower+= lower
    R_2 = 1-sum_upper/sum_lower
    return R_2

X = [12+i*24 for i in range(9)]
FigureName_list=['CT_distribution_'+str(X[i])+'h' for i in range(len(X))]

for i in range(len(X)):
    print(len(X)-i)
    plotx = filtered_X_CT[int(X[i])]
    ploty = filtered_nCT[int(X[i])]
    k,b=optimize.curve_fit(f_1,plotx,ploty)[0]
    X_line = np.arange(0,600,0.01)
    Y_line = X_line*k+b
    R_2 = R2_func(plotx,ploty)
    FigureName = FigureName_list[i]
    fig, ax = plt.subplots(figsize=(2.33, 2.33),constrained_layout=False)
    ax.scatter([i+5 for i in plotx], ploty, s=3, color='red', alpha=0.8) #,label = 'shell')
    ax.plot(X_line,Y_line,color = 'black',linestyle = '-.')
#ax.axvline(x=linefitmin, color='grey', linestyle='--', linewidth=1)
    ax.set_title(str(X[i])+'h', fontsize = 8 ,loc = 'left')
    ax.set_ylabel( r'$\mathrm{ln(n(CT))}$', fontsize=7,labelpad = 0)
    ax.set_xlabel(r'$\mathrm{CT}$ (nm)', fontsize=7,labelpad = 0)
    ax.text(400,-1.5,'k = '+str(abs(round(k,3))),fontsize=7,horizontalalignment='left', verticalalignment='center')
    ax.text(400,-2.5,'$\mathrm{R^2}$ = '+str(abs(round(R_2,3))),fontsize=7,horizontalalignment='left', verticalalignment='center')
    ax.set_xlim([0,600])
    ax.set_ylim([-16,0])
    ax.set_yticks([-16,-12, -8,-4, 0])
#    ax.set_yticklabels([-10, -8,-6, -4, -2])
    ax.set_xticks([0,150,300,450,600])
    ax.set_xticklabels([0,150,300,450,600])
    ax.tick_params(axis='y', which='major', direction='in', width=0.5, length=2, labelsize=7, pad=2, rotation=0)
    ax.tick_params(axis='x', which='major', direction='in', width=0.5, length=2, labelsize=7, pad=2, rotation=0)
    ax.tick_params(top=True, bottom=True, left=True, right=True)
#ax.autoscale(tight=True)
#plt.legend(loc='upper right',fontsize=10,frameon=False)#,bbox_to_anchor = (1,0.5))
#fig.savefig(FigurePath + FigureName+'.pdf')
    fig.savefig(FigurePath+FigureName,dpi=1000)

