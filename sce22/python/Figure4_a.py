#Dp_Dc_distribution
#show Dc distribution and Dp distribution and k
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import csv
import os
import NC_X
from scipy import optimize
from matplotlib.ticker import MultipleLocator

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

Dp_Dc_Conc = []
#Read 
with open(DataPath+DataName, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        float_row = [float(x) for x in row]
        Dp_Dc_Conc.append(float_row)
print('validate the len of list')
print(len(Dp_Dc_Conc))
Dp = Dp_Dc_Conc[:240]
Dc = Dp_Dc_Conc[240:480]
conc = Dp_Dc_Conc[480:]
print(len(Dp))
print(len(Dc))
print(len(conc))
print()
File_number = len(Dp)


bin_gap=10
n_bin=60
#data cut
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
nCT10_40nm = [[] for i in range(File_number)] #save dp 60-80
nCT70_70nm = [[] for i in range(File_number)] #save dp 80-100
nCT70_100nm = [[] for i in range(File_number)] #save dp 100-120
nCT100_130nm = [[] for i in range(File_number)] #save dp 120-140

nCT = filter_fun(Dp,Dc,conc,bin_gap,n_bin,File_number,Dc_min=0,Dc_max=10000)
nCT10_40nm = filter_fun(Dp,Dc,conc,bin_gap,n_bin,File_number,Dc_min=50,Dc_max=60)
nCT40_70nm = filter_fun(Dp,Dc,conc,bin_gap,n_bin,File_number,Dc_min=60,Dc_max=70)
nCT70_100nm = filter_fun(Dp,Dc,conc,bin_gap,n_bin,File_number,Dc_min=70,Dc_max=80)
nCT100_130nm = filter_fun(Dp,Dc,conc,bin_gap,n_bin,File_number,Dc_min=80,Dc_max=90)


list_nCT = []
list_nCT10_40nm = []
list_nCT40_70nm = []
list_nCT70_100nm = []
list_nCT100_130nm = []
# time add
def nDp_to_filtered_nD(nDp,timeafter=48):
    timeafter=48
    nDp = list(map(list, zip(*nDp)))#do T process keep Dp-Bin as hang
    list_nDp = []
    for i in range(len(nDp)):
        list_nDp.append(sum(nDp[i][timeafter:])) # do Sum(n(dp[time]))
    BC_conc_SUM2 = sum(list_nDp) # do Sum n[dp]
    norm_nDp = [list_nDp[i]/BC_conc_SUM2 for i in range(len(list_nDp))] #do normalized process
    return norm_nDp
mean_CT_bin = 0
norm_nCT = nDp_to_filtered_nD(nCT,timeafter=48)

norm_nCT10_40nm = nDp_to_filtered_nD(nCT10_40nm,timeafter=48)
norm_nCT40_70nm = nDp_to_filtered_nD(nCT40_70nm,timeafter=48)
norm_nCT70_100nm = nDp_to_filtered_nD(nCT70_100nm,timeafter=48)
norm_nCT100_130nm = nDp_to_filtered_nD(nCT100_130nm,timeafter=48)

print('validation of norm_sum=1?')
print(sum(norm_nCT)) #validation
print()

nCT_ln = list(map(logF, norm_nCT)) #log process
nCT_ln10_40nm = list(map(logF, norm_nCT10_40nm))
nCT_ln40_70nm = list(map(logF, norm_nCT40_70nm))
nCT_ln70_100nm = list(map(logF, norm_nCT70_100nm))
nCT_ln100_130nm = list(map(logF, norm_nCT100_130nm))

linefitmin = 0
filtered_X_CT = [i for i, j in zip(x,nCT_ln) if j != 'a' and i>=linefitmin ] # do filter n>0 and Dp>=140nm
filtered_nCT = [j for i, j in zip(x,nCT_ln) if j != 'a' and i>=linefitmin]

filtered_X_CT10_40nm = [i for i, j in zip(x,nCT_ln10_40nm) if j != 'a' ] # do filter n>0 and Dp>=140nm
filtered_nCT10_40nm = [j for i, j in zip(x,nCT_ln10_40nm) if j != 'a' ]

filtered_X_CT40_70nm = [i for i, j in zip(x,nCT_ln40_70nm) if j != 'a' ] # do filter n>0 and Dp>=140nm
filtered_nCT40_70nm = [j for i, j in zip(x,nCT_ln40_70nm) if j != 'a' ]

filtered_X_CT70_100nm = [i for i, j in zip(x,nCT_ln70_100nm) if j != 'a' ] # do filter n>0 and Dp>=140nm
filtered_nCT70_100nm = [j for i, j in zip(x,nCT_ln70_100nm) if j != 'a' ]

filtered_X_CT100_130nm = [i for i, j in zip(x,nCT_ln100_130nm) if j != 'a' ] # do filter n>0 and Dp>=140nm
filtered_nCT100_130nm = [j for i, j in zip(x,nCT_ln100_130nm) if j != 'a' ]

filtered_X_CT10_40nm=[i+5 for i in filtered_X_CT10_40nm]
filtered_X_CT40_70nm=[i+5 for i in filtered_X_CT40_70nm]
filtered_X_CT70_100nm=[i+5 for i in filtered_X_CT70_100nm]
filtered_X_CT100_130nm=[i+5 for i in filtered_X_CT100_130nm]
filtered_X_CT=[i+5 for i in filtered_X_CT]



#coculate B and K
k,b=optimize.curve_fit(f_1,filtered_X_CT,filtered_nCT)[0]
k10_40nm,b10_40nm=optimize.curve_fit(f_1,filtered_X_CT10_40nm,filtered_nCT10_40nm)[0]
k40_70nm,b40_70nm=optimize.curve_fit(f_1,filtered_X_CT40_70nm,filtered_nCT40_70nm)[0]
k70_100nm,b70_100nm=optimize.curve_fit(f_1,filtered_X_CT70_100nm,filtered_nCT70_100nm)[0]
k100_130nm,b100_130nm=optimize.curve_fit(f_1,filtered_X_CT100_130nm,filtered_nCT100_130nm)[0]
#data: y=kx+b
K_list = [k10_40nm,k40_70nm,k70_100nm,k100_130nm]

print('sep_bin k')
print(K_list)
print(max(K_list))
print(min(K_list))
print()
#x1 = np.arange(0,600,0.01)
x1 = np.arange(min(filtered_X_CT),max(filtered_X_CT),0.01)
y1 = x1*k+b
print('1/k:',end='')
print(1/k)
print()
x10_40nm = np.arange(min(filtered_X_CT10_40nm),max(filtered_X_CT10_40nm),0.01)
x40_70nm = np.arange(min(filtered_X_CT40_70nm),max(filtered_X_CT40_70nm),0.01)
x70_100nm = np.arange(min(filtered_X_CT70_100nm),max(filtered_X_CT70_100nm),0.01)
x100_130nm = np.arange(min(filtered_X_CT100_130nm),max(filtered_X_CT100_130nm),0.01)

y10_40nm = x10_40nm*k10_40nm+b10_40nm
y40_70nm = x40_70nm*k40_70nm+b40_70nm
y70_100nm = x70_100nm*k70_100nm+b70_100nm
y100_130nm = x100_130nm*k100_130nm+b100_130nm

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
print('R2:',end='')
print(R_2)
print()
#plot1
FigureName = 'Figure4_a'
fig, ax = plt.subplots(figsize=(3, 2.9),constrained_layout=False)
fig.subplots_adjust(left=0.15)
fig.subplots_adjust(bottom=0.15)
#fig.subplots_adjust(top=0.90)
#fig.subplots_adjust(right=0.90)
ax.scatter(filtered_X_CT, filtered_nCT, s=5, color='red', alpha=0.8) #,label = 'shell')
ax.plot(x1,y1,color = 'black',linestyle = '-.')
#ax.axvline(x=linefitmin, color='grey', linestyle='--', linewidth=1)
ax.set_title('(a)', fontsize = 12 ,loc = 'left')
ax.set_ylabel( r'$\mathrm{ln(n(CT))}$', fontsize=10,labelpad = 0)
ax.set_xlabel(r'$\mathrm{CT}$ (nm)', fontsize=10,labelpad = 5)
ax.text(400,-1.5,'k = '+str(abs(round(k,3))),fontsize=9,horizontalalignment='left', verticalalignment='center')
ax.text(400,-2.5,'$\mathrm{R^2}$ = '+str(abs(round(R_2,3))),fontsize=9,horizontalalignment='left', verticalalignment='center')
ax.set_ylim([-14,0])
ax.set_yticks([-12,-8,-4,0])
ax.set_xlim([0,600])
ax.set_xticks([0,150,300,450,600])
ax.set_xticklabels([0,150,300,450,600])
ax.tick_params(axis='y', which='major', direction='in', width=0.8, length=4, labelsize=10, pad=3, rotation=0)
ax.tick_params(axis='x', which='major', direction='in', width=0.8, length=4, labelsize=10, pad=3, rotation=0)
ax.tick_params(top=True, bottom=True, left=True, right=True)
fig.savefig(FigurePath+FigureName,dpi=1000)


