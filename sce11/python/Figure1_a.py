from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import os
import NC_X
import matplotlib
import numpy as np
import move as mv
import matplotlib.font_manager as font_manager
from scipy.optimize import curve_fit 

font_manager.fontManager.addfont('/data/home/zzy/.fonts/ARIAL.TTF')
plt.rcParams['font.sans-serif']= 'ARIAL'

DataPath = '../out/'
FigurePath = '../Figure/'

lists = [] #savd filenames
BC_containing_conc = []
BC_containing_mass = []
Y_Mix_X = []
Y_mean_DcDp = []
X = []
for file in os.listdir(DataPath):
        #print(file)
   if os.path.splitext(file)[1].lower() in '.nc':
            lists.append(file)
lists = sorted(lists,key=lambda keys:[ord(i) for i in keys],reverse=False)
lists = lists[1:]

for i in range(len(lists)):
    DataName = lists[i]
    Y_Mix_X.append(NC_X.coculateC_Opt_X_func(DataPath,DataName))
    X.append(int(DataName[-6:-3]))
X = [i-1 for i in X]
#do move average
Y_Mix_Xmove = mv.smoothmoveavg(Y_Mix_X,24,'valid')
X1 = mv.smoothmoveavg(X,24,'valid')
X1 = [i-0.5 for i in X1]


#fitting
Y0 = Y_Mix_X[0]
Yinf = sum(Y_Mix_X[216:])/len(Y_Mix_X[216:])
def function(t,K):
    y = Y0*np.exp(-K*t)+Yinf*(1-np.exp(-K*t))
    return y
x_fit = np.linspace(1,240,2400)

#Y_init =np.array([(i - Yinf) for i in Y_Mix_X])

X_init = np.array(X)
Y_init = np.array(Y_Mix_X)
p_est_l,err_est_l = curve_fit(function,X_init,Y_init,method='lm')
print(p_est_l)
print(Y0)
print(Yinf)
y_fit = function(x_fit,*p_est_l)
tau = 1/p_est_l[-1]
taulist = [tau*(1+i) for i in range(1)]
taulist = np.array(taulist)
tauYlist = function(taulist,*p_est_l)
#print(T_steady)

X_bias = X
Y_bias = [abs(Y_Mix_X[i]-function(X_bias[i],*p_est_l)) for i in range(len(X_bias))]

#draw code
FigureName = 'Figure1_aopt'
fig, ax = plt.subplots(figsize=(7, 3),constrained_layout=False)
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.1)
fig.subplots_adjust(right=0.95)
ax.plot(X, Y_Mix_X, color='#1c1464',label = '$\mathrm{\chi}$')
ax.plot(X1, Y_Mix_Xmove, color='#1c1464',linestyle = '--',alpha=0.6,label = 'moving averaged $\mathrm{\chi}$')
ax.plot(x_fit,[abs(i-Yinf) for i in y_fit], color='red', linestyle='-.',label = r'$\mathrm{\chi_{fit}}$ - $\mathrm{\chi_{\infty}}$')
ax.scatter(tau,abs(tauYlist[0]-Yinf),marker='*',color = 'red',s=10)
ax.axvline(x=tau, ymin=0, ymax=1, color='black', linestyle='--', linewidth=1)
ax.text(tau+13,0.25,'$\mathrm{\\tau}}$ = '+str(round(tau,1))+'h',fontsize=10,horizontalalignment='center', verticalalignment='center')
ax.legend(loc='upper right',fontsize=9,frameon=False)#,bbox_to_anchor = (1,0.5))
ax.set_title(r'(a)', fontsize = 12 ,loc = 'left')
ax.set_xlabel('Time of simulation(h)', fontsize=10)
ax.set_ylabel(r'Mixing state metric $\mathrm{\chi}$(%)', fontsize=10)
ax.set_xlim([-5,245])
ax.set_ylim([0,1])
ax.set_yticks([0, 0.25,0.5,0.75,1])
ax.set_yticklabels([0,25,50,75,100])
ax.set_xticks([24*i for i in range(11)])
ax.set_xticklabels([24*i for i in range(11)])
ax.tick_params(axis='y', which='major', direction='in', width=0.8, length=4, labelsize=10, pad=2, rotation=0)
ax.tick_params(axis='x', which='major', direction='in', width=0.8, length=4, labelsize=10, pad=2, rotation=0)
ax.tick_params(top=True, bottom=True, left=True, right=True)
fig.savefig(FigurePath+FigureName,dpi=1000)
