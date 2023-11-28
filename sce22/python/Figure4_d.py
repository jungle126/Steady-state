#Dp_Dc_distribution
#show Dc distribution and Dp distribution and k 
import csv
import math
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os
import NC_X
from scipy import optimize
import matplotlib.font_manager as font_manager
font_manager.fontManager.addfont('/data/home/zzy/.fonts/ARIAL.TTF')
plt.rcParams['font.sans-serif']= 'ARIAL'


#log10 func
def log10F(va):
    if va == 0:
        return 0
    else:
        return math.log10(va)

#set DataPath
DataPath = '../out/'
#set FigsavePath
FigurePath = '../Figure/'
DataName = 'urban_plume_0001_00000013.nc'
#set nf name
nf = nc.Dataset(DataPath + DataName,'r')
aero_conc = np.array(nf.variables['aero_num_conc'])
mass = np.array(nf.variables['aero_particle_mass'])
#mass=20*5000
conc = aero_conc.tolist()
BC_mass = mass[18] # BC mass
BC_mass = BC_mass.tolist() #np to list
POA_mass = mass[17]  #POA(OC) mass
POA_mass = POA_mass.tolist()
SOA_mass = mass[5]+mass[6]+mass[7]+mass[8]+mass[9]+mass[10]+mass[11]+mass[12] 
SOA_mass = SOA_mass.tolist()
SO4_mass = mass[0]
SO4_mass = SO4_mass.tolist()
NH4_mass = mass[3]
NH4_mass = NH4_mass.tolist()
NO3_mass = mass[1]
NO3_mass = NO3_mass.tolist()

massT = mass.T #mass = 5000*20
P_mass = [] #mass of particle
for i in range(len(massT)):
    P_mass.append(sum(massT[i]))

##filter the BC-containing particle
#new list to save BC-containing paticle
conc_BC = []
BC_BC = []
POA_BC = []
SOA_BC = []
SO4_BC = []
NH4_BC = []
NO3_BC = []
pg_BC = []
for i in range(len(BC_mass)):
    if BC_mass[i]!= 0:
        conc_BC.append(conc[i])
        BC_BC.append(BC_mass[i])
        POA_BC.append(POA_mass[i])
        SOA_BC.append(SOA_mass[i])
        SO4_BC.append(SO4_mass[i])
        NH4_BC.append(NH4_mass[i])
        NO3_BC.append(NO3_mass[i])
        pg_BC.append(BC_mass[i]*1e15)
#test sum_mass ug/m3
BC_SUM = sum([conc_BC[i]*BC_BC[i]*1e9 for i in range(len(conc_BC))])
POA_SUM = sum([conc_BC[i]*POA_BC[i]*1e9 for i in range(len(conc_BC))])
SOA_SUM = sum([conc_BC[i]*SOA_BC[i]*1e9 for i in range(len(conc_BC))])
SO4_SUM = sum([conc_BC[i]*SO4_BC[i]*1e9 for i in range(len(conc_BC))])
print(BC_SUM)
print(POA_SUM)
print(SOA_SUM)
print(SO4_SUM)



#for i in range(len(BC_mass)):
pg_log = list(map(log10F, pg_BC)) # x-pg to x-logpg
#log -5to-1 0.2bin  data cut
binlen = 20
BC = [0 for i in range(binlen)]
POA = [0 for i in range(binlen)]
SOA = [0 for i in range(binlen)]
SO4 = [0 for i in range(binlen)]
NH4 = [0 for i in range(binlen)]
NO3 = [0 for i in range(binlen)]
X_pg = [-5+4/binlen*i for i in range(binlen)]
for i in range(len(pg_log)):
    for j in range(binlen):
        pg_log_min = -5+4/binlen*j 
        pg_log_max = -5+4/binlen*(j+1)
        if pg_log_min <= pg_log[i] < pg_log_max:
            BC[j] += conc_BC[i]*BC_BC[i]*1e9   #1e9 is kg/m3 to ug/m3
            POA[j] += conc_BC[i]*POA_BC[i]*1e9
            SOA[j] += conc_BC[i]*SOA_BC[i]*1e9
            SO4[j] += conc_BC[i]*SO4_BC[i]*1e9
            NH4[j] += conc_BC[i]*NH4_BC[i]*1e9
            NO3[j] += conc_BC[i]*NO3_BC[i]*1e9
VF_BC = []
for i in range(len(X_pg)):
    V_BC = BC[i]/1800
    V_POA = POA[i]/1000
    V_SOA = SOA[i]/1400
    V_SO4 = SO4[i]/1800
    V_NO3 = NO3[i]/1800
    V_NH4 = NH4[i]/1800
    if (V_BC+V_POA+V_SOA+V_SO4+V_NO3+V_NH4) == 0:
        VF_BC_i = 0
    else:VF_BC_i = V_BC/(V_BC+V_POA+V_SOA+V_SO4+V_NO3+V_NH4)
    VF_BC.append(VF_BC_i)

SavePath = '../Data/'
CsvName = 'VF_massBC.csv'
datalist = [VF_BC,[10**i for i in X_pg]]
#save csv file
with open(SavePath+CsvName, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    for lst in datalist:
        csvwriter.writerow(lst)
T1 = max([BC[i]+POA[i]+SOA[i]+SO4[i]+NH4[i]+NO3[i] for i in range(len(X_pg))])
BC = [i/T1 for i in BC]
POA = [i/T1 for i in POA]
SOA = [i/T1 for i in SOA]
SO4 = [i/T1 for i in SO4]
NH4 = [i/T1 for i in NH4]
NO3 = [i/T1 for i in NO3]

MF_POA = sum(POA)/sum(BC) 
MF_SOA = sum(SOA)/sum(BC)
MF_SO4 = sum(SO4)/sum(BC)
MF_NH4 = sum(NH4)/sum(BC)
MF_NO3 = sum(NO3)/sum(BC)

POA_U = [i*MF_POA for i in BC]
SOA_U = [i*MF_SOA for i in BC]
SO4_U = [i*MF_SO4 for i in BC]
NH4_U = [i*MF_NH4 for i in BC]
NO3_U = [i*MF_NO3 for i in BC]

T1 = max([BC[i]+POA[i]+SOA[i]+SO4[i]+NH4[i]+NO3[i] for i in range(len(X_pg))])
BC = [i/T1 for i in BC]

FigureName = 'Figure4_d.png'
fig, ax = plt.subplots(figsize=(7, 3.25),constrained_layout=False) #7;3.25
#fig.subplots_adjust(top=0.9)
#fig.subplots_adjust(bottom=0.125)
#fig.subplots_adjust(left=0.15)
#fig.subplots_adjust(right=0.95

color_list = ['#281A90','#0B53DB','#1788C4','#32AA88','#95B457','#FBB232']
color_list_U = [color + '80' for color in color_list]
label_list = ['BC','POA','SOA','$\mathrm{SO_4}$','$\mathrm{NH_4}$','$\mathrm{NO_3}$']
ax.tick_params(top='on', right='on', which='both')
ax.stackplot(X_pg, BC, POA_U, SOA_U,SO4_U,NH4_U,NO3_U, colors=color_list_U)
ax.stackplot(X_pg, BC, POA, SOA,SO4,NH4,NO3, colors=color_list, labels=label_list)
ax.plot(X_pg, BC, linestyle = '--',color = 'white')
ax.plot(X_pg, [BC[i]+POA[i]+SOA[i]+SO4[i]+NH4[i]+NO3[i] for i in range(len(X_pg))], linestyle = '-',color = 'black')
ax.plot(X_pg, [BC[i]+POA_U[i]+SOA_U[i]+SO4_U[i]+NH4_U[i]+NO3_U[i] for i in range(len(X_pg))], linestyle = '-',color = 'grey')
plt.legend(loc='upper right',fontsize=9,frameon=False)
ax.set_title('(d)',fontsize = 12,loc = 'left')
plt.ylim(0,1.2)
ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
ax.set_xlim(-5,-2-4/binlen)
ax.set_xticks([-4,-3,-2])
ax.set_xticklabels(['1$0^{-4}$','1$0^{-3}$','1$0^{-2}$'])
ax.set_xlabel('Mass of BC contained in each particle(pg)',fontsize = 10,labelpad = 2)
#ax.set_ylabel('Distribution of aerosol components with\nrespect to per-particle BC mass($\mathrm{\mu g/m^3}$)',fontsize =10,labelpad = 2)
ax.set_ylabel('Distribution of aerosol components with\nrespect to per-particle BC mass',fontsize =10,labelpad = 2)
ax.tick_params(axis='y', which='major', direction='in', width=0.8, length=4, labelsize=10, pad=2, rotation=0)
ax.tick_params(axis='x', which='major', direction='in', width=0.8, length=4, labelsize=10, pad=3, rotation=0)
fig.savefig(FigurePath+FigureName,dpi=1000)
exit()




