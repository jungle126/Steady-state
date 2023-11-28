import numpy as np
import matplotlib.pyplot as plt
import os
import NC_X as NC
import matplotlib.colors as mcolors
DataPath = '../out/'
FigurePath = '../Figure/'
FigureName = 'aero_mass_pie'
a = 'SO4,NO3,Cl,NH4,MSA,ARO1,ARO2,ALK1,OLE1,APl1,APl2,LIM1,LIM2,CO3,Na,Ca,OIN,OC,BC,H2O'
aero_list = a.split(',')
#print(gas_list)
#exit()

lists = []
X = []

for file in os.listdir(DataPath):
        #print(file)
   if os.path.splitext(file)[1].lower() in '.nc':
            lists.append(file)
lists = sorted(lists,key=lambda keys:[ord(i) for i in keys],reverse=False)
lists = lists[1:]

SO4 = []
NO3 = []
Cl = []
NH4 = []
MSA = []
ARO1 = []
ARO2 = []
ALK1 = []
OLE1 = []
API1 = []
API2 = []
LIM1 = []
LIM2 = []
CO3 = []
Na = []
Ca = []
OIN = []
OC = []
BC = []
H2O = []
Aero = []

X = [] #time label
for i in range(len(lists)):
    DataName = lists[i]
    SO4.append(NC.coculateC_SO4_masssum_func(DataPath,DataName))    #0
    NO3.append(NC.coculateC_NO3_masssum_func(DataPath,DataName))    #1
    Cl.append(NC.coculateC_Cl_masssum_func(DataPath,DataName))      #2
    NH4.append(NC.coculateC_NH4_masssum_func(DataPath,DataName))    #3
    MSA.append(NC.coculateC_MSA_masssum_func(DataPath,DataName))    #4
    ARO1.append(NC.coculateC_ARO1_masssum_func(DataPath,DataName))  #5
    ARO2.append(NC.coculateC_ARO2_masssum_func(DataPath,DataName))  #6
    ALK1.append(NC.coculateC_ALK1_masssum_func(DataPath,DataName))  #7
    OLE1.append(NC.coculateC_OLE1_masssum_func(DataPath,DataName))  #8
    API1.append(NC.coculateC_API1_masssum_func(DataPath,DataName))  #9
    API2.append(NC.coculateC_API2_masssum_func(DataPath,DataName))  #10
    LIM1.append(NC.coculateC_LIM1_masssum_func(DataPath,DataName))  #11
    LIM2.append(NC.coculateC_LIM2_masssum_func(DataPath,DataName))  #12
    CO3.append(NC.coculateC_CO3_masssum_func(DataPath,DataName))    #13
    Na.append(NC.coculateC_Na_masssum_func(DataPath,DataName))      #14
    Ca.append(NC.coculateC_Ca_masssum_func(DataPath,DataName))      #15
    OIN.append(NC.coculateC_OIN_masssum_func(DataPath,DataName))    #16
    OC.append(NC.coculateC_OC_masssum_func(DataPath,DataName))      #17
    BC.append(NC.coculateC_BC_masssum_func(DataPath,DataName))      #18
    H2O.append(NC.coculateC_H20_masssum_func(DataPath,DataName))    #19
    Aero.append(NC.coculate_particle_summass_func(DataPath,DataName)) #sum
    X.append(100*eval(DataName[-6])+10*eval(DataName[-5])+eval(DataName[-4])-1)

#OM
OC_avg = sum(OC)/len(OC)*1e9
MSA_avg = sum(MSA)/len(MSA)*1e9
ARO1_avg = sum(ARO1)/len(ARO1)*1e9
ARO2_avg = sum(ARO2)/len(ARO2)*1e9
ALK1_avg = sum(ALK1)/len(ALK1)*1e9
OLE1_avg = sum(OLE1)/len(OLE1)*1e9
API1_avg = sum(API1)/len(API1)*1e9
API2_avg = sum(API2)/len(API2)*1e9
LIM1_avg = sum(LIM1)/len(LIM1)*1e9
LIM2_avg = sum(LIM2)/len(LIM2)*1e9
OIN_avg = sum(OIN)/len(OIN)*1e9
OM_avg = OC_avg+MSA_avg+ARO1_avg+ARO2_avg+ALK1_avg+OLE1_avg+API1_avg+API2_avg+LIM1_avg+LIM2_avg+OIN_avg

#BC SO4 NO3 NH3
BC_avg = sum(BC)/len(BC)*1e9
SO4_avg = sum(SO4)/len(SO4)*1e9
NO3_avg = sum(NO3)/len(NO3)*1e9
NH4_avg = sum(NH4)/len(NH4)*1e9

#aero_Sum
Aero_avg =sum(Aero)/len(Aero)*1e9
other_avg = Aero_avg-OM_avg-BC_avg-SO4_avg-NO3_avg-NH4_avg

##test
#print(Aero_avg)
#print(BC_avg)
#print(OM_avg)
#print(SO4_avg)
#print(NO3_avg)
#print(NH4_avg)
#print(other_avg)
#exit()

# data process
pie_data = [OM_avg,BC_avg,NO3_avg,SO4_avg,NH4_avg]
name_list1 = ['OM','BC','NO3','SO4','NH4']
name_list = ['OM avg','BC avg','NO$_3^-$ avg','SO$_4^{2-}$ avg','NH$_4^+$ avg']
Sum_pie=sum(pie_data)
for i in range(len(pie_data)):
    print(name_list1[i],end=':')
    print("%.1f"%(pie_data[i]/Sum_pie*100),end = '%')
    print()
exit()
col_list = ['#A1A9D0','#73bad6','#F0988C','#F6CAE5','#96CCCB']

#plt
plt.figure(figsize=(15,10),dpi=300)
patches, texts, autotexts=plt.pie(x=pie_data,labels=name_list,colors=col_list,center=(0,3),autopct='%.2f%%',labeldistance=None,pctdistance=0.6,textprops={'fontsize':15,'color':'black'})
plt.legend(patches,name_list,fontsize='15',loc='lower center',frameon=False,ncol=6)#,bbox_to_anchor=(0,0,0.5,1))
plt.savefig(FigurePath+FigureName)
