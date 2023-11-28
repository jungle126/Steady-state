# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 15:12:13 2021

@author: Jungle
"""
#flog
#modify 2022-4-12 !fix the problem about BC_conc_sum
#modify 2022-4-14 !add gas_background variables in 1.4
#modify 2022-5-11 !fix the problem about BC_conc_sum again(with conc)
#modify 2022-5-14 !fix the problem about Opt_X(with)
#modify 2022-7-10 !fix the problem about X>1 (log2->log)
#modify 2022-7-10 !fix the problem about X>1 (opt_X)
#updata 2023-5-22 !improve the coculation of Dp Dc in 1.1.3 (core-shell mode)
#updata 2023-5-26 !change the mean_Dc_Dp(all particles) to mean_Dc_Dp(BC-containing)in coculateC_DcDp_func
#updata 2023-5-29 !fix the bug meandc_dp of a time = 0 
#plan
#BC_conc_sum(without conc)  and  Opt_X(without conc) need to repair

import math
import numpy as np
import netCDF4 as nc
# %% 1  function used to process the NCfile with 'aero_num_conc'

# %%% 1.1  function used to coculate the value about mix state


# %%%% 1.1.1  function used to coculate the X of particle (contain conc value)
# input DataPath,DataName
# output X
def coculateC_X_func(DataPath,DataName):
    """coculate the X of particle (contain conc value) 
    input:DataPath,DataName
    output:X"""
    nf = nc.Dataset(DataPath + DataName,'r')
    aero_conc = np.array(nf.variables['aero_num_conc'])
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    conc = aero_conc.tolist()
    Ha = 0 
    U = 0  
    for i in range(len(mass)):
        U+=sum(mass[i])*conc[i]
    for i in range(len(mass)):
        Pi = sum(mass[i])/U
        Hi = 0
        for a in range(len(mass_T)):
            if mass[i][a]!=0:
                Pia = mass[i][a]/sum(mass[i])
                Hi -= Pia*math.log(Pia)
        Ha += Hi*conc[i]*Pi
    Da = math.exp(Ha)   
    Hr = 0
    for a in range(len(mass_T)):
        Ua = 0
        for i in range(len(mass)):
            Ua += mass[i][a]*conc[i]
        Pa = Ua/U
        if Pa != 0:
            Hr -= Pa*math.log(Pa)
    Dr = math.exp(Hr)
    if Da ==1 or Dr ==1:
        return 0
    else:
        X = (Da-1)/(Dr-1)
        return X

# %%%%  1.1.2  function used to coculate the optical_X of particle (contain conc value)
# input DataPath,DataName
# output Opt_X
def coculateC_Opt_X_func(DataPath,DataName):
    """coculate the optical_X of particle (contain conc value)
    input:DataPath,DataName
    output:Opy_X"""
    nf = nc.Dataset(DataPath + DataName,'r')
    aero_conc = np.array(nf.variables['aero_num_conc'])
    mass = np.array(nf.variables['aero_particle_mass'])
    mass = mass.T
    mass = mass.tolist()
    conc = aero_conc.tolist()

    Opt_Ha = 0
    Opt_U = 0
    Opt_U_BC = 0
    for i in range(len(mass)):
        Opt_U+=sum(mass[i])*conc[i]

    for i in range(len(mass)):
        Opt_Pi = sum(mass[i])/Opt_U
        Opt_Pi_BC = mass[i][18]/sum(mass[i])
        Opt_Pi_NBC = 1-Opt_Pi_BC
        Opt_Hi = 0
        for j in [Opt_Pi_NBC,Opt_Pi_BC]:
            if j != 0:
                Opt_Hi = Opt_Hi-j*math.log(j)
        #modify here
        Opt_Ha += Opt_Hi*conc[i]*Opt_Pi
        Opt_U_BC+=mass[i][18]*conc[i]
    Opt_P_BC = Opt_U_BC/Opt_U
    Opt_Da = math.exp(Opt_Ha)
    Opt_P_NBC = 1-Opt_P_BC
    Opt_Hr = 0
    for k in [Opt_P_NBC,Opt_P_BC]:
        if k != 0:
            Opt_Hr -= k*math.log(k)
    Opt_Dr = math.exp(Opt_Hr)
    if Opt_Da == 1 or Opt_Dr == 1:
        return 0
    else:
        Opt_X = (Opt_Da-1)/(Opt_Dr-1)
        return Opt_X

# %%%% 1.1.3  function used to coculate the values of diameter of particle (contain conc value)
# input DataPath,DataName
# output Dc_Dp_mean,Dc,Dp,Dc_Dp
def coculateC_DcDp_func(DataPath,DataName):
    """coculate the values of diameter of particle (contain conc value) 
    input:DataPath,DataName
    output:Dc_Dp_mean,Dc,Dp,Dc_Dp"""
    nf = nc.Dataset(DataPath + DataName,'r')
    aero_den = np.array(nf.variables['aero_density'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    mass = np.array(nf.variables['aero_particle_mass'])
    pi = 3.14159
    den = aero_den.tolist()
    conc = aero_conc.tolist()
    BC_mass = mass[18] # BC质量
    BC_mass = BC_mass.tolist()
    
    
    P_mass = [] #单个particle总质量
    rho_shell_list = [] # particle密度
    Dc_list = [] # BC核直径
    Dp_list = [] # particle的直径
    Dc_Dp_list = []
    mass = mass.T
    rho_core = den[18]
    for i in range(len(BC_mass)):
        P_mass.append(sum(mass[i]))
        rho_shell = 0
        for j in range(len(den)):
            if j!=18:
                rho_shell += den[j]*mass[i,j]/(sum(mass[i])-BC_mass[i])
        rho_shell_list.append(rho_shell)
        Dc_m = math.pow(BC_mass[i]/den[18]*6/pi,1/3)
        Dc_nm = Dc_m*1e9
        if rho_shell==0: #BC_mass[i]==P_mass[i]:
            Dp_nm=Dc_nm
        else:
            Dp_m = math.pow((6*P_mass[i]-pi*Dc_m**3*(rho_core-rho_shell))/(pi*rho_shell),1/3)
            Dp_nm = Dp_m*1e9
        Dc_list.append(Dc_nm)
        Dp_list.append(Dp_nm)
        Dc_Dp_list.append(Dc_nm/Dp_nm)
    Dc_Dp_sum = 0
    conc_sum = 0
    for i in range(len(Dc_Dp_list)):
        if Dc_Dp_list[i] != 0:
            Dc_Dp_sum += Dc_Dp_list[i]*conc[i]
            conc_sum += conc[i]
    Dc_Dp_mean = Dc_Dp_sum/conc_sum
    return Dc_Dp_mean,Dc_list,Dp_list,Dc_Dp_list

# %%% 1.2 function used to coculate the values of mass
# %%%% 1.2.1  function used to coculate the sum_mass of SO4
def coculateC_SO4_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of SO4 (contain conc value)
    input:DataPath,DataName
    output:SO4_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    SO4 = mass_T[0]
    func = lambda x,y:x*y
    SO4_sum = sum(list(map(func,SO4,conc))) 
    return SO4_sum
# %%%% 1.2.2  function used to coculate the list_mass of SO4
def coculate_SO4_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of SO4
    input:DataPath,DataName
    output:SO4"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    SO4 = mass_T[0]
    return SO4
# %%%% 1.2.3  function used to coculate the sum_mass of NO3
def coculateC_NO3_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of NO3 (contain conc value)
    input:DataPath,DataName
    output:NO3_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass_T = mass_T.tolist()
    mass = mass.tolist()
    NO3 = mass_T[1]
    func = lambda x,y:x*y
    NO3_sum = sum(list(map(func,NO3,conc)))        
    return NO3_sum
# %%%% 1.2.4  function used to coculate the list_mass of NO3
def coculate_NO3_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of NO3
    input:DataPath,DataName
    output:NO3"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass_T = mass_T.tolist()
    mass = mass.tolist()
    NO3 = mass_T[1]
    return NO3
# %%%% 1.2.5  function used to coculate the sum_mass of Cl
def coculateC_Cl_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of Cl (contain conc value)
    input:DataPath,DataName
    output:Cl_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Cl = mass_T[2]
    func = lambda x,y:x*y
    Cl_sum = sum(list(map(func,Cl,conc))) 
    return Cl_sum
# %%%% 1.2.6  function used to coculate the list_mass of Cl
def coculate_Cl_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of Cl
    input:DataPath,DataName
    output:Cl"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Cl = mass_T[2]
    return Cl
# %%%% 1.2.7  function used to coculate the sum_mass of NH4
def coculateC_NH4_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of NH4 (contain conc value)
    input:DataPath,DataName
    output:NH4_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    NH4 = mass_T[3]
    func = lambda x,y:x*y
    NH4_sum = sum(list(map(func,NH4,conc))) 
    return NH4_sum
# %%%% 1.2.8  function used to coculate the sum_mass of NH4
def coculate_NH4_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of NH4
    input:DataPath,DataName
    output:NH4"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    NH4 = mass_T[3]
    return NH4
# %%%% 1.2.9  function used to coculate the sum_mass of MSA
def coculateC_MSA_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of MSA (contain conc value)
    input:DataPath,DataName
    output:MSA_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    MSA = mass_T[4]
    func = lambda x,y:x*y
    MSA_sum = sum(list(map(func,MSA,conc))) 
    return MSA_sum
# %%%% 1.2.10  function used to coculate the list_mass of MSA
def coculate_MSA_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of MSA
    input:DataPath,DataName
    output:MSA"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    MSA = mass_T[4]
    return MSA
# %%%% 1.2.11  function used to coculate the sum_mass of ARO1
def coculateC_ARO1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of ARO1 (contain conc value)
    input:DataPath,DataName
    output:ARO1_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ARO1 = mass_T[5]
    func = lambda x,y:x*y
    ARO1_sum = sum(list(map(func,ARO1,conc))) 
    return ARO1_sum
# %%%% 1.2.12  function used to coculate the list_mass of ARO1
def coculate_ARO1_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of ARO1
    input:DataPath,DataName
    output:ARO1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ARO1 = mass_T[5]
    return ARO1
# %%%% 1.2.13  function used to coculate the sum_mass of ARO2
def coculateC_ARO2_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of ARO2 (contain conc value)
    input:DataPath,DataName
    output:ARO2_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ARO2 = mass_T[6]
    func = lambda x,y:x*y
    ARO2_sum = sum(list(map(func,ARO2,conc))) 
    return ARO2_sum
# %%%% 1.2.14  function used to coculate the list_mass of ARO2
def coculate_ARO2_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of ARO2
    input:DataPath,DataName
    output:ARO2"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ARO2 = mass_T[6]
    return ARO2
# %%%% 1.2.15  function used to coculate the sum_mass of ALK1
def coculateC_ALK1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of ALK1 (contain conc value)
    input:DataPath,DataName
    output:ALK1_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ALK1 = mass_T[7]
    func = lambda x,y:x*y
    ALK1_sum = sum(list(map(func,ALK1,conc))) 
    return ALK1_sum
# %%%% 1.2.16  function used to coculate the list_mass of ALK1
def coculate_ALK1_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of ALK1
    input:DataPath,DataName
    output:ALK1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ALK1 = mass_T[7]
    return ALK1
# %%%% 1.2.17 function used to coculate the sum_mass of OLE1
def coculateC_OLE1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of OLE1 (contain conc value)
    input:DataPath,DataName
    output:OLE1_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OLE1 = mass_T[8]
    func = lambda x,y:x*y
    OLE1_sum = sum(list(map(func,OLE1,conc))) 
    return OLE1_sum
# %%%% 1.2.18  function used to coculate the list_mass of OLE1
def coculate_OLE1_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of OLE1
    input:DataPath,DataName
    output:OLE1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OLE1 = mass_T[8]
    return OLE1
# %%%% 1.2.19  function used to coculate the sum_mass of API1
def coculateC_API1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of API1 (contain conc value)
    input:DataPath,DataName
    output:API1_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    API1 = mass_T[9]
    func = lambda x,y:x*y
    API1_sum = sum(list(map(func,API1,conc))) 
    return API1_sum
# %%%% 1.2.20  function used to coculate the list_mass of API1
def coculate_API1_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of API1
    input:DataPath,DataName
    output:API1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    API1 = mass_T[9]
    return API1
# %%%% 1.2.21  function used to coculate the sum_mass of API2
def coculateC_API2_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of API2 (contain conc value)
    input:DataPath,DataName
    output:API2_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    API2 = mass_T[10]
    func = lambda x,y:x*y
    API2_sum = sum(list(map(func,API2,conc))) 
    return API2_sum
# %%%% 1.2.22  function used to coculate the list_mass of API2
def coculate_API2_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of API2
    input:DataPath,DataName
    output:API2"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    API2 = mass_T[10]
    return API2
# %%%% 1.2.23  function used to coculate the sum_mass of LIM1
def coculateC_LIM1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of LIM1 (contain conc value)
    input:DataPath,DataName
    output:LIM1_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    LIM1 = mass_T[11]
    func = lambda x,y:x*y
    LIM1_sum = sum(list(map(func,LIM1,conc))) 
    return LIM1_sum
# %%%% 1.2.24  function used to coculate the list_mass of LIM1
def coculate_LIM1_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of LIM1
    input:DataPath,DataName
    output:LIM1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    LIM1 = mass_T[11]
    return LIM1
# %%%% 1.2.25  function used to coculate the sum_mass of LIM2
def coculateC_LIM2_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of LIM2 (contain conc value)
    input:DataPath,DataName
    output:LIM2_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    LIM2 = mass_T[12]
    func = lambda x,y:x*y
    LIM2_sum = sum(list(map(func,LIM2,conc))) 
    return LIM2_sum
# %%%% 1.2.26  function used to coculate the list_mass of LIM2
def coculate_LIM2_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of LIM2
    input:DataPath,DataName
    output:LIM2"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    LIM2 = mass_T[12]
    return LIM2
# %%%% 1.2.27  function used to coculate the sum_mass of CO3
def coculateC_CO3_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of CO3 (contain conc value)
    input:DataPath,DataName
    output:CO3_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    CO3 = mass_T[13]
    func = lambda x,y:x*y
    CO3_sum = sum(list(map(func,CO3,conc))) 
    return CO3_sum
# %%%% 1.2.28  function used to coculate the list_mass of CO3
def coculate_CO3_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of CO3
    input:DataPath,DataName
    output:CO3"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    CO3 = mass_T[13]
    return CO3
# %%%% 1.2.29  function used to coculate the sum_mass of Na
def coculateC_Na_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of Na (contain conc value)
    input:DataPath,DataName
    output:Na_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Na = mass_T[14]
    func = lambda x,y:x*y
    Na_sum = sum(list(map(func,Na,conc))) 
    return Na_sum
# %%%% 1.2.30  function used to coculate the list_mass of Na
def coculate_Na_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of Na
    input:DataPath,DataName
    output:Na"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Na = mass_T[14]
    return Na
# %%%% 1.2.31  function used to coculate the sum_mass of Ca
def coculateC_Ca_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of Ca (contain conc value)
    input:DataPath,DataName
    output:Ca_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Ca = mass_T[15]
    func = lambda x,y:x*y
    Ca_sum = sum(list(map(func,Ca,conc))) 
    return Ca_sum
# %%%% 1.2.32  function used to coculate the list_mass of Ca
def coculate_Ca_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of Ca
    input:DataPath,DataName
    output:Ca"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Ca = mass_T[15]
    return Ca
# %%%% 1.2.33  function used to coculate the sum_mass of OIN
def coculateC_OIN_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of OIN (contain conc value)
    input:DataPath,DataName
    output:OIN_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OIN = mass_T[16]
    func = lambda x,y:x*y
    OIN_sum = sum(list(map(func,OIN,conc))) 
    return OIN_sum
# %%%% 1.2.34  function used to coculate the list_mass of OIN
def coculate_OIN_masslist_func(DataPath,DataName):
    """function:coculate the sum_mass of OIN
    input:DataPath,DataName
    output:OIN"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OIN = mass_T[16]
    return OIN
# %%%% 1.2.35  function used to coculate the sum_mass of OC
def coculateC_OC_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of OC (contain conc value)
    input:DataPath,DataName
    output:OC_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OC = mass_T[17]
    func = lambda x,y:x*y
    OC_sum = sum(list(map(func,OC,conc))) 
    return OC_sum
# %%%% 1.2.36  function used to coculate the list_mass of OC
def coculate_OC_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of OC
    input:DataPath,DataName
    output:OC"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OC = mass_T[17]
    return OC
# %%%% 1.2.37  function used to coculate the sum_mass of BC
def coculateC_BC_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of BC (contain conc value)
    input:DataPath,DataName
    output:BC_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    BC = mass_T[18]
    func = lambda x,y:x*y
    BC_sum = sum(list(map(func,BC,conc))) 
    return BC_sum
# %%%% 1.2.38  function used to coculate the list_mass of BC
def coculate_BC_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of BC
    input:DataPath,DataName
    output:BC"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    BC = mass_T[18]
    return BC
# %%%% 1.2.39  function used to coculate the sum_mass of H2O
def coculateC_H20_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of H2O (contain conc value)
    input:DataPath,DataName
    output:H2O_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    H20 = mass_T[19]
    func = lambda x,y:x*y
    H20_sum = sum(list(map(func,H20,conc))) 
    return H20_sum
# %%%% 1.2.40  function used to coculate the sum_mass of H2O
def coculate_H20_masslist_func(DataPath,DataName):
    """function:coculate the list_mass of H2O
    input:DataPath,DataName
    output:H2O"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    H20 = mass_T[19]
    return H20
# %%%% 1.2.41  function used to coculate the list_mass of particle
def coculate_particle_listmass_func(DataPath,DataName):
    """function:coculate the list_mass of Particle
    input:DataPath,DataName
    output:particle_mass_list"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass']) 
    mass_T = mass
    mass = mass.T
    mass_T = mass_T.tolist()
    mass = mass.tolist()
    particle_mass_list = []
    for i in range(len(mass)):
        particle_mass_list.append(sum(mass[i]))
    return particle_mass_list

# %%%% 1.2.42  function used to coculate the sum_mass of particle
def coculate_particle_summass_func(DataPath,DataName):
    """function:coculate the list_mass of Particle
    input:DataPath,DataName
    output:particle_mass_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass_T = mass_T.tolist()
    mass = mass.tolist()
    particle_mass_sum = 0
    for i in range(len(mass)):
        particle_mass_sum += conc[i] * sum(mass[i])
    return particle_mass_sum

# %%%% 1.2.43  function used to coculate the sum_mass of particle(contain BC)
def coculate_particle_containBC_summass_func(DataPath,DataName):
    """function:coculate the list_mass of Particle
    input:DataPath,DataName
    output:particle_mass_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass_T = mass
    mass = mass.T
    mass_T = mass_T.tolist()
    mass = mass.tolist()
    particle_containBC_mass_sum = 0
    for i in range(len(mass)):
        if mass[i][18] > 0:
            particle_containBC_mass_sum += conc[i] * sum(mass[i])
    return particle_containBC_mass_sum



# %%% 1.3  function used to process the concentration 

# %%%% 1.3.1  function used to process the number concentration of particle and particle containing BC
def coculateC_conc_func(DataPath,DataName):
    """function:process the number concentration of particle and particle containing BC (contain conc value)
    input:DataPath,DataName
    output:conc,conc_sum,BC_conc_sum"""
    nf = nc.Dataset(DataPath + DataName,'r')
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass = np.array(nf.variables['aero_particle_mass'])
    BC_mass = mass[18] # BC质量
    BC_mass = BC_mass.tolist()
   # BC_conc_sum= 0
    #for i in range(len(mass)):
    #    if mass[i][18]!=0:
    #        BC_conc_sum += conc[i]
    BC_conc = [conc[i] for i in range(len(conc)) if BC_mass[i] > 0]
    BC_conc_sum = sum(BC_conc)
    conc_sum = sum(conc)
    return conc,conc_sum,BC_conc_sum

# %%%% 1.3.2  function used to process the concentration of 2 kinds of gass(S&Cl)
def coculateC_S_Cl_mass_func(DataPath,DataName):
    """function:process the concentration of S & Cl
    input:DataPath,DataName
    output:gas_S_mass, gas_Cl_mass"""
    nf = nc.Dataset(DataPath + DataName,'r')
    gas_mixing_ratio = np.array(nf.variables['gas_mixing_ratio'])

    func = lambda x:x*1e-3
    gas_conc = list(map(func,gas_mixing_ratio))
    gas_S_mass = gas_conc[0]*32/98+gas_conc[17]
    gas_Cl_mass = gas_conc[2]*35.5/36.5
    return gas_S_mass, gas_Cl_mass
# %%% 1.4 function used to process background gass(ppb)
def coculate_backgass_func(DataPath,DataName):
    """function:process the concentration of S & Cl
    input:DataPath,DataName
    output:gas_S_mass, gas_Cl_mass"""
    nf = nc.Dataset(DataPath + DataName,'r')
    gas_mixing_ratio = np.array(nf.variables['gas_mixing_ratio'])
    gas_mixing_ratio = gas_mixing_ratio.tolist()
    return gas_mixing_ratio    

# %%% 1.5 function used to coculate average of mass_coating to mass_BCcore
def coculate_mass_coating_BCcore_func(DataPath,DataName,mass_core = 1e-18):
    """function:show the growth rate of coating
    input:DataPath, DataName
    output:value of mass_coating/mass_BCcore"""
    nf = nc.Dataset(DataPath + DataName,'r')
    aero_conc = np.array(nf.variables['aero_num_conc'])
    conc = aero_conc.tolist()
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    BC_mass = mass_T[18] # BC质量
    list_mass_coating_BCcore = []
    BC_conc = []
    for i in range(len(BC_mass)):
        if mass_core*1.1 > BC_mass[i] > mass_core*0.9:
            k = (sum(mass[i])-BC_mass[i])/BC_mass[i]
            list_mass_coating_BCcore.append(k)
            BC_conc.append(conc[i])
    print(sum(BC_conc))
    mass_coating_BCcore =sum([BC_conc[i]*list_mass_coating_BCcore[i] for i in range(len(BC_conc))])/sum(BC_conc)
    return mass_coating_BCcore 

# %% 2  function used to process the NCfile without 'aero_num_conc'

# %%% 2.1  function used to coculate the value about mix state
# %%%% 2.1.1  This function is used for coculate the X of particle 
# input DataPath,DataName
# output X
def coculateN_X_func(DataPath,DataName):
    """coculate the X of particle (without conc value) 
    input:DataPath,DataName
    output:X"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    
    Ha = 0 
    U = 0
    
    for i in range(len(mass)):
        U+=sum(mass[i])
        
    
    for i in range(len(mass)):
        Pi = sum(mass[i])/U
        Hi = 0
        for a in range(len(mass_T)):
            if mass[i][a]!=0:
                Pia = mass[i][a]/sum(mass[i])
                Hi -= Pia*math.log(Pia)
        Ha += Hi*Pi
    Da = math.exp(Ha)
    
    Hr = 0
    for a in range(len(mass_T)):
        Ua = 0
        for i in range(len(mass)):
            Ua += mass[i][a]
        Pa = Ua/U
        if Pa != 0:
            Hr -= Pa*math.log(Pa)
    Dr = math.exp(Hr)
    X = (Da-1)/(Dr-1)
    return X
# %%%% 2.1.2  This function is used for coculate the Opt_X of particle 
# input DataPath,DataName
# output Opt_X
def coculateN_Opt_X_func(DataPath,DataName):
    """coculate the optical_X of particle (without conc value) 
    input:DataPath,DataName
    output:Opt_X"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass = mass.T
    mass = mass.tolist()
    
    Opt_Ha = 0 
    Opt_U = 0
    Opt_U_BC = 0
    for i in range(len(mass)):
        Opt_U+=sum(mass[i])
        
    for i in range(len(mass)):
        Opt_Pi = sum(mass[i])/Opt_U
        Opt_Pi_BC = mass[i][18]/sum(mass[i])
        Opt_Pi_NBC = 1-Opt_Pi_BC
        if Opt_Pi_BC!=0:
            Opt_Hi = -Opt_Pi_BC*math.log(Opt_Pi_BC)-Opt_Pi_NBC*math.log(Opt_Pi_NBC)
            Opt_Ha += Opt_Hi*Opt_Pi
        Opt_U_BC+=mass[i][18]
    Opt_P_BC = Opt_U_BC/Opt_U
    Opt_Da = math.exp(Opt_Ha)
    Opt_P_NBC = 1-Opt_P_BC
    Opt_Hr = -Opt_P_BC*math.log(Opt_P_BC)-Opt_P_NBC*math.log(Opt_P_NBC)
    Opt_Dr = math.exp(Opt_Hr)
    Opt_X = (Opt_Da-1)/(Opt_Dr-1)
    return Opt_X



# %%%% 2.1.3  function used to coculate the values of diameter of particle 
# input DataPath,DataName
# output Dc_Dp_mean,Dc,Dp,Dc_Dp
def coculateN_DcDp_conc_func(DataPath,DataName):
    """coculate the values of diameter of particle (without conc value) 
    input:DataPath,DataName
    output:Dc_Dp_mean,Dc,Dp,Dc_Dp"""
    nf = nc.Dataset(DataPath + DataName,'r')
    aero_den = np.array(nf.variables['aero_density'])
    mass = np.array(nf.variables['aero_particle_mass'])
    
    den = aero_den.tolist()
    BC_mass = mass[18] # BC质量
    BC_mass = BC_mass.tolist()
    
    
    P_mass = [] #单个particle总质量
    aero_den = [] # particle密度
    Dc = [] # BC核直径
    Dp = [] # particle的直径
    Dc_Dp = []
    mass = mass.T
    for i in range(len(BC_mass)):
        P_mass.append(sum(mass[i]))
        density = 0
        for j in range(len(den)):
            density = density + den[j]*mass[i,j]/sum(mass[i])
        aero_den.append(density)
        Dc.append(math.pow(BC_mass[i]/den[18]*3/4/3.14,1/3)*1e9*2)
        Dp.append(math.pow(P_mass[i]/aero_den[i]*3/4/3.14, 1/3)*1e9*2)
        Dc_Dp.append((math.pow(BC_mass[i]/den[18]*3/4/3.14,1/3)*1e9*2)/(math.pow(P_mass[i]/aero_den[i]*3/4/3.14, 1/3)*1e9*2))
    Dc_Dp_sum = 0
    for i in range(len(Dc_Dp)):
        Dc_Dp_sum += Dc_Dp[i]
    Dc_Dp_mean = Dc_Dp_sum/len(Dc_Dp)
    return Dc_Dp_mean,Dc,Dp,Dc_Dp
# %%% 2.2 function used to coculate the values of mass 
# %%%% 2.2.1  function used to coculate the sum_mass of SO4
def coculateN_SO4_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of SO4 (without conc value)
    input:DataPath,DataName
    output:SO4"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    SO4 = mass_T[0]
    SO4_sum = sum(SO4) 
    return SO4_sum
# %%%% 2.2.2  function used to coculate the sum_mass of NO3
def coculateN_NO3_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of NO3 (without conc value)
    input:DataPath,DataName
    output:NO3"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    NO3 = mass_T[1]
    NO3_sum = sum(NO3) 
    return NO3_sum
# %%%% 2.2.3  function used to coculate the sum_mass of Cl
def coculateN_Cl_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of Cl (without conc value)
    input:DataPath,DataName
    output:Cl"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Cl = mass_T[2]
    Cl_sum = sum(Cl) 
    return Cl_sum
# %%%% 2.2.4  function used to coculate the sum_mass of NH4
def coculateN_NH4_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of Cl (without conc value)
    input:DataPath,DataName
    output:Cl"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    NH4 = mass_T[3]
    NH4_sum = sum(NH4) 
    return NH4_sum
# %%%% 2.2.5  function used to coculate the sum_mass of MSA
def coculateN_MSA_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of MSA (without conc value)
    input:DataPath,DataName
    output:MSA"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    MSA = mass_T[4]
    MSA_sum = sum(MSA) 
    return MSA_sum
# %%%% 2.2.6  function used to coculate the sum_mass of ARO1
def coculateN_ARO1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of ARO1 (without conc value)
    input:DataPath,DataName
    output:ARO1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ARO1 = mass_T[5]
    ARO1_sum = sum(ARO1) 
    return ARO1_sum
# %%%% 2.2.7  function used to coculate the sum_mass of ARO2
def coculateN_ARO2_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of ARO2 (without conc value)
    input:DataPath,DataName
    output:ARO2"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ARO2 = mass_T[6]
    ARO2_sum = sum(ARO2) 
    return ARO2_sum
# %%%% 2.2.8  function used to coculate the sum_mass of ALK1
def coculateN_ALK1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of ALK1 (without conc value)
    input:DataPath,DataName
    output:ALK1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    ALK1 = mass_T[7]
    ALK1_sum = sum(ALK1) 
    return ALK1_sum
# %%%% 2.2.9  function used to coculate the sum_mass of OLE1
def coculateN_OLE1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of OLE1 (without conc value)
    input:DataPath,DataName
    output:OLE1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OLE1 = mass_T[8]
    OLE1_sum = sum(OLE1) 
    return OLE1_sum
# %%%% 2.2.10  function used to coculate the sum_mass of API1
def coculateN_API1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of API1 (without conc value)
    input:DataPath,DataName
    output:API1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    API1 = mass_T[9]
    API1_sum = sum(API1) 
    return API1_sum
# %%%% 2.2.11  function used to coculate the sum_mass of API2
def coculateN_API2_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of API2 (without conc value)
    input:DataPath,DataName
    output:API2"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    API2 = mass_T[10]
    API2_sum = sum(API2) 
    return API2_sum
# %%%% 2.2.12  function used to coculate the sum_mass of LIM1
def coculateN_LIM1_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of LIM1 (without conc value)
    input:DataPath,DataName
    output:LIM1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    LIM1 = mass_T[11]
    LIM1_sum = sum(LIM1) 
    return LIM1_sum
# %%%% 2.2.13  function used to coculate the sum_mass of LIM2
def coculateN_LIM2_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of LIM2 (without conc value)
    input:DataPath,DataName
    output:LIM2"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    LIM2 = mass_T[12]
    LIM2_sum = sum(LIM2) 
    return LIM2_sum
# %%%% 2.2.14  function used to coculate the sum_mass of CO3
def coculateN_CO3_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of CO3 (without conc value)
    input:DataPath,DataName
    output:CO3"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    CO3 = mass_T[13]
    CO3_sum = sum(CO3) 
    return CO3_sum
# %%%% 2.2.15  function used to coculate the sum_mass of Na
def coculateN_Na_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of Na (without conc value)
    input:DataPath,DataName
    output:Na"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Na = mass_T[14]
    Na_sum = sum(Na) 
    return Na_sum
# %%%% 2.2.16  function used to coculate the sum_mass of Ca
def coculateN_Ca_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of Ca (without conc value)
    input:DataPath,DataName
    output:Ca"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    Ca = mass_T[15]
    Ca_sum = sum(Ca) 
    return Ca_sum
# %%%% 2.2.17  function used to coculate the sum_mass of OIN
def coculateN_OIN_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of OIN (without conc value)
    input:DataPath,DataName
    output:OIN"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OIN = mass_T[16]
    OIN_sum = sum(OIN) 
    return OIN_sum
# %%%% 2.2.18  function used to coculate the sum_mass of OC
def coculateN_OC_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of OC (without conc value)
    input:DataPath,DataName
    output:OC"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    OC = mass_T[17]
    OC_sum = sum(OC) 
    return OC_sum
# %%%% 2.2.19  function used to coculate the sum_mass of BC
def coculateN_BC_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of BC (without conc value)
    input:DataPath,DataName
    output:BC"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    BC = mass_T[18]
    BC_sum = sum(BC) 
    return BC_sum
# %%%% 2.2.20  function used to coculate the sum_mass of H2O
def coculateN_H2O_masssum_func(DataPath,DataName):
    """function:coculate the sum_mass of H2O (without conc value)
    input:DataPath,DataName
    output:H2O"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass'])
    mass_T = mass
    mass = mass.T
    mass = mass.tolist()
    mass_T = mass_T.tolist()
    H2O = mass_T[19]
    H2O_sum = sum(H2O) 
    return H2O_sum

# %%%% 2.3  function used to coculate the number concentration 
# %%%% 2.3.1  function used to process the number concentration of particle and particle containing BC

def coculateN_conc_func(DataPath,DataName):
    """function:process the number concentration of particle and particle containing BC (without conc value)
    input:DataPath,DataName
    output:API1"""
    nf = nc.Dataset(DataPath + DataName,'r')
    mass = np.array(nf.variables['aero_particle_mass']) 
    mass = mass.T
    mass = mass.tolist()
    BC_conc = 0
    for i in range(len(mass)):
        if mass[i][18]!=0:
            BC_conc += 1
    conc = len(mass)
    return conc,BC_conc


        
