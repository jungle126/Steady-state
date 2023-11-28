import numpy as np
import csv
import os
import NC_X

#set datapath
DataPath = '../out/'

#set file save path and filename
SavePath = '../Data/'
CsvName = 'Dp_Dc_Conc.csv'

#set filename in three paths,share the same name
lists = []
for file in os.listdir(DataPath):
   if os.path.splitext(file)[1].lower() in '.nc':
       lists.append(file)
lists = sorted(lists,key=lambda keys:[ord(i) for i in keys],reverse=False)
listtime = [i+1 for i in range(240)] #set time 1-240h
DataName = lists[1:]  #name list 240

#set blank list [ScA_Dp,ScA_Dc,ScA_Conc,ScB_Dp,ScB_Dc,ScB_Conc,ScC_Dp,ScC_Dc,ScC_Conc]
Sc_Dp = [[] for i in range(len(DataName))]  # 1*240 scenario Dp
Sc_Dc = [[] for i in range(len(DataName))]  # 1*240 scenario Dc
Sc_Conc = [[] for i in range(len(DataName))]# 1*240 scenario Conc

#get the Dp Dc Conc
for i in range(len(Sc_Dp)):
    Dc_Dp_mean,Sc_Dc[i],Sc_Dp[i],Dc_Dp = NC_X.coculateC_DcDp_func(DataPath, DataName[i])
    Sc_Conc[i],conc_sum,BC_conc_sum = NC_X.coculateC_conc_func(DataPath, DataName[i])

#put all data together into fomation like
#SC_Dp1
#Sc_Dc1
#Sc_Conc1
#...
datalist = []
datalist.extend(Sc_Dp)
datalist.extend(Sc_Dc)
datalist.extend(Sc_Conc)
print(len(datalist))

#save csv file
with open(SavePath+CsvName, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    for lst in datalist:
        csvwriter.writerow(lst)

