#Dp_Dc_distribution
#show Dc distribution and Dp distribution and k
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import os
import NC_X
from scipy import optimize
import csv
font_manager.fontManager.addfont('/data/home/zzy/.fonts/ARIAL.TTF')
plt.rcParams['font.sans-serif']= 'ARIAL'

#f_1
def f_1(x,A,B):
    return A * x + B
#log func
def logF(va):
    if va == 0:
        return 0
    else:
        return math.log(va)


DataPath = '../Data/'
DataName = 'Dp_Dc_Conc.csv'
FigurePath = '../Fig/'

Dp_Dc_Conc = []
#Read
with open(DataPath+DataName, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        float_row = [float(x) for x in row]
        Dp_Dc_Conc.append(float_row)
print('validation of list_len')
print(len(Dp_Dc_Conc))
Dp = Dp_Dc_Conc[:240]
Dc = Dp_Dc_Conc[240:480]
conc = Dp_Dc_Conc[480:]
print(len(Dp))
print(len(Dc))
print(len(conc))
print()
File_number = len(Dp)


CT = []
CONC_CT=[]
conc_nDp_sum = 0
conc_nDp_sum40_120 = 0

#filter 96h after
conc = conc[48:]
Dc = Dc[48:]
Dp = Dp[48:]
for i in range(len(conc)):
    for j in range(len(conc[i])):
        if Dc[i][j]>0:
            CT.append(Dp[i][j] - Dc[i][j])
            CONC_CT.append(conc[i][j])


sum_CT=0
sum_CT_600 = 0
conc_CT_600 = 0
for i in range(len(CONC_CT)):
    sum_CT+=CONC_CT[i]*CT[i]
    if CT[i]<=600:
        sum_CT_600+=CONC_CT[i]*CT[i]
        conc_CT_600+=CONC_CT[i]
mean_CT = sum_CT/sum(CONC_CT)
mean_CT_600 = sum_CT_600/conc_CT_600
print('mean_CT:',end='')
print(mean_CT)
print('mean_CT_600:',end='')
print(mean_CT_600)
print()
