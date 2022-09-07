# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 12:26:20 2022

The purpose of this script is to process/plot the data collected from the 
'ablation stake garden' on the Kaskawulsh in August 2022 (GF, KR, TH)

@author: katierobinson
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

stakedata ='F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/AblationStakeData_2022.csv'

df = pd.read_csv(stakedata)
    #df['Tributary'].astype('Int64')
arr = np.array(df)

sitename = arr[:,0]
initial_height = arr[:,1]
initial_uncertainty = arr[:,2]
final_height_min = arr[:,3]
final_height_max = arr[:,4]
final_uncertainty = arr[:,5]

# average final height min and max to get abs. final height
#subtract initial height from final height to get abs. change 

final_height = np.zeros(final_height_max.shape)
for i in range(0,len(final_height_max)):
    final = np.nanmean([final_height_max[i],final_height_min[i]])
    final_height[i] = final

height_change = (final_height - initial_height)/100


#load debris thickness data

debdata ='F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/Stakedebrisdata_2022.csv'

df2 = pd.read_csv(debdata)
    #df['Tributary'].astype('Int64')
arr2 = np.array(df2)

# find avg debris thickness at each site by taking the mean of the 12 deb measurements

CI00_i = (arr2[:,1])
DI00_i = (arr2[:,3])
DB01_i = (arr2[:,5])
DB02_i = (arr2[:,7])
DB03_i = (arr2[:,9])
DB04_i = (arr2[:,11])
MR07_i = (arr2[:,13])
KV04_i = (arr2[:,15])

debris_initial = np.zeros(len(sitename))
initial_debris_lists = [CI00_i,DI00_i,DB01_i,DB02_i,DB03_i,DB04_i,MR07_i,KV04_i]
for i in range(0,len(initial_debris_lists)):
    debris_initial[i] = np.mean(initial_debris_lists[i])

CI00_f = (arr2[:,18])
DI00_f = (arr2[:,20])
DB01_f = (arr2[:,22])
DB02_f = (arr2[:,24])
DB03_f = (arr2[:,26])
DB04_f = (arr2[:,28])
MR07_f = (arr2[:,30])
KV04_f = (arr2[:,32])


debris_final = np.zeros(len(sitename))
final_debris_lists = [CI00_f,DI00_f,DB01_f,DB02_f,DB03_f,DB04_f,MR07_f,KV04_f]
for i in range(0,len(final_debris_lists)):
    debris_final[i] = np.mean(final_debris_lists[i])

debris_average = np.zeros(len(debris_final))
for i in range(0,len(debris_final)):
    avg = np.nanmean([debris_final[i],debris_initial[i]])
    debris_average[i] = avg


#calculate uncertainties on STAKE HEIGHT: 
#upper bound:
#smallest possible initial and biggest possible final:
min_init = initial_height - initial_uncertainty
max_final = final_height_max + final_uncertainty
upper_bound = (max_final - min_init)/100

#lower bound:
#tallest possible intial and smallest possible final
max_init = initial_height + initial_uncertainty
min_final = final_height_min - final_uncertainty
lower_bound = (min_final - max_init)/100

stake_uncertainty_range = (upper_bound - lower_bound)/2

#OR calculate it based on rules of addition/subtraction
stake_uncertainty = np.zeros(initial_uncertainty.shape)
for i in range(0,len(stake_uncertainty)):
    uncertainty = (np.sqrt((initial_uncertainty[i]**2)+(final_uncertainty[i]**2)))/100
    stake_uncertainty[i] = uncertainty

#calculate uncertainties on DEBRIS THICKNESS:
# one uncertainty range per site: 
#uncertainties on initial debris thicknesses:
uCI00_i = (arr2[:,2])
uDI00_i = (arr2[:,4])
uDB01_i = (arr2[:,6])
uDB02_i = (arr2[:,8])
uDB03_i = (arr2[:,10])
uDB04_i = (arr2[:,12])
uMR07_i = (arr2[:,14])
uKV04_i = (arr2[:,16])

initial_deb_uncertainties_list = [uCI00_i,uDI00_i,uDB01_i,uDB02_i,uDB03_i,uDB04_i,uMR07_i,uKV04_i]

uCI00_f = (arr2[:,19])
uDI00_f = (arr2[:,21])
uDB01_f = (arr2[:,23])
uDB02_f = (arr2[:,25])
uDB03_f = (arr2[:,27])
uDB04_f = (arr2[:,29])
uMR07_f = (arr2[:,31])
uKV04_f = (arr2[:,33])

final_deb_uncertainties_list = [uCI00_f,uDI00_f,uDB01_f,uDB02_f,uDB03_f,uDB04_f,uMR07_f,uKV04_f]

for i in range(0,len(initial_deb_uncertainties_list)):
    min_deb_init = np.mean(initial_debris_lists[i] - initial_deb_uncertainties_list[i])
    max_deb_init = np.mean(initial_debris_lists[i] + initial_deb_uncertainties_list[i])
    deb_init_range = (max_deb_init - min_deb_init)/2
    
    max_deb_init = np.mean(final_debris_lists[i] - final_deb_uncertainties_list[i])
    

#calculate the uncertainty on the debris for EACH SITE:
def calculate_debristhickness_uncertainty(debris_site_thicknesses,meanthicknesses):
    '''
    calculate the uncertainty on the average debris thickness (from the 
    12 debris thickness measurements taken in field) using 
    the formula from: 
    https://www.sml-inc.com/uncertainy.htmhttp://science.clemson.edu/physics/labs/tutorials/errorp/index.html
    
    debris_site_thicknesses is a list of 8 lists containing the 12 deb measurements from each site
    meanthicknesses is a list of the mean deb thickness at each site (8 values total)
    '''
    debris_uncertainty = np.zeros(len(debris_site_thicknesses))
    for i in range(0,len(debris_site_thicknesses)):
        x1 = debris_site_thicknesses[i][0]
        x2 = debris_site_thicknesses[i][1]
        x3 = debris_site_thicknesses[i][2]
        x4 = debris_site_thicknesses[i][3]
        x5 = debris_site_thicknesses[i][4]
        x6 = debris_site_thicknesses[i][5]
        x7 = debris_site_thicknesses[i][6]
        x8 = debris_site_thicknesses[i][7]
        x9 = debris_site_thicknesses[i][8]
        x10 = debris_site_thicknesses[i][9] 
        x11 = debris_site_thicknesses[i][10]
        x12 = debris_site_thicknesses[i][11]
        xbar = meanthicknesses[i]
        N = len(debris_site_thicknesses[i])
        
        debris_uncertainty[i] = np.sqrt(((x1-xbar)**2 + (x2-xbar)**2 + (x3-xbar)**2 + (x4-xbar)**2 + (x5-xbar)**2 + (x6-xbar)**2 + (x7-xbar)**2 + (x8-xbar)**2 + (x9-xbar)**2 + (x10-xbar)**2 + (x11-xbar)**2 + (x12-xbar)**2)/(N-1))
        
    return debris_uncertainty

initial_debris_uncertainty = calculate_debristhickness_uncertainty(initial_debris_lists,debris_initial)
final_debris_uncertainty = calculate_debristhickness_uncertainty(final_debris_lists,debris_final)

interp_initdebris = interp1d(debris_initial,height_change,kind='linear')
x_init = np.arange(0,12.3,0.1)
y_init = interp_initdebris(x_init)

interp_findebris = interp1d(debris_final,height_change,kind='linear')
x_fin = np.arange(0,5,0.1)
y_fin = interp_findebris(x_fin)

interp_avgdebris = interp1d(debris_average,height_change,kind='linear')
x_avg = np.arange(0,8.7,0.1)
y_avg = interp_avgdebris(x_avg)

#OSTREM CURVE FROM LOOMIS (1970)
loomis ='F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/Loomis1970_data.csv'
loomis_avg ='F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/Loomis1970_data_averaged.csv'

df = pd.read_csv(loomis)
arr3 = np.array(df)
loomis_abl = arr3[:,0]
loomis_till = arr3[:,1]

df = pd.read_csv(loomis_avg)
arr4 = np.array(df)
loomis_abl_avg = arr4[:,0]
loomis_till_avg = arr4[:,1]

interploomis = interp1d(loomis_till_avg,loomis_abl_avg,kind='quadratic')

x_loomis = np.arange(0,37.5,0.1)
y_loomis = interploomis(x_loomis)
######################################################################################################
########################################PLOTS#########################################################
######################################################################################################

plt.figure(figsize=(8,8))
plt.title('Debris Thickness vs. Ablation',fontsize=12)
#plt.plot(debris_initial,height_change,'royalblue')
#plt.plot(debris_final,height_change,'red')
plt.plot(x_init,y_init,c='royalblue')
plt.plot(x_fin,y_fin,c='red')
#plt.plot(debris_average,height_change,'orange')
#plt.scatter(debris_initial,height_change,c='royalblue')
#plt.scatter(debris_final,height_change,c='red')
#plt.scatter(debris_average,height_change,c='orange')
plt.errorbar(debris_initial,height_change,stake_uncertainty,c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_final,height_change,stake_uncertainty,c='red',fmt="o",capsize=5)
plt.errorbar(debris_initial,height_change,yerr=None,xerr=initial_debris_uncertainty,c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_final,height_change,yerr=None,xerr=final_debris_uncertainty,c='red',fmt="o",capsize=5)
#plt.errorbar(debris_average,height_change,uncertainty_range,c='orange',fmt="o",capsize=5)
plt.legend(['Initial debris thickness (July)','Final debris thickness (Aug)'],fontsize=12)
plt.xlabel('Debris Thickness (cm)',fontsize=12)
plt.ylabel('Change in Stake Height (m)',fontsize=12)
plt.xlim(0,15)
#plt.savefig('stakegarden_data.png',bbox_inches = 'tight')


plt.figure(figsize=(8,6))
plt.plot(x_init,y_init)
plt.scatter(debris_initial,height_change)
plt.plot(x_avg,y_avg)
plt.scatter(debris_average,height_change)
plt.plot(x_fin,y_fin)
plt.scatter(debris_final,height_change)
plt.ylim(1,2.5)

#ablation stakes were out for 43 DAYS
#convert stake heights to ablation in cm/day
height_change_cmday = height_change/43*100

interp_cmday_kmr_init = interp1d(debris_initial,height_change_cmday,kind='quadratic')
x_kmr_i = np.arange(0,12.3,0.1)
y_kmr_i = interp_cmday_kmr_init(x_kmr_i)

interp_cmday_kmr_fin = interp1d(debris_final,height_change_cmday,kind='quadratic')
x_kmr_f = np.arange(0,4.95,0.01)
y_kmr_f = interp_cmday_kmr_fin(x_kmr_f)

loomis4 = np.poly1d(np.polyfit(loomis_till,loomis_abl, 4))
polylineloomis = np.linspace(0, 37.5, 100)

plt.figure(figsize=(7,7))
plt.scatter(loomis_till,loomis_abl,c='orange')
plt.plot(polylineloomis, loomis4(polylineloomis), color='orange')
plt.scatter(debris_initial,height_change_cmday,c='royalblue')
plt.plot(x_kmr_i,y_kmr_i,'royalblue')
plt.scatter(debris_final,height_change_cmday,c='red')
plt.plot(x_kmr_f,y_kmr_f,'red')
plt.xlim(0,40)
plt.xlabel('Debris Thickness (cm)',fontsize=12)
plt.ylabel('Change in Stake Height (cm/day)',fontsize=12)
plt.legend(['Loomis (1970)','2022 Field Data (July thicknesses)','2022 Field Data (Aug thicknesses)'],fontsize=12)
#plt.savefig('2022vs1970_curves.png',bbox_inches = 'tight')

#FIND PEAK MELT AND TRANSITION THICKNESS:
#PEAK MELT:
peakmelt_i = np.where(y_init == np.max(y_init))
peakmelt_thickness_i = x_init[peakmelt_i]
print('Jul obs: debris thickness resulting in peak melt is: ' + str(peakmelt_thickness_i))
peakmelt_f = np.where(y_fin == np.max(y_fin))
peakmelt_thickness_f = x_fin[peakmelt_f]
print('Aug obs: debris thickness resulting in peak melt is: ' + str(peakmelt_thickness_f))

#transition thickness:
diff_from_CI00_i = [99]
for i in y_init[1:]:
    diff_from_CI00_i.append(i-y_init[0])
diff_from_cleanice_i = np.abs(diff_from_CI00_i)
    
#minafterpeak = np.absmin(diff_from_cleanice[peakmelt_i[0][0]:])
CI00_eq_melt_i = np.where(diff_from_cleanice_i == np.min(diff_from_cleanice_i[peakmelt_i[0][0]:]))
transition_thickness_i = x_init[CI00_eq_melt_i]

diff_from_CI00_f = [99]
for i in y_fin[1:]:
    diff_from_CI00_f.append(i-y_fin[0])
diff_from_cleanice_f = np.abs(diff_from_CI00_f)
    
#minafterpeak = np.absmin(diff_from_cleanice[peakmelt_i[0][0]:])
CI00_eq_melt_f = np.where(diff_from_cleanice_f == np.min(diff_from_cleanice_f[peakmelt_f[0][0]:]))
transition_thickness_f = x_fin[CI00_eq_melt_f]    
print('Jul obs: transition thickness = ' + str(transition_thickness_i))
print('Aug obs: transition thickness = ' + str(transition_thickness_f))
