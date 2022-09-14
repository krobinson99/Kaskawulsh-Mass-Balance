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
import matplotlib.dates as mdates
#import os
#from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.lines as mlines
#from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import splrep
from scipy.interpolate import splev
#from sklearn.metrics import mean_squared_error
import datetime as dt

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

height_change = np.array((final_height - initial_height)/100,dtype=np.float)


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
#plt.xlim(0,15)
plt.margins(x=0.01)
#plt.savefig('stakegarden_data.png',bbox_inches = 'tight')


plt.figure(figsize=(8,6))
plt.plot(x_init,y_init)
plt.scatter(debris_initial,height_change)
plt.plot(x_avg,y_avg)
plt.scatter(debris_average,height_change)
plt.plot(x_fin,y_fin)
plt.scatter(debris_final,height_change)
plt.ylim(1,2.5)

#average debris uncertainty
average_debris_uncertainty = np.zeros(initial_debris_uncertainty.shape)
for i in range(0,len(average_debris_uncertainty)):
    average_debris_uncertainty[i] = np.mean([initial_debris_uncertainty[i],final_debris_uncertainty[i]])

blue_line = mlines.Line2D([], [], color='royalblue', marker='_', linestyle='None',
                          markersize=20, linewidth = 100, label='Initial debris thickness (July)')

orange_line = mlines.Line2D([], [], color='orange', marker='_', linestyle='None',
                          markersize=20, label='Average debris thickness')

red_line = mlines.Line2D([], [], color='red', marker='_', linestyle='None',
                          markersize=20, label='Final debris thickness (Aug)')

white_dot = mlines.Line2D([], [], color='white',mec='k', marker='.', linestyle='None',
                          markersize=12, label='Kovacs Drill')


cluster = np.array([1,1,1,1,1,1,1,2]) 
plt.figure(figsize=(10,7))
plt.title('Debris Thickness vs. Ablation',fontsize=14)
plt.plot(debris_initial[:-1],height_change[:-1],c='royalblue')
plt.plot(debris_average[:-1],height_change[:-1],c='orange')
plt.plot(x_fin[:30],y_fin[:30],c='red')

plt.errorbar(debris_initial[cluster==1],height_change[cluster==1],stake_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1],stake_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1],stake_uncertainty[cluster==1],c='red',fmt="o",capsize=5)
plt.errorbar(debris_initial[cluster==1],height_change[cluster==1],yerr=None,xerr=initial_debris_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1],yerr=None,xerr=average_debris_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1],yerr=None,xerr=final_debris_uncertainty[cluster==1],c='red',fmt="o",capsize=5)

plt.errorbar(debris_initial[cluster==2],height_change[cluster==2],stake_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2],stake_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2],stake_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
plt.errorbar(debris_initial[cluster==2],height_change[cluster==2],yerr=None,xerr=initial_debris_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2],yerr=None,xerr=average_debris_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2],yerr=None,xerr=final_debris_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
#plt.legend(['Initial debris thickness (July)','Average debris thickness','Final debris thickness (Aug)'],fontsize=12)
plt.legend(handles=[blue_line,orange_line,red_line,white_dot],fontsize=14)
plt.xlabel('Debris Thickness (cm)',fontsize=14)
plt.ylabel('Change in Stake Height (m)',fontsize=14)
#plt.xlim(0,15)
plt.margins(x=0.01)
#plt.savefig('stakegarden_data_linearfit.png',bbox_inches = 'tight')

#ablation stakes were out for 43 DAYS
#convert stake heights to ablation in cm/day
height_change_cmday = np.array((height_change/43*100),dtype=np.float)

interp_cmday_kmr_init = interp1d(debris_initial,height_change_cmday,kind='linear')
x_kmr_i = np.arange(0,12.3,0.1)
y_kmr_i = interp_cmday_kmr_init(x_kmr_i)

interp_cmday_kmr_fin = interp1d(debris_final,height_change_cmday,kind='linear')
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

plt.figure(figsize=(7,7))
#plt.scatter(loomis_till,loomis_abl,c='orange')
plt.plot(polylineloomis, loomis4(polylineloomis), color='orange')
plt.scatter(debris_initial,height_change_cmday,c='royalblue')
#plt.plot(x_kmr_i,y_kmr_i,'royalblue')
plt.scatter(debris_final,height_change_cmday,c='red')
plt.scatter(loomis_till,loomis_abl,c='orange')
#plt.plot(x_kmr_f,y_kmr_f,'red')
plt.margins(x=0.01)
#plt.xlim(0,40)
plt.xlabel('Debris Thickness (cm)',fontsize=12)
plt.ylabel('Change in Stake Height (cm/day)',fontsize=12)
plt.legend(['Loomis (1970)','2022 Field Data (July thicknesses)','2022 Field Data (Aug thicknesses)'],fontsize=12)
#plt.savefig('2022vs1970_curves_v2.png',bbox_inches = 'tight')

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


###############################################################################
#############FIT A CURVE TO THE DATA POINTS (EXCLUDING KV04) ##################
###############################################################################
fit_initdeb =  np.poly1d(np.polyfit(debris_initial[:-1],height_change[:-1], 3)) #problem is with the height_change dtype = object, fixed by changing dtype to np.float
polyline_initdeb = np.linspace(0,6.85,100)

fit_avgdeb =  np.poly1d(np.polyfit(debris_average[:-1],height_change[:-1], 3)) #use deg=3 for cubic polynomial
polyline_avgdeb = np.linspace(0,4.84,100)

fit_findeb =  np.poly1d(np.polyfit(debris_final[:-1],height_change[:-1], 3))
polyline_findeb = np.linspace(0,2.85,100)

cluster = np.array([1,1,1,1,1,1,1,2]) 
plt.figure(figsize=(10,7))
plt.title('Debris Thickness vs. Ablation (July 19 - Aug 31 2022)',fontsize=14)
plt.plot(polyline_initdeb,fit_initdeb(polyline_initdeb),c='royalblue')
plt.plot(polyline_avgdeb,fit_avgdeb(polyline_avgdeb),c='orange')
plt.plot(polyline_findeb,fit_findeb(polyline_findeb),c='red')

plt.errorbar(debris_initial[cluster==1],height_change[cluster==1],stake_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1],stake_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1],stake_uncertainty[cluster==1],c='red',fmt="o",capsize=5)
plt.errorbar(debris_initial[cluster==1],height_change[cluster==1],yerr=None,xerr=initial_debris_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1],yerr=None,xerr=average_debris_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1],yerr=None,xerr=final_debris_uncertainty[cluster==1],c='red',fmt="o",capsize=5)

plt.errorbar(debris_initial[cluster==2],height_change[cluster==2],stake_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2],stake_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2],stake_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
plt.errorbar(debris_initial[cluster==2],height_change[cluster==2],yerr=None,xerr=initial_debris_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2],yerr=None,xerr=average_debris_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2],yerr=None,xerr=final_debris_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
#plt.legend(['Initial debris thickness (July)','Average debris thickness','Final debris thickness (Aug)'],fontsize=12)
plt.legend(handles=[blue_line,orange_line,red_line,white_dot],fontsize=14)
plt.xlabel('Debris Thickness (cm)',fontsize=14)
plt.ylabel('Change in Stake Height (m)',fontsize=14)
#plt.xlim(0,15)
plt.margins(x=0.01)
#plt.savefig('stakegarden_data_curvefit.png',bbox_inches = 'tight')

###############################################################################
##################COMPARE OBSERVED CURVES TO LOOMIS:###########################
###############################################################################
# normalize curves to cm/day (for 43 days) by multiplying heightchange by (100/43)

plt.figure(figsize=(10,7))
#plt.scatter(loomis_till,loomis_abl,c='orange')
plt.plot(polylineloomis, loomis4(polylineloomis), color='mediumseagreen',linewidth=2.5)
plt.plot(polyline_findeb,fit_findeb(polyline_findeb)*(100/43),c='red',linewidth=2.5)
plt.plot(polyline_avgdeb,fit_avgdeb(polyline_avgdeb)*(100/43),c='orange',linewidth=2.5)
plt.plot(polyline_initdeb,fit_initdeb(polyline_initdeb)*(100/43),c='royalblue',linewidth=2.5)
plt.scatter(loomis_till,loomis_abl,c='mediumseagreen')
plt.scatter(debris_final[cluster==1],height_change_cmday[cluster==1],c='red')
plt.scatter(debris_average[cluster==1],height_change_cmday[cluster==1],c='orange')
plt.scatter(debris_initial[cluster==1],height_change_cmday[cluster==1],c='royalblue')
plt.scatter(debris_final[cluster==2],height_change_cmday[cluster==2],c='white',edgecolors='red')
plt.scatter(debris_average[cluster==2],height_change_cmday[cluster==2],c='white',edgecolors='orange')
plt.scatter(debris_initial[cluster==2],height_change_cmday[cluster==2],c='white',edgecolors='royalblue')
#plt.plot(x_kmr_f,y_kmr_f,'red')
plt.margins(x=0.01)
#plt.xlim(0,40)
plt.xlabel('Debris Thickness (cm)',fontsize=14)
plt.ylabel('Change in Stake Height (cm/day)',fontsize=14)
plt.legend(['Loomis (1970)','2022 Field Data (Aug thicknesses)','Average debris thickness','2022 Field Data (July thicknesses)'],fontsize=12)
#plt.savefig('loomiscomparison.png',bbox_inches = 'tight')

###############################################################################
######################PLOT MELT FACTORS########################################
###############################################################################
#get melt factors by dividing the curved by the clean ice melt (cim)
cim = height_change[0] #clean ice melt

plt.figure(figsize=(10,7))
plt.title('Melt Factors',fontsize=14)
plt.plot(polyline_initdeb,(fit_initdeb(polyline_initdeb))/cim,c='royalblue')
plt.plot(polyline_avgdeb,(fit_avgdeb(polyline_avgdeb))/cim,c='orange')
plt.plot(polyline_findeb,(fit_findeb(polyline_findeb))/cim,c='red')

plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]/cim,stake_uncertainty[cluster==1]/cim,c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]/cim,stake_uncertainty[cluster==1]/cim,c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]/cim,stake_uncertainty[cluster==1]/cim,c='red',fmt="o",capsize=5)
plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]/cim,yerr=None,xerr=initial_debris_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]/cim,yerr=None,xerr=average_debris_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]/cim,yerr=None,xerr=final_debris_uncertainty[cluster==1],c='red',fmt="o",capsize=5)

plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]/cim,stake_uncertainty[cluster==2]/cim,c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]/cim,stake_uncertainty[cluster==2]/cim,c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]/cim,stake_uncertainty[cluster==2]/cim,c='red',fmt="o",capsize=5,mfc='white',mec='red')
plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]/cim,yerr=None,xerr=initial_debris_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]/cim,yerr=None,xerr=average_debris_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]/cim,yerr=None,xerr=final_debris_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
#plt.legend(['Initial debris thickness (July)','Average debris thickness','Final debris thickness (Aug)'],fontsize=12)
plt.legend(handles=[blue_line,orange_line,red_line,white_dot],fontsize=14)
plt.xlabel('Debris Thickness (cm)',fontsize=14)
plt.ylabel('$\dfrac{Total\,Melt}{Clean\,Ice\,Melt}$',fontsize=14)
#plt.xlim(0,15)
plt.margins(x=0.01)
#plt.savefig('stakegarden_meltfactors.png',bbox_inches = 'tight')

####### PLOT STAKE HEIGHT IN M.W.E.############################################
#check: what ice density did EY use in model conversions? ??????
#convert stake height to m.w.e (1 kg/m2 = 1 mm w.e.)
p_ice = 900 #kg/m3
plt.figure(figsize=(10,7))
plt.title('Debris Thickness vs. Ablation (July 19 - Aug 31 2022)',fontsize=14)
plt.plot(polyline_initdeb,(fit_initdeb(polyline_initdeb))*(p_ice/1000),c='royalblue')
plt.plot(polyline_avgdeb,(fit_avgdeb(polyline_avgdeb))*(p_ice/1000),c='orange')
plt.plot(polyline_findeb,(fit_findeb(polyline_findeb))*(p_ice/1000),c='red')

plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='red',fmt="o",capsize=5)
plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=initial_debris_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=average_debris_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=final_debris_uncertainty[cluster==1],c='red',fmt="o",capsize=5)

plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='red',fmt="o",capsize=5,mfc='white',mec='red')
plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=initial_debris_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=average_debris_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=final_debris_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
#plt.legend(['Initial debris thickness (July)','Average debris thickness','Final debris thickness (Aug)'],fontsize=12)
plt.legend(handles=[blue_line,orange_line,red_line,white_dot],fontsize=14)
plt.xlabel('Debris Thickness (cm)',fontsize=14)
plt.ylabel('Melt (m w.e.)',fontsize=14)
#plt.xlim(0,15)
plt.margins(x=0.01)
#plt.savefig('stakegarden_data_curvefit_mwe.png',bbox_inches = 'tight')

###############################################################################
#In the above plots the peak of curve is not fully captured by the curve
#next: try fitting some splines to the data to see if they capture the peak better
###############################################################################
# METHOD 1) from scipy.interpolate import CubicSpline
#NOTE: univariatespline with k=3 gives same curve as np.polyfit with k =3 and scipy.curve_fit with a cubic polynomial
debris_final_increasingx = np.zeros(debris_final.shape)
debris_final_increasingx[:] = debris_final[:]
debris_final_increasingx[4] = debris_final[5]
debris_final_increasingx[5] = debris_final[4]

heightchange_debrisfinal = np.zeros(debris_final.shape)
heightchange_debrisfinal[:] = height_change[:]
heightchange_debrisfinal[4] = height_change[5]
heightchange_debrisfinal[5] = height_change[4]


def univariatespline(debris_list,degree,heightchange=height_change[cluster==1]*(p_ice/1000)):
    x = debris_list[cluster==1]
    x[1] = 0.0001
    y = heightchange
    
    spline = UnivariateSpline(x,y,k=degree) #gives exact same curve as np.polyfit with degree = 3 
    x_spline = np.linspace(0, np.max(x)+0.01, 100)
    y_spline = spline(x_spline)
    
    return x_spline, y_spline #returns debris thickness (cm) vs melt (mw.e.) curves

#try the scipy optimize.curve_fit function
def cubicpolynomial(x,a,b,c,d):
    return a*(x**3) + b*(x**2) + c*x  + d

popt, pope = curve_fit(cubicpolynomial,debris_initial[cluster==1], height_change[cluster==1]*(p_ice/1000))
a,b,c,d = popt

x_curvefit = np.linspace(0, np.max(debris_initial[cluster==1])+0.01, 100)
y_curvefit = cubicpolynomial(x_curvefit,a,b,c,d)

#try andrew's suggestion: scipy.interpolate.splrep
def slprepcurve(debris,k,melt=height_change*(p_ice/1000)):
    knots = np.array(debris[np.where(melt == np.max(melt))])
    tck = splrep(debris[cluster==1], melt[cluster==1],k=k,s=0.01,t=knots)
    xnew = np.linspace(0, np.max(debris[cluster==1])+0.01, 100)
    ynew = splev(xnew,tck)
    
    return xnew, ynew, tck

#PLOT COMPARISON OF ALL CURVE FITTING METHODS    
plt.figure(figsize =(12,6))
plt.subplot(1,3,1)
plt.scatter(debris_initial[cluster==1], height_change[cluster==1]*(p_ice/1000), c='white',edgecolors='k',zorder=10,s=70)
plt.plot(univariatespline(debris_initial,3)[0],univariatespline(debris_initial,3)[1], 'blue',linewidth=1.5)
plt.plot(univariatespline(debris_initial,4)[0],univariatespline(debris_initial,4)[1], 'deeppink',linewidth=1.5)
plt.plot(univariatespline(debris_initial,5)[0],univariatespline(debris_initial,5)[1], 'orange',linewidth=1.5)
#plt.plot(x_curvefit,y_curvefit,linestyle='--')
#plt.plot(x_new,y_new)
plt.plot(slprepcurve(debris_initial,k=2)[0],slprepcurve(debris_initial,k=2)[1],linestyle='--',linewidth=1.5,c='blueviolet')
plt.plot(slprepcurve(debris_initial,k=3)[0],slprepcurve(debris_initial,k=3)[1],linestyle='--',linewidth=1.5,c='darkturquoise')
plt.plot(slprepcurve(debris_initial,k=4)[0],slprepcurve(debris_initial,k=4)[1],linestyle='--',linewidth=1.5,c='k')
#plt.plot(polyline_initdeb,(fit_initdeb(polyline_initdeb))*(p_ice/1000),c='k',linestyle='--')
plt.ylim(1,2.3)
plt.xlim(0,8)
plt.legend(['univariate spline, k=3','univariate spline, k=4','univariate spline, k=5','splrep, k=2','splrep, k=3','splrep, k=4','Data'])
plt.title('2022 Field Data (July thicknesses)')
plt.xlabel('Debris Thickness (cm)',fontsize=12)
plt.ylabel('Melt (m w.e.)',fontsize=12)

plt.subplot(1,3,2)
plt.scatter(debris_average[cluster==1], height_change[cluster==1]*(p_ice/1000), c='white',edgecolors='k',zorder=10,s=70)
plt.plot(univariatespline(debris_average,3)[0],univariatespline(debris_average,3)[1], 'blue',linewidth=1.5)
plt.plot(univariatespline(debris_average,4)[0],univariatespline(debris_average,4)[1], 'deeppink',linewidth=1.5)
plt.plot(univariatespline(debris_average,5)[0],univariatespline(debris_average,5)[1], 'orange',linewidth=1.5)
#plt.plot(polyline_initdeb,(fit_initdeb(polyline_initdeb))*(p_ice/1000),c='k',linestyle='--')
plt.plot(slprepcurve(debris_average,k=2)[0],slprepcurve(debris_average,k=2)[1],linestyle='--',linewidth=1.5,c='blueviolet')
plt.plot(slprepcurve(debris_average,k=3)[0],slprepcurve(debris_average,k=3)[1],linestyle='--',linewidth=1.5,c='darkturquoise')
plt.plot(slprepcurve(debris_average,k=4)[0],slprepcurve(debris_average,k=4)[1],linestyle='--',linewidth=1.5,c='k')
plt.ylim(0.5,2.5)
plt.xlim(0,6)
plt.legend(['univariate spline, k=3','univariate spline, k=4','univariate spline, k=5','splrep, k=2','splrep, k=3','splrep, k=4','Data'])
plt.title('Average debris thickness')
plt.xlabel('Debris Thickness (cm)',fontsize=12)
plt.ylabel('Melt (m w.e.)',fontsize=12)

plt.subplot(1,3,3)
newheights = heightchange_debrisfinal[cluster==1]*(p_ice/1000)
newheights2 = heightchange_debrisfinal*(p_ice/1000) #no cluster index needed for slprep func
plt.scatter(debris_final_increasingx[cluster==1], heightchange_debrisfinal[cluster==1]*(p_ice/1000), c='white',edgecolors='k',zorder=10,s=70)
plt.plot(univariatespline(debris_final_increasingx,3,newheights)[0],univariatespline(debris_final_increasingx,3,newheights)[1], 'blue',linewidth=1.5)
plt.plot(univariatespline(debris_final_increasingx,4,newheights)[0],univariatespline(debris_final_increasingx,4,newheights)[1],'deeppink',linewidth=1.5)
plt.plot(univariatespline(debris_final_increasingx,5,newheights)[0],univariatespline(debris_final_increasingx,5,newheights)[1], 'orange',linewidth=1.5)
#plt.plot(polyline_initdeb,(fit_initdeb(polyline_initdeb))*(p_ice/1000),c='k',linestyle='--')
plt.plot(slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[0],slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[1],linestyle='--',linewidth=1.5,c='blueviolet')
plt.plot(slprepcurve(debris_final_increasingx,k=3,melt=newheights2)[0],slprepcurve(debris_final_increasingx,k=3,melt=newheights2)[1],linestyle='--',linewidth=1.5,c='darkturquoise')
plt.plot(slprepcurve(debris_final_increasingx,k=4,melt=newheights2)[0],slprepcurve(debris_final_increasingx,k=4,melt=newheights2)[1],linestyle='--',linewidth=1.5,c='k')
plt.ylim(1,2.5)
plt.xlim(0,4)
plt.legend(['univariate spline, k=3','univariate spline, k=4','univariate spline, k=5','splrep, k=2','splrep, k=3','splrep, k=4','Data'])
plt.title('2022 Field Data (Aug thicknesses)')
plt.xlabel('Debris Thickness (cm)',fontsize=12)
plt.ylabel('Melt (m w.e.)',fontsize=12)
plt.tight_layout()
#plt.savefig('splinecomparisons.png',bbox_inches = 'tight')

#CALCULATE METRICS TO DETERMINE WHICH CURVE HAS BEST FIT
#calculating the adjusted r^2 value based on this formula: https://www.statology.org/curve-fitting-python/
def adjR_univariate(debris, k, heightchange = height_change[cluster==1]):
    results = {}
    x = debris[cluster==1]
    x[1] = 0.0001
    y = heightchange*(p_ice/1000)
    
    model = UnivariateSpline(x,y,k=k)
    yhat = model(x)

    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y - ybar)**2)
    results['r_squared'] = 1- (((1-(ssreg/sstot))*(len(y)-1))/(len(y)-k-1))

    return results

def adjR_splrep(debris, k, heightchange = height_change):
    results = {}
    x = debris
    y = heightchange*(p_ice/1000)
    tck = (slprepcurve(x,k=2,melt=y))[2]
    
    yhat = splev(x[cluster==1],tck)

    ybar = np.sum(y[cluster==1])/len(y[cluster==1])
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y[cluster==1] - ybar)**2)
    results['r_squared'] = 1- (((1-(ssreg/sstot))*(len(y[cluster==1])-1))/(len(y[cluster==1])-k-1))

    return results

#initial debris
for k in 3,4,5:
    result = adjR_univariate(debris_initial,k)
    print('initial debris... univariate ... k=' +str(k) + ' ... ' + str(result))
for k in 2,3,4:
    result = adjR_splrep(debris_initial,k)
    print('initial debris... splrep ....... k=' +str(k) + ' ... ' + str(result))
print('\n')
#avg debris
for k in 3,4,5:
    result = adjR_univariate(debris_average,k)
    print('avg debris... univariate ... k=' +str(k) + ' ... ' + str(result))
for k in 2,3,4:
    result = adjR_splrep(debris_average,k)
    print('avg debris... splrep ....... k=' +str(k) + ' ... ' + str(result))
print('\n')
#final debris
for k in 3,4,5:
    result = adjR_univariate(debris_final_increasingx,k,heightchange=heightchange_debrisfinal[cluster==1])
    print('final debris... univariate ... k=' +str(k) + ' ... ' + str(result))
for k in 2,3,4:
    result = adjR_splrep(debris_final_increasingx,k,heightchange=heightchange_debrisfinal)
    print('final debris... splrep ....... k=' +str(k) + ' ... ' + str(result))

#CONCLUSION: MOVE FORWARD WITH THE SPLREP, K =2 CURVE FOR ALL 3 INSTANCES OF DEBRIS
    
#re-plot the ostrem curves with the chosen spline;
plt.figure(figsize=(10,7))
plt.title('Debris Thickness vs. Ablation (July 19 - Aug 31 2022) \n data fitted with spline of k = 2',fontsize=14)
plt.plot(slprepcurve(debris_initial,k=2)[0],slprepcurve(debris_initial,k=2)[1],linestyle='-',linewidth=1.5,c='royalblue')
plt.plot(slprepcurve(debris_average,k=2)[0],slprepcurve(debris_average,k=2)[1],linestyle='-',linewidth=1.5,c='orange')
plt.plot(slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[0],slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[1],linestyle='-',linewidth=1.5,c='red')

plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='red',fmt="o",capsize=5)
plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=initial_debris_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=average_debris_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=final_debris_uncertainty[cluster==1],c='red',fmt="o",capsize=5)

plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='red',fmt="o",capsize=5,mfc='white',mec='red')
plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=initial_debris_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=average_debris_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=final_debris_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
#plt.legend(['Initial debris thickness (July)','Average debris thickness','Final debris thickness (Aug)'],fontsize=12)
plt.legend(handles=[blue_line,orange_line,red_line,white_dot],fontsize=14)
plt.xlabel('Debris Thickness (cm)',fontsize=14)
plt.ylabel('Melt (m w.e.)',fontsize=14)
#plt.xlim(0,15)
plt.margins(x=0.01)
#plt.savefig('stakegarden_data_k2splines_mwe.png',bbox_inches = 'tight')

#plot the ostrem curves with the splines that VISUALLY seem like the best fit:
plt.figure(figsize=(10,7))
plt.title('Debris Thickness vs. Ablation (July 19 - Aug 31 2022) \n visual best fit splines',fontsize=14)
plt.plot(univariatespline(debris_initial,5)[0],univariatespline(debris_initial,5)[1], 'royalblue',linewidth=1.5)
plt.plot(slprepcurve(debris_average,k=3)[0],slprepcurve(debris_average,k=3)[1],linestyle='-',linewidth=1.5,c='orange')
plt.plot(slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[0],slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[1],linestyle='-',linewidth=1.5,c='red')

plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]*(p_ice/1000),stake_uncertainty[cluster==1]*(p_ice/1000),c='red',fmt="o",capsize=5)
plt.errorbar(debris_initial[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=initial_debris_uncertainty[cluster==1],c='royalblue',fmt="o",capsize=5)
plt.errorbar(debris_average[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=average_debris_uncertainty[cluster==1],c='orange',fmt="o",capsize=5)
plt.errorbar(debris_final[cluster==1],height_change[cluster==1]*(p_ice/1000),yerr=None,xerr=final_debris_uncertainty[cluster==1],c='red',fmt="o",capsize=5)

plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]*(p_ice/1000),stake_uncertainty[cluster==2]*(p_ice/1000),c='red',fmt="o",capsize=5,mfc='white',mec='red')
plt.errorbar(debris_initial[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=initial_debris_uncertainty[cluster==2],c='royalblue',fmt="o",capsize=5,mfc='white',mec='royalblue')
plt.errorbar(debris_average[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=average_debris_uncertainty[cluster==2],c='orange',fmt="o",capsize=5,mfc='white',mec='orange')
plt.errorbar(debris_final[cluster==2],height_change[cluster==2]*(p_ice/1000),yerr=None,xerr=final_debris_uncertainty[cluster==2],c='red',fmt="o",capsize=5,mfc='white',mec='red')
#plt.legend(['Initial debris thickness (July)','Average debris thickness','Final debris thickness (Aug)'],fontsize=12)
plt.legend(handles=[blue_line,orange_line,red_line,white_dot],fontsize=14)
plt.xlabel('Debris Thickness (cm)',fontsize=14)
plt.ylabel('Melt (m w.e.)',fontsize=14)
#plt.xlim(0,15)
plt.margins(x=0.01)
#plt.savefig('stakegarden_data_bestfitsplines_mwe.png',bbox_inches = 'tight')

###############################################################################
###############################################################################
######### CALCULATE REFERENCE CURVE BASED ON PDD WEIGHTED AVERAGE #############
###############################################################################
###############################################################################

#Step 1. Calculate the true "average"/"reference" curve based on a P.D.D. weighted average

#import temperature data from outpost station for july19-aug31 2022
outpostJulA = np.genfromtxt('D:/FieldData/2022/July/Campbell/Data/Outpost/24July2022_A/TOA5_5216.FiveMin.dat', skip_header=1, delimiter=',')
outpostJulB = np.genfromtxt('D:/FieldData/2022/July/Campbell/Data/Outpost/24July2022_B/TOA5_5216.FiveMin.dat', skip_header=1, delimiter=',')
outpostAug = np.genfromtxt('D:/FieldData/2022/August/Campbell/Data/Outpost/TOA5_5216.FiveMin.dat', skip_header=1, delimiter=',')

airtemp_julA = outpostJulA[3:,6] #starts in Jul 2021, need to cut so that it starts Jul 2022
airtemp_julB = outpostJulB[3:,6]
airtemp_aug = outpostAug[3:,6]
#no missing data WITHIN the arrays, but there may be gaps in time BETWEEN the arrays so need to be careful of that 

#timestamps are showing up as NaN here, need to manually create timestamps with datetime (?)
dates_julA = pd.date_range(start="2021-07-19 17:40:00",end="2022-07-24 17:20:00",freq='5min')
dates_julB = pd.date_range(start="2022-07-24 17:25:00",end="2022-07-24 17:35:00",freq='5min') #5 min time gap between julB and aug records
dates_aug = pd.date_range(start="2022-07-24 17:40:00",end="2022-09-01 10:45:00",freq='5min')

#append the missing time with temp = nan to the jul b record, then concatenate the whole thing
np.append(airtemp_julB,np.nan)
fulltemprecord = np.concatenate((airtemp_julA,airtemp_julB,airtemp_aug))
alldates = np.concatenate((dates_julA,dates_julB,dates_aug))

#cut the temp and dates into the exact days needed (jul 19 to aug 31)
#jul19midnight 
jul19 = np.where(alldates == np.datetime64(dt.datetime(2022,7,19)))
#aug31 11:55pm
aug31 = np.where(alldates == np.datetime64(dt.datetime(2022,8,31,23,55)))

fielddates = alldates[jul19[0][0]:aug31[0][0]+1]
fieldtemps = fulltemprecord[jul19[0][0]:aug31[0][0]+1]


#get daily mean temperatures
#convert the data to a dataframe
df_5min = pd.DataFrame(data={'timestamp': fielddates, 'temp': fieldtemps})
df_daily = df_5min.groupby(by=pd.Grouper(freq='D', key='timestamp')).mean()
#print(df_daily)
arr = np.array(df_daily)
dailytemp = arr[:,0]
#make datetime array to match the daily temps
aws_days = pd.date_range(start="2022-07-19",end="2022-08-31",freq='D')

#plot daily temps:
#plot full 5 min temp record
plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.title('Outpost AWS 5 min record',fontsize=12)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=7))
plt.plot(fielddates,fieldtemps,'k')
plt.gcf().autofmt_xdate()
plt.ylim(-5,30)
plt.ylabel('Temperature (\u00B0 C)',fontsize=12)
plt.subplot(1,2,2)
plt.title('Outpost AWS Daily Averages',fontsize=12)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=7))
plt.plot(aws_days,dailytemp,'k')
plt.gcf().autofmt_xdate()
plt.ylim(-5,30)
plt.ylabel('Temperature (\u00B0 C)',fontsize=12)
#plt.savefig('OutpostAWSTemp.png',bbox_inches = 'tight')

#calculate PDD in first half of time period vs 2nd half 
firsthalf_days = aws_days[0:22]
secondhalf_days = aws_days[22:]
firsthalf_PDD = np.sum(dailytemp[0:22]) #no days below zero so can just do a straight sum
secondhalf_PDD = np.sum(dailytemp[22:])

firsthalfPDD_cumu = np.cumsum(dailytemp[0:22])
secondhalfPDD_cumu = np.cumsum(dailytemp[22:])

#calculate the PDD weighted average debris thicknesses:
debris_PDDaverage = np.zeros(len(debris_final))
for i in range(0,len(debris_final)):
    avg = ((debris_initial[i]*firsthalf_PDD) + (debris_final[i]*secondhalf_PDD))/(firsthalf_PDD+secondhalf_PDD)
    debris_PDDaverage[i] = avg
    

#FINAL CURVES USED TO DEFINE THE RELATIONSHIP B/W PEAK
deb_init_x, deb_init_y = univariatespline(debris_initial,5)[0], univariatespline(debris_initial,5)[1]
deb_avg_x, deb_avg_y = slprepcurve(debris_PDDaverage,k=3)[0],slprepcurve(debris_PDDaverage,k=3)[1] 
deb_fin_x, deb_fin_y = slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[0],slprepcurve(debris_final_increasingx,k=2,melt=newheights2)[1]

#plt.plot(deb_init_x,deb_init_y)
#plt.plot(deb_avg_x,deb_avg_y)
#plt.plot(deb_fin_x,deb_fin_y)

#get reference values (peak melt thickness and transition thickness):
peakmelt_ref = np.where(deb_avg_y == np.max(deb_avg_y))
peakmelt_thickness_ref = deb_avg_x[peakmelt_ref]

peakmelt_i = np.where(deb_init_y == np.max(deb_init_y))
peakmelt_thickness_i = deb_init_x[peakmelt_i]

peakmelt_f = np.where(deb_fin_y == np.max(deb_fin_y))
peakmelt_thickness_f = deb_fin_x[peakmelt_f]

print('Jul obs: debris thickness resulting in peak melt is: ' + str(peakmelt_thickness_i))
print('PDD weighted average debris thickness resulting in peak melt is: ' + str(peakmelt_thickness_ref))
print('Aug obs: debris thickness resulting in peak melt is: ' + str(peakmelt_thickness_f))

#transition thickness:
cleanicemelt = height_change[0]*(p_ice/1000)

diff_from_CI00_i = []
for i in deb_init_y:
    diff_from_CI00_i.append(i-cleanicemelt)
diff_from_cleanice_i = np.abs(diff_from_CI00_i)    
CI00_eq_melt_i = np.where(diff_from_cleanice_i == np.min(diff_from_cleanice_i[peakmelt_i[0][0]:]))
transition_thickness_i = deb_init_x[CI00_eq_melt_i]

diff_from_CI00_ref = []
for i in deb_avg_y:
    diff_from_CI00_ref.append(i-cleanicemelt)
diff_from_cleanice_ref = np.abs(diff_from_CI00_ref)    
CI00_eq_melt_ref = np.where(diff_from_cleanice_ref == np.min(diff_from_cleanice_ref[peakmelt_i[0][0]:]))
transition_thickness_ref = deb_avg_x[CI00_eq_melt_ref]

diff_from_CI00_f = []
for i in deb_fin_y:
    diff_from_CI00_f.append(i-cleanicemelt)
diff_from_cleanice_f = np.abs(diff_from_CI00_f)    
CI00_eq_melt_f = np.where(diff_from_cleanice_f == np.min(diff_from_cleanice_f[peakmelt_i[0][0]:]))
transition_thickness_f = deb_fin_x[CI00_eq_melt_f]

print('Jul obs: transition thickness = ' + str(transition_thickness_i))
print('PDD weighted average transition thickness = ' + str(transition_thickness_ref))
print('Aug obs: transition thickness = ' + str(transition_thickness_f))




