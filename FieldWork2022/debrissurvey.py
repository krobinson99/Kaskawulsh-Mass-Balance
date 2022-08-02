# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 13:24:35 2022

script to look at the debris transect data from field work completed by
GF, KR, and TH in July 2022

@author: katierobinson
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats.stats import pearsonr  

debrisdata = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/KaskawulshDebris_RawData_NoHeaders.csv'

df = pd.read_csv(debrisdata)
arr = np.array(df)

cellname = arr[:,0]
debthickness_obs = arr[:,1] # cm
uncertainty = arr[:,2] # cm
reachedice = arr[:,3] # boolean 0 (false) or 1 (true). If 0, ice surface was not reached, therefore true deb thickness is greater than the recorded thickness
dr_estimate = arr[:,4] #debris thickness estimate for the grid cell from Rounce et al. (2021)

#replace uncertainties where we did NOT reach the ice surface with 0, and append them to a new list for creating arrows in the plot
no_ice_surface_observed = []
DRthickness_where_noiceobs = []
arrow_height = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice)):
    if reachedice[i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty[i] = np.nan
        no_ice_surface_observed.append(debthickness_obs[i])
        DRthickness_where_noiceobs.append(dr_estimate[i])
        arrowheight = 90 - debthickness_obs[i]
        arrow_height.append(arrowheight)
        
        

plt.figure()
plt.plot(dr_estimate)
plt.plot(debthickness_obs)

blue_arrow = mlines.Line2D([], [], color='blue', marker="^", linestyle='None',
                          markersize=10, label='Measurements\nwhere ice\nsurface was\nnot reached')

red_line = mlines.Line2D([], [], color='crimson', marker='|', linestyle='None',
                          markersize=10, label='Debris thickness measurement uncertainty')

black_arrow = mlines.Line2D([], [], color='k', marker="^", linestyle='None',
                          markersize=10, label='Arrows = Measurements\nwhere ice surface\nwas not reached')

blue_line = mlines.Line2D([], [], color='mediumblue', marker='|', linestyle='None',
                          markersize=10, label='Debris thickness measurement uncertainty')

orange_line = mlines.Line2D([], [], color='darkorange', marker='|', linestyle='None',
                          markersize=10, label='Debris thickness measurement uncertainty')

black_line = mlines.Line2D([], [], color='k', marker='|', linestyle='None',
                          markersize=15, label='Error bars =\nUncertainty on thickness\nmeasurements that\nreached ice surface')

bluedot = mlines.Line2D([], [], color='mediumblue', marker='.', linestyle='None',
                          markersize=10, label='Transect A')

reddot = mlines.Line2D([], [], color='crimson', marker='.', linestyle='None',
                          markersize=10, label='Transect C')

orangedot = mlines.Line2D([], [], color='darkorange', marker='.', linestyle='None',
                          markersize=10, label='Transect B')

#ALL DATA POINTS: NO AVERAGING
plt.figure(figsize=(8,8))
plt.errorbar(dr_estimate,debthickness_obs,uncertainty,color='k',fmt="o",ecolor='red',capsize=5)
for x,y,z in zip(DRthickness_where_noiceobs,no_ice_surface_observed,arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='blue',length_includes_head = True)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.title('Rounce et al. (2021) Debris Thickness Estimate vs Observed Debris Thickness',fontsize=14)
plt.legend(handles=[red_line,blue_arrow,],loc='lower left', bbox_to_anchor=(0.85,0.7),fontsize=12)
#plt.savefig('KW_debris_survey_scatter.png',bbox_inches = 'tight')


#TRANSECT A,B,C SEPERATELY
A_no_ice_surface_observed = []
A_DRthickness_where_noiceobs = []
A_arrow_height = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice[0:28])):
    if reachedice[0:28][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty[0:28][i] = np.nan
        A_no_ice_surface_observed.append(debthickness_obs[0:28][i])
        A_DRthickness_where_noiceobs.append(dr_estimate[0:28][i])
        arrowheight = 90 - debthickness_obs[0:28][i]
        A_arrow_height.append(arrowheight)
        
B_no_ice_surface_observed = []
B_DRthickness_where_noiceobs = []
B_arrow_height = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice[29:50])):
    if reachedice[29:50][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty[29:50][i] = np.nan
        B_no_ice_surface_observed.append(debthickness_obs[29:50][i])
        B_DRthickness_where_noiceobs.append(dr_estimate[29:50][i])
        arrowheight = 90 - debthickness_obs[29:50][i]
        B_arrow_height.append(arrowheight)

C_no_ice_surface_observed = []
C_DRthickness_where_noiceobs = []
C_arrow_height = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice[50:])):
    if reachedice[50:][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty[50:][i] = np.nan
        C_no_ice_surface_observed.append(debthickness_obs[50:][i])
        C_DRthickness_where_noiceobs.append(dr_estimate[50:][i])
        arrowheight = 90 - debthickness_obs[50:][i]
        C_arrow_height.append(arrowheight)

plt.figure(figsize=(20,6))
plt.subplot(1,3,1)
plt.errorbar(dr_estimate[0:28],debthickness_obs[0:28],uncertainty[0:28],color='k',fmt="o",ecolor='mediumblue',capsize=5)
for x,y,z in zip(A_DRthickness_where_noiceobs,A_no_ice_surface_observed,A_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='mediumblue',length_includes_head = True)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.title('Transect A',fontsize=14)

plt.subplot(1,3,2)
plt.errorbar(dr_estimate[29:50],debthickness_obs[29:50],uncertainty[29:50],color='k',fmt="o",ecolor='darkorange',capsize=5)
for x,y,z in zip(B_DRthickness_where_noiceobs,B_no_ice_surface_observed,B_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='darkorange',length_includes_head = True)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.title('Transect B',fontsize=14)

plt.subplot(1,3,3)
plt.errorbar(dr_estimate[50:],debthickness_obs[50:],uncertainty[50:],color='k',fmt="o",ecolor='crimson',capsize=5)
for x,y,z in zip(C_DRthickness_where_noiceobs,C_no_ice_surface_observed,C_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='crimson',length_includes_head = True)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,35)
plt.ylim(0,35)
plt.title('Transect C',fontsize=14)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.legend(handles=[black_arrow,black_line],loc='lower left', bbox_to_anchor=(0.85,0.7),fontsize=12)
#plt.savefig('KW_debris_survey_3transects.png',bbox_inches = 'tight')


plt.figure(figsize=(8,8))
plt.errorbar(dr_estimate[0:28],debthickness_obs[0:28],uncertainty[0:28],color='mediumblue',fmt="o",ecolor='mediumblue',capsize=7)
plt.errorbar(dr_estimate[29:50],debthickness_obs[29:50],uncertainty[29:50],color='darkorange',fmt="o",ecolor='darkorange',capsize=7)
plt.errorbar(dr_estimate[50:],debthickness_obs[50:],uncertainty[50:],color='crimson',fmt="o",ecolor='crimson',capsize=7)
for x,y,z in zip(A_DRthickness_where_noiceobs,A_no_ice_surface_observed,A_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='mediumblue',length_includes_head = True)
for x,y,z in zip(B_DRthickness_where_noiceobs,B_no_ice_surface_observed,B_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='darkorange',length_includes_head = True)
for x,y,z in zip(C_DRthickness_where_noiceobs,C_no_ice_surface_observed,C_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='crimson',length_includes_head = True)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.title('Rounce et al. (2021) Debris Thickness Estimate vs Observed Debris Thickness',fontsize=14)
plt.legend(handles=[black_arrow,bluedot,orangedot,reddot],loc='lower left', bbox_to_anchor=(0.85,0.7),fontsize=12)
#plt.savefig('KW_debris_survey_scatter_coloured.png',bbox_inches = 'tight')



#CALCULATE CORRELATION COEFFICIENTS AMONGST THE DATA SETS
print('All data correlation coefficient = ' + str(pearsonr(dr_estimate,debthickness_obs)[0]))
print('Transect A correlation coefficient = ' + str(pearsonr(dr_estimate[0:28],debthickness_obs[0:28])[0]))
print('Transect B correlation coefficient = ' + str(pearsonr(dr_estimate[29:50],debthickness_obs[29:50])[0]))
print('Transect C correlation coefficient = ' + str(pearsonr(dr_estimate[50:],debthickness_obs[50:])[0]))



Alldata_corrcoef = np.corrcoef(np.array(debthickness_obs),np.array(dr_estimate))
TransectA_corrcoef = np.corrcoeff(debthickness_obs,dr_estimate)
