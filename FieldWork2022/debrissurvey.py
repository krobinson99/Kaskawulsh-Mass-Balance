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
red_line = mlines.Line2D([], [], color='red', marker='|', linestyle='None',
                          markersize=10, label='Uncertainty\non measured\ndebris thickness')

#ALL DATA POINTS: NO AVERAGING
plt.figure(figsize=(8,8))
plt.errorbar(dr_estimate,debthickness_obs,uncertainty,color='k',fmt="o",ecolor='red',capsize=5)
for x,y,z in zip(DRthickness_where_noiceobs,no_ice_surface_observed,arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='blue',length_includes_head = True)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.title('Rounce et al. (2021) Debris Thickness Estimate vs Observed Debris Thickness',fontsize=14)
plt.legend(handles=[red_line,blue_arrow,],loc='lower left', bbox_to_anchor=(0.85,0.7),fontsize=12)
plt.savefig('KW_debris_survey_scatter.png',bbox_inches = 'tight')