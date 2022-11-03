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
ice_surface_observed = []
DR_thickness_where_iceobs = []
uncertainties_iceobs = []
cellnames_whereiceobs = []
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
    else:
        ice_surface_observed.append(debthickness_obs[i])
        DR_thickness_where_iceobs.append(dr_estimate[i])
        uncertainties_iceobs.append(uncertainty[i])
        cellnames_whereiceobs.append(cellname[i])
        
        

plt.figure()
plt.plot(dr_estimate)
plt.plot(debthickness_obs)

blue_arrow = mlines.Line2D([], [], color='blue', marker="^", linestyle='None',
                          markersize=10, label='Measurements\nwhere ice\nsurface was\nnot reached')

red_line = mlines.Line2D([], [], color='crimson', marker='|', linestyle='None',
                          markersize=10, label='Debris thickness measurement uncertainty')

black_arrow = mlines.Line2D([], [], color='k', marker="^", linestyle='None',
                          markersize=10, label='Measurement\nuncertainties\nwhere ice surface\nwas not reached')

blue_line = mlines.Line2D([], [], color='mediumblue', marker='|', linestyle='None',
                          markersize=10, label='Debris thickness measurement uncertainty')

orange_line = mlines.Line2D([], [], color='darkorange', marker='|', linestyle='None',
                          markersize=10, label='Debris thickness measurement uncertainty')

black_line = mlines.Line2D([], [], color='k', marker='|', linestyle='None',
                          markersize=15, label='Uncertainty on thickness\nmeasurements that did\nreach ice surface')

bluedot = mlines.Line2D([], [], color='mediumblue', marker='.', linestyle='None',
                          markersize=10, label='Transect A')

reddot = mlines.Line2D([], [], color='crimson', marker='.', linestyle='None',
                          markersize=10, label='Transect C')

orangedot = mlines.Line2D([], [], color='darkorange', marker='.', linestyle='None',
                          markersize=10, label='Transect B')

blackstar = mlines.Line2D([], [], color='k', marker='*', linestyle='None',
                          markersize=10, label='Measurements\nwhere ice surface\nwas not reached')


#ALL DATA POINTS: NO AVERAGING
plt.figure(figsize=(8,8))
plt.errorbar(dr_estimate,debthickness_obs,uncertainty,color='k',fmt="o",ecolor='red',capsize=5,zorder=10,markersize=8)
for x,y,z in zip(DRthickness_where_noiceobs,no_ice_surface_observed,arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='blue',length_includes_head = True,zorder=5)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.scatter(DRthickness_where_noiceobs,no_ice_surface_observed,color='w',marker='*',zorder=15)
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.title('Rounce et al. (2021) Debris Thickness Estimate vs Observed Debris Thickness',fontsize=14)
plt.legend(handles=[red_line,blue_arrow,blackstar],loc='lower left', bbox_to_anchor=(0.85,0.7),fontsize=12)
plt.savefig('KW_debris_survey_scatter.png',bbox_inches = 'tight')


#TRANSECT A,B,C SEPERATELY
A_no_ice_surface_observed = []
A_DRthickness_where_noiceobs = []
A_arrow_height = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice[0:29])):
    if reachedice[0:29][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty[0:29][i] = np.nan
        A_no_ice_surface_observed.append(debthickness_obs[0:29][i])
        A_DRthickness_where_noiceobs.append(dr_estimate[0:29][i])
        arrowheight = 90 - debthickness_obs[0:29][i]
        A_arrow_height.append(arrowheight)
        
B_no_ice_surface_observed = []
B_DRthickness_where_noiceobs = []
B_arrow_height = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice[29:51])):
    if reachedice[29:51][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty[29:51][i] = np.nan
        B_no_ice_surface_observed.append(debthickness_obs[29:51][i])
        B_DRthickness_where_noiceobs.append(dr_estimate[29:51][i])
        arrowheight = 90 - debthickness_obs[29:51][i]
        B_arrow_height.append(arrowheight)

C_no_ice_surface_observed = []
C_DRthickness_where_noiceobs = []
C_arrow_height = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice[51:])):
    if reachedice[51:][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty[51:][i] = np.nan
        C_no_ice_surface_observed.append(debthickness_obs[51:][i])
        C_DRthickness_where_noiceobs.append(dr_estimate[51:][i])
        arrowheight = 90 - debthickness_obs[51:][i]
        C_arrow_height.append(arrowheight)

plt.figure(figsize=(20,6))
plt.subplot(1,3,1)
plt.errorbar(dr_estimate[0:29],debthickness_obs[0:29],uncertainty[0:29],color='k',fmt="o",ecolor='mediumblue',capsize=5,zorder=10)
for x,y,z in zip(A_DRthickness_where_noiceobs,A_no_ice_surface_observed,A_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='mediumblue',length_includes_head = True,zorder=5)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.scatter(A_DRthickness_where_noiceobs,A_no_ice_surface_observed,color='w',marker='*',zorder=15)
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.title('Transect A',fontsize=14)

plt.subplot(1,3,2)
plt.errorbar(dr_estimate[29:51],debthickness_obs[29:51],uncertainty[29:51],color='k',fmt="o",ecolor='darkorange',capsize=5,zorder=10)
for x,y,z in zip(B_DRthickness_where_noiceobs,B_no_ice_surface_observed,B_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='darkorange',length_includes_head = True,zorder=5)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.scatter(B_DRthickness_where_noiceobs,B_no_ice_surface_observed,color='w',marker='*',zorder=15)
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.title('Transect B',fontsize=14)

plt.subplot(1,3,3)
plt.errorbar(dr_estimate[51:],debthickness_obs[51:],uncertainty[51:],color='k',fmt="o",ecolor='crimson',capsize=5,zorder=10)
for x,y,z in zip(C_DRthickness_where_noiceobs,C_no_ice_surface_observed,C_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='crimson',length_includes_head = True,zorder=5)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.scatter(C_DRthickness_where_noiceobs,C_no_ice_surface_observed,color='w',marker='*',zorder=15)
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,35)
plt.ylim(0,35)
plt.title('Transect C',fontsize=14)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
#plt.legend(handles=[black_arrow,blackstar,black_line],loc='lower center', bbox_to_anchor=(0.85,0.7),fontsize=12)
plt.figlegend(handles=[black_arrow,blackstar,black_line],loc='center right',fontsize=12)
plt.tight_layout()
#plt.savefig('KW_debris_survey_3transects.png',bbox_inches = 'tight')


plt.figure(figsize=(8,8))
plt.errorbar(dr_estimate[0:29],debthickness_obs[0:29],uncertainty[0:29],color='mediumblue',fmt="o",ecolor='mediumblue',capsize=7,zorder=10)
plt.errorbar(dr_estimate[29:51],debthickness_obs[29:51],uncertainty[29:51],color='darkorange',fmt="o",ecolor='darkorange',capsize=7,zorder=10)
plt.errorbar(dr_estimate[51:],debthickness_obs[51:],uncertainty[51:],color='crimson',fmt="o",ecolor='crimson',capsize=7,zorder=10)
for x,y,z in zip(A_DRthickness_where_noiceobs,A_no_ice_surface_observed,A_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='mediumblue',length_includes_head = True,zorder=5)
for x,y,z in zip(B_DRthickness_where_noiceobs,B_no_ice_surface_observed,B_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='darkorange',length_includes_head = True,zorder=5)
for x,y,z in zip(C_DRthickness_where_noiceobs,C_no_ice_surface_observed,C_arrow_height):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='crimson',length_includes_head = True,zorder=5)
plt.scatter(DRthickness_where_noiceobs,no_ice_surface_observed,color='w',marker='*',zorder=15)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.title('Rounce et al. (2021) Debris Thickness Estimate vs Observed Debris Thickness',fontsize=14)
plt.legend(handles=[black_arrow,blackstar,bluedot,orangedot,reddot],loc='lower left', bbox_to_anchor=(0.85,0.7),fontsize=12)
#plt.savefig('KW_debris_survey_scatter_coloured.png',bbox_inches = 'tight')



#CALCULATE CORRELATION COEFFICIENTS AMONGST THE DATA SETS
print('All data correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs,ice_surface_observed)[0]))
print('Transect A correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs[0:21],ice_surface_observed[0:21])[0]))
print('Transect B correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs[21:29],ice_surface_observed[21:29])[0]))
print('Transect C correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs[29:],ice_surface_observed[29:])[0]))














# AVERAGED DATA
debrisdata_avg = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/KaskawulshDebris_AveragedData_NoHeaders.csv'

df = pd.read_csv(debrisdata_avg)
arr = np.array(df)


cellname_avg = arr[:,0]
debthickness_obs_avg = arr[:,1] # cm
uncertainty_avg = arr[:,2] # cm
reachedice_avg = arr[:,3] # boolean 0 (false) or 1 (true). If 0, ice surface was not reached, therefore true deb thickness is greater than the recorded thickness
dr_estimate_avg = arr[:,4] #debris thickness estimate for the grid cell from Rounce et al. (2021)

#replace uncertainties where we did NOT reach the ice surface with 0, and append them to a new list for creating arrows in the plot
ice_surface_observed_avg = []
DR_thickness_where_iceobs_avg = []
uncertainties_iceobs_avg = []
cellnames_whereiceobs_avg = []
no_ice_surface_observed_avg = []
DRthickness_where_noiceobs_avg = []
arrow_height_avg = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice_avg)):
    if reachedice_avg[i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty_avg[i] = np.nan
        no_ice_surface_observed_avg.append(debthickness_obs_avg[i])
        DRthickness_where_noiceobs_avg.append(dr_estimate_avg[i])
        arrowheight_avg = 90 - debthickness_obs_avg[i]
        arrow_height_avg.append(arrowheight_avg)
    else:
        ice_surface_observed_avg.append(debthickness_obs_avg[i])
        DR_thickness_where_iceobs_avg.append(dr_estimate_avg[i])
        uncertainties_iceobs_avg.append(uncertainty_avg[i])
        cellnames_whereiceobs_avg.append(cellname_avg[i])
        

#TRANSECT A,B,C SEPERATELY
A_no_ice_surface_observed_avg = []
A_DRthickness_where_noiceobs_avg = []
A_arrow_height_avg = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice_avg[0:18])):
    if reachedice_avg[0:18][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty_avg[0:18][i] = np.nan
        A_no_ice_surface_observed_avg.append(debthickness_obs_avg[0:18][i])
        A_DRthickness_where_noiceobs_avg.append(dr_estimate_avg[0:18][i])
        arrowheight_avg = 90 - debthickness_obs_avg[0:18][i]
        A_arrow_height_avg.append(arrowheight_avg)
        
B_no_ice_surface_observed_avg = []
B_DRthickness_where_noiceobs_avg = []
B_arrow_height_avg = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice_avg[18:38])):
    if reachedice_avg[18:38][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty_avg[18:38][i] = np.nan
        B_no_ice_surface_observed_avg.append(debthickness_obs_avg[18:38][i])
        B_DRthickness_where_noiceobs_avg.append(dr_estimate_avg[18:38][i])
        arrowheight_avg = 90 - debthickness_obs_avg[18:38][i]
        B_arrow_height_avg.append(arrowheight_avg)

C_no_ice_surface_observed_avg = []
C_DRthickness_where_noiceobs_avg = []
C_arrow_height_avg = [] #height needed for the arrow to reach top of grid (y = 90)
for i in range(0,len(reachedice_avg[51:])):
    if reachedice_avg[51:][i] == False: #aka if that measurement did NOT reach the ice..
        uncertainty_avg[51:][i] = np.nan
        C_no_ice_surface_observed_avg.append(debthickness_obs_avg[38:][i])
        C_DRthickness_where_noiceobs_avg.append(dr_estimate_avg[38:][i])
        arrowheight_avg = 90 - debthickness_obs_avg[38:][i]
        C_arrow_height_avg.append(arrowheight_avg)

plt.figure(figsize=(20,6))
plt.subplot(1,3,1)
plt.errorbar(dr_estimate_avg[0:18],debthickness_obs_avg[0:18],uncertainty_avg[0:18],color='k',fmt="o",ecolor='mediumblue',capsize=5,zorder=10)
for x,y,z in zip(A_DRthickness_where_noiceobs_avg,A_no_ice_surface_observed_avg,A_arrow_height_avg):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='mediumblue',length_includes_head = True,zorder=5)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.scatter(A_DRthickness_where_noiceobs_avg,A_no_ice_surface_observed_avg,color='w',marker='*',zorder=15)
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.title('Transect A',fontsize=14)

plt.subplot(1,3,2)
plt.errorbar(dr_estimate_avg[18:38],debthickness_obs_avg[18:38],uncertainty_avg[18:38],color='k',fmt="o",ecolor='darkorange',capsize=5,zorder=10)
for x,y,z in zip(B_DRthickness_where_noiceobs_avg,B_no_ice_surface_observed_avg,B_arrow_height_avg):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='darkorange',length_includes_head = True,zorder=5)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.scatter(B_DRthickness_where_noiceobs_avg,B_no_ice_surface_observed_avg,color='w',marker='*',zorder=15)
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,90)
plt.ylim(0,90)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.title('Transect B',fontsize=14)

plt.subplot(1,3,3)
plt.errorbar(dr_estimate_avg[38:],debthickness_obs_avg[38:],uncertainty_avg[38:],color='k',fmt="o",ecolor='crimson',capsize=5,zorder=10)
for x,y,z in zip(C_DRthickness_where_noiceobs_avg,C_no_ice_surface_observed_avg,C_arrow_height_avg):
    plt.arrow(x,y,0,z,width=0.1,head_width=1.5,color='crimson',length_includes_head = True,zorder=5)
#plt.arrow(75,80,0,5,width=0.1,head_width=1.5,color='red')
plt.scatter(C_DRthickness_where_noiceobs_avg,C_no_ice_surface_observed_avg,color='w',marker='*',zorder=15)
plt.xlabel('Rounce et al. (2021) Debris Thickness Estimate (cm)',fontsize=14)
plt.ylabel('Measured debris thickness (cm)',fontsize=14)
plt.xlim(0,35)
plt.ylim(0,35)
plt.title('Transect C',fontsize=14)
plt.plot([0,100],[0,100],color='grey',linestyle='--')
plt.legend(handles=[black_arrow,blackstar,black_line],loc='lower left', bbox_to_anchor=(0.85,0.7),fontsize=12)
#plt.savefig('KW_debris_survey_3transects_averageddata.png',bbox_inches = 'tight')



#CALCULATE CORRELATION COEFFICIENTS AMONGST THE AVERAGED DATA SET
print('All data correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs_avg,ice_surface_observed_avg)[0]))
print('Transect A correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs_avg[0:10],ice_surface_observed_avg[0:10])[0]))
print('Transect B correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs_avg[10:16],ice_surface_observed_avg[10:16])[0]))
print('Transect C correlation coefficient = ' + str(pearsonr(DR_thickness_where_iceobs_avg[16:],ice_surface_observed_avg[16:])[0]))

###############################################################################
#  MAKE FIGURES FOR GRAIN SIZE DISTRIBUTION OF IMPORTED/NATURAL DEBRIS SAMPLES
###############################################################################

particlediameter = np.array([75, 37.5, 19, 12.5, 4.75, 2, 0.85, 0.425, 0.25, 0.15, 0.075, 0]) #units = mm (from GrainSizeAnalysisKW_Imported.xls)
ppassing_imported = np.array([100.0, 100.0, 100.0, 99.6, 97.4, 92.7, 69.2, 10.0, 1.2, 0.5, 0.3, 0.0])
ppassing_natural = np.array([100.0, 100.0, 97.7, 89.7, 49.6, 15.1, 5.3, 3.2, 2.4, 1.8, 1.0, 0.0])

plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(particlediameter,ppassing_imported)
plt.scatter(particlediameter,ppassing_imported)
plt.xlabel('Particle Diameter (mm)',fontsize=12)
plt.ylabel('% Finer by Weight',fontsize=12)
plt.xlim(0,30)
plt.subplot(1,2,2)
plt.plot(particlediameter,ppassing_natural)
plt.scatter(particlediameter,ppassing_natural)
plt.xlabel('Particle Diameter (mm)',fontsize=12)
plt.ylabel('% Finer by Weight',fontsize=12)
plt.xlim(0,30)


plt.figure(figsize=(6,6))
plt.title('Grain Size Distribution \n for Debris from Kaskawulsh Ablation Experiment',fontsize=12)
plt.plot(particlediameter,ppassing_imported)
plt.scatter(particlediameter,ppassing_imported)
plt.xlabel('Particle Diameter (mm)',fontsize=12)
plt.ylabel('% Finer by Weight',fontsize=12)
plt.xlim(0,30)
plt.plot(particlediameter,ppassing_natural)
plt.scatter(particlediameter,ppassing_natural)
plt.xlabel('Particle Diameter (mm)',fontsize=12)
plt.ylabel('% Finer by Weight',fontsize=12)
plt.grid()
plt.legend(['DB01/DB02/DB03/DB04 Debris','MR07/KV04 Debris'],fontsize=12)
#plt.savefig('KWgrainsizedistribution.png',bbox_inches = 'tight')

###############################################################################
#  MAKE FIGURES FOR DEBRIS ALBEDO (data from Matthew Sturm at UAlaska - in ASD SAND folder)
###############################################################################
albedodata = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/ASD SAND 8-11-22/Sand_Albedo_Plottable.csv'

df = pd.read_csv(albedodata)
arr = np.array(df)

wavelength = arr[:,0]
drysand1 = arr[:,1]
drysand2 = arr[:,2]
drysand3 = arr[:,3]
wetsand1 = arr[:,4]
wetsand2 = arr[:,5]

plt.figure(figsize=(8,6))
plt.plot(wavelength,drysand1)
plt.plot(wavelength,drysand2)
plt.plot(wavelength,drysand3)
plt.plot(wavelength,wetsand1)
plt.plot(wavelength,wetsand2)
plt.ylim(0,1)
plt.margins(x=0)
plt.grid()
plt.xlabel('Wavelength (nm)',fontsize=14)
plt.ylabel('Albedo',fontsize=14)
plt.legend(['Dry Debris 1','Dry Debris 2','Dry Debris 3','Wet Debris 1','Wet Debris 2'],fontsize=12)
plt.savefig('debrisalbedo.png',bbox_inches = 'tight')


