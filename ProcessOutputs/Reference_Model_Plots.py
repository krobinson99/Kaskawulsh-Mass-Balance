# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 15:09:24 2023

@author: katierobinson
"""

import numpy as np
from netCDF4 import Dataset
import pandas as pd
import matplotlib
import cmocean
import matplotlib.pyplot as plt
import sys
import os
import random

# Import functions for model:
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import write_config_file, save_to_netcdf, shiftedColorMap
from Model_functions_ver4 import get_meanSP, Calculate_Pmean, MassBalance, max_superimposed_ice

sys.path.insert(2,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\ProcessOutputs')
from MBM_plottingfunctions import load_hydrologic_year, calculate_mb_components_timeseries, \
calculate_mb_components_distributed, runoff_timeseries_12years, runoff_piecharts_12years, \
distributed_average_mass_balance, massbalance_timeseries_12years, \
massbalance_timeseries_average, snowrunoff_timeseries_12years, snowrunoff_timeseries_average, \
annualrunoff_stackedbar, date_of_zero_balance, runoff_timeseries_average_discharge_withstddev, calculate_stddev_timeseries, \
runoff_timeseries_average_SLRformat, distributed_runoff

MODEL_OUTPUTS = 'D:/Model Runs/REF_MODEL/Sim_99999_v2'
NARR_INPUTS = 'D:/BiasCorrected_files/KRH' 
sim = 99999

years = np.arange(1979,2021+1)
Glacier_ID = 'KRH'
R2S = 1
Precip_inputs = 'D:/BiasCorrected_files/KRH'                                 # Path to folder where downscaled & bias-corrected data is.
Temp_inputs = 'D:/BiasCorrected_files/KRH'

# Input geometry
Easting_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Xgrid.txt' # Paths to text files defining Easting/Northing coords of every model gridcell
Northing_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Ygrid.txt'
Sfc_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_SfcType.txt'      # Path to text file where 1 = off-glacier, 0 = on-glacier, NaN = not in the domain.

Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
Sfc = np.loadtxt(Sfc_grid)

Catchmentoutline = np.array(Sfc)
Catchmentoutline[np.where(np.isfinite(Sfc))] = 0
Catchmentoutline[np.where(np.isnan(Sfc))] = 1

KRH_tributaries = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Tributaries.txt')

All_glacierized_area = np.array(Sfc)
All_glacierized_area[np.where(Sfc==1)] = np.nan

# =============================================================================
#  Get timeseries (Kaskawulsh average) and distributed totals for all years:
# =============================================================================
massbal_kw, totalsnowmelt_kw, refreezing_kw, netsnowmelt_kw, gl_icemelt_kw, superimp_icemelt_kw, rain_kw, refrozen_rain_kw, rain_runoff_kw, accumulation_kw, \
snowdepth_kw, Ptau_kw, SI_kw, temp_kw = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

massbal_dist, totalsnowmelt_dist, refreezing_dist, netsnowmelt_dist, gl_icemelt_dist, superimp_icemelt_dist, rain_dist, refrozen_rain_dist, rain_runoff_dist, accumulation_dist, \
snowdepth_dist, Ptau_dist, SI_dist, temp_dist = calculate_mb_components_distributed(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)

# =============================================================================
#  Get timeseries (catchment average) 
# =============================================================================
massbal_krh, totalsnowmelt_krh, refreezing_krh, netsnowmelt_krh, gl_icemelt_krh, superimp_icemelt_krh, rain_krh, refrozen_rain_krh, rain_runoff_krh, accumulation_krh, \
snowdepth_krh, Ptau_krh, SI_krh, temp_krh = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

massbal_krh_std, totalsnowmelt_krh_std, refreezing_krh_std, netsnowmelt_krh_std, gl_icemelt_krh_std, superimp_icemelt_krh_std, rain_krh_std, refrozen_rain_krh_std, rain_runoff_krh_std, accumulation_krh_std, \
snowdepth_krh_std, Ptau_krh_std, SI_krh_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

# =============================================================================
#  Get timeseries (all glacierized area) 
# =============================================================================
massbal_allgl, totalsnowmelt_allgl, refreezing_allgl, netsnowmelt_allgl, gl_icemelt_allgl, superimp_icemelt_allgl, rain_allgl, refrozen_rain_allgl, rain_runoff_allgl, accumulation_allgl, \
snowdepth_allgl, Ptau_allgl, SI_allgl, temp_allgl = calculate_mb_components_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

massbal_allgl_std, totalsnowmelt_allgl_std, refreezing_allgl_std, netsnowmelt_allgl_std, gl_icemelt_allgl_std, superimp_icemelt_allgl_std, rain_allgl_std, refrozen_rain_allgl_std, rain_runoff_allgl_std, accumulation_allgl_std, \
snowdepth_allgl_std, Ptau_allgl_std, SI_allgl_std = calculate_stddev_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

# =============================================================================
#  Get timeseries at specific points:
# Ablation zone: 57, 216
# Accumulation xone: 210, 231
# ELA: 87,67
# =============================================================================

# =============================================================================
#  KW wide average plots: 
# =============================================================================
runoff_timeseries_12years('Glacier-wide average runoff',2010,years,0.031,2.2,gl_icemelt_kw,netsnowmelt_kw, rain_runoff_kw, superimp_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

runoff_timeseries_average_discharge_withstddev('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,250,2,gl_icemelt_krh, netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh,gl_icemelt_std, netsnowmelt_std, rain_runoff_std, superimp_icemelt_std,Sfc)

runoff_piecharts_12years('Glacier-wide average',2007,years,gl_icemelt_kw,netsnowmelt_kw, rain_runoff_kw, superimp_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_12years('Glacier-wide average mass balance',2007, years, 0.055, 1.5, accumulation_kw, refrozen_rain_kw ,netsnowmelt_kw, superimp_icemelt_kw, gl_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('Glacier-wide average mass balance: ',np.arange(1980,2021+1),years,0.025,0.8, accumulation_kw, refrozen_rain_kw, netsnowmelt_kw, superimp_icemelt_kw, gl_icemelt_kw)

snowrunoff_timeseries_12years('Glacier-wide average snow melt runoff',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_kw, refreezing_kw)

snowrunoff_timeseries_average('Glacier-wide average snow melt runoff: ',np.arange(2007,2018+1),years,-0.005, 0.01, -0.233333, 0.7, totalsnowmelt_kw, refreezing_kw)

annualrunoff_stackedbar('Glacier-wide total annual runoff',years,netsnowmelt_kw,gl_icemelt_kw,superimp_icemelt_kw,rain_runoff_kw)
    
date_of_zero_balance(years,massbal_kw)

# =============================================================================
#  Catchment wide average plots: 
# =============================================================================

#RUNOFF

runoff_timeseries_12years('Catchment-wide average runoff',2007,years,0.03,2,gl_icemelt_krh,netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

# avg runoff with standard deviations:
runoff_timeseries_average_discharge_withstddev('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,300,4,gl_icemelt_krh, netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh,gl_icemelt_krh_std, netsnowmelt_krh_std, rain_runoff_krh_std, superimp_icemelt_krh_std,Sfc)

runoff_timeseries_average_SLRformat('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_krh, netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh,gl_icemelt_krh_std, netsnowmelt_krh_std, rain_runoff_krh_std, superimp_icemelt_krh_std,Sfc)

runoff_piecharts_12years('Catchment-wide average',2007,years,gl_icemelt_krh,netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

snowrunoff_timeseries_12years('Catchment-wide average snow melt runoff',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_krh, refreezing_krh)

snowrunoff_timeseries_average('Catchment-wide average snow melt runoff: ',np.arange(2007,2018+1),years,-0.005, 0.015, -0.233333, 0.7, totalsnowmelt_krh, refreezing_krh)

annualrunoff_stackedbar('Catchment-wide annual runoff',years[1:],netsnowmelt_krh,gl_icemelt_krh,superimp_icemelt_krh,rain_runoff_krh,'discharge',np.where(np.isfinite(Sfc))[0].shape[0]*(200*200)/1e9)
    

#MASS BALANCE
massbalance_timeseries_12years('Catchment-wide average mass balance',2007, years, 0.055, 1.5, accumulation_krh, refrozen_rain_krh ,netsnowmelt_krh, superimp_icemelt_krh, gl_icemelt_krh)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('Catchment-wide average mass balance: ',np.arange(2007,2017+1),years,0.025,0.8, accumulation_krh, refrozen_rain_krh, netsnowmelt_krh, superimp_icemelt_krh, gl_icemelt_krh, accumulation_krh_std, refrozen_rain_krh_std, netsnowmelt_krh_std, superimp_icemelt_krh_std, gl_icemelt_krh_std)

# =============================================================================
#  All glacierized area plots: 
# =============================================================================
runoff_timeseries_12years('Glacier-wide average runoff',2010,years,0.031,2.2,gl_icemelt_kw,netsnowmelt_kw, rain_runoff_kw, superimp_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

runoff_timeseries_average_discharge_withstddev('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,300,4,gl_icemelt_allgl, netsnowmelt_allgl, rain_runoff_allgl, superimp_icemelt_allgl, gl_icemelt_allgl_std, netsnowmelt_allgl_std, rain_runoff_allgl_std, superimp_icemelt_allgl_std,All_glacierized_area)

runoff_timeseries_average_SLRformat('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_allgl, netsnowmelt_allgl, rain_runoff_allgl, superimp_icemelt_allgl, gl_icemelt_allgl_std, netsnowmelt_allgl_std, rain_runoff_allgl_std, superimp_icemelt_allgl_std,All_glacierized_area)


runoff_piecharts_12years('Glacier-wide average',2007,years,gl_icemelt_kw,netsnowmelt_kw, rain_runoff_kw, superimp_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_12years('Glacier-wide average mass balance',2007, years, 0.055, 1.5, accumulation_kw, refrozen_rain_kw ,netsnowmelt_kw, superimp_icemelt_kw, gl_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('Glacier-wide average mass balance: ',np.arange(2007,2017+1),years,0.025,0.8, accumulation_allgl, refrozen_rain_allgl, netsnowmelt_allgl, superimp_icemelt_allgl, gl_icemelt_allgl, accumulation_allgl_std, refrozen_rain_allgl_std, netsnowmelt_allgl_std, superimp_icemelt_allgl_std, gl_icemelt_allgl_std)

snowrunoff_timeseries_12years('Glacier-wide average snow melt runoff',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_kw, refreezing_kw)

snowrunoff_timeseries_average('Glacier-wide average snow melt runoff: ',np.arange(2007,2018+1),years,-0.005, 0.01, -0.233333, 0.7, totalsnowmelt_kw, refreezing_kw)

annualrunoff_stackedbar('Glacier-wide total annual runoff',years,netsnowmelt_kw,gl_icemelt_kw,superimp_icemelt_kw,rain_runoff_kw)
    
date_of_zero_balance(years,massbal_allgl)

sim = 25
time=np.arange(0,365)
plt.plot(time,np.cumsum(gl_icemelt_allgl[sim]))
plt.fill_between(time,np.cumsum(gl_icemelt_allgl[sim]-gl_icemelt_allgl_std[sim]),np.cumsum(gl_icemelt_allgl[sim]+gl_icemelt_allgl_std[sim]),alpha=0.3)

plt.plot(time,np.cumsum(netsnowmelt_allgl[sim]))
plt.fill_between(time,np.cumsum(netsnowmelt_allgl[sim]-netsnowmelt_allgl_std[sim]),np.cumsum(netsnowmelt_allgl[sim]+netsnowmelt_allgl_std[sim]),alpha=0.3)

plt.plot(time,np.cumsum(superimp_icemelt_allgl[sim]))
plt.fill_between(time,np.cumsum(superimp_icemelt_allgl[sim]-superimp_icemelt_allgl_std[sim]),np.cumsum(superimp_icemelt_allgl[sim]+superimp_icemelt_allgl_std[sim]),alpha=0.3)

plt.plot(time,np.cumsum(rain_runoff_allgl[sim]))
plt.fill_between(time,np.cumsum(rain_runoff_allgl[sim]-rain_runoff_allgl_std[sim]),np.cumsum(rain_runoff_allgl[sim]+rain_runoff_allgl_std[sim]),alpha=0.3)

plt.grid()
plt.ylim(0,1.7)
# distributed runoff components:
# =============================================================================
# ice_sum = np.zeros((Sfc.shape))
# snow_sum = np.zeros((Sfc.shape))
# rain_sum = np.zeros((Sfc.shape))
# SI_sum = np.zeros((Sfc.shape))
# 
# for year in avg_years:
#     i = year - all_years[0]
#     ice_sum += np.array(dist_gl_icemelt[i][:365])
#     snow_sum += np.array(dist_netsnowmelt[i][:365])
#     rain_sum += np.array(dist_rain_runoff[i][:365])
#     SI_sum += np.array(dist_superimp_icemelt[i][:365])
#     
# ice_mean = np.array(ice_sum/len(avg_years))
# snow_mean = np.array(snow_sum/len(avg_years))
# rain_mean = np.array(rain_sum/len(avg_years))
# SI_mean = np.array(SI_sum/len(avg_years))
# 
# plt.figure()
# =============================================================================


# =============================================================================
# Distributed plots:
# =============================================================================
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.8333333333333334,stop=1,name='massbal')
distributed_average_mass_balance(np.arange(2007,2017+1),years,massbal_dist,np.linspace(-10,2.05,17),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
#plt.savefig(os.path.join(MODEL_OUTPUTS,'2007-2018_massbalance.png'),bbox_inches='tight')

distributed_runoff(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,10,17),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_glaciericemelt(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,10,17),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_snowrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,0.7,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_rainrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,0.3,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_SImelt(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,0.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)


# =============================================================================
# Practice uncertainty plots:
# =============================================================================   
noise = random.uniform(0.6, 1.0)
print(noise)

rounce_debris_synthetic = []
for i in range(0,len(massbal_kw[0])):
    noise = random.uniform(0.6, 1.0)
    rounce_debris_synthetic.append(massbal_kw[0][i]*noise)
    
plt.figure()
plt.plot(np.cumsum(massbal_kw[25]))
plt.plot(np.cumsum(rounce_debris_synthetic))

plt.plot(np.cumsum(massbal_kw[25])-np.cumsum(rounce_debris_synthetic)[:-1])

fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')       
ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=14)

total_accumulation = snowfall_mean + refrain_mean
total_ablation = snowmelt_mean + SImelt_mean + gl_melt_mean

ax.plot(np.arange(0,len(dates)),total_accumulation,c='mediumblue',label='Accumulation')    
ax.plot(np.arange(0,len(dates)),-total_ablation,c='red',label='Ablation')    
ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
ax.grid()
ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
ax.tick_params(axis='y',labelsize=14)
ax.margins(x=0)

massbal = total_accumulation - total_ablation

transition_date = dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]]
ax0 = ax.twinx()
ax0.plot(np.arange(0,len(dates)),np.cumsum(massbal),c='k',label='Cumulative balance\n$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
#ax0.plot(np.arange(0,len(dates)),np.cumsum(gl_icemelt[year-years[0]]),c='turquoise',label='cumulative runoff')
#ax0.plot(np.arange(0,len(dates)),np.cumsum(snowmelt[year-years[0]]),c='royalblue',label='cumulative runoff')
#ax0.plot(np.arange(0,len(dates)),np.cumsum(superimp_icemelt[year-years[0]]),c='orange')
ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
ax0.set_ylabel('Cumulative Mass Balance (m w.e.)',fontsize=14)
ax0.margins(x=0)
ax0.tick_params(axis='y',labelsize=14)
ax0.legend(fontsize=14,loc='lower left')
    
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
#plt.legend(by_label.values(), by_label.keys())
ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
fig.tight_layout()





# =============================================================================
# total annual glacier contribution to total catchment runoff:
# =============================================================================   
yrly_snowmelt = []
yrly_icemelt = []
yrly_SImelt = []
yrly_rain = []
for year in years:
    i = year-years[0]
    # get total of each component
    yrly_snowmelt.append(np.nansum(netsnowmelt_dist[i]))
    yrly_icemelt.append(np.nansum(gl_icemelt_dist[i]))
    yrly_SImelt.append(np.nansum(superimp_icemelt_dist[i]))
    yrly_rain.append(np.nansum(rain_runoff_dist[i]))

total_runoff_krh = np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt) + np.array(yrly_icemelt)  

yrly_snowmelt = []
yrly_icemelt = []
yrly_SImelt = []
yrly_rain = []
for year in years:
    i = year-years[0]
    # get total of each component
    yrly_snowmelt.append(np.nansum(netsnowmelt_dist[i][np.where(Sfc==0)]))
    yrly_icemelt.append(np.nansum(gl_icemelt_dist[i][np.where(Sfc==0)]))
    yrly_SImelt.append(np.nansum(superimp_icemelt_dist[i][np.where(Sfc==0)]))
    yrly_rain.append(np.nansum(rain_runoff_dist[i][np.where(Sfc==0)]))

total_runoff_gl = np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt) + np.array(yrly_icemelt)  

gl_fraction = total_runoff_gl/total_runoff_krh*100

m, b = np.polyfit(years,gl_fraction,deg=1)

plt.figure()
plt.plot(years,gl_fraction)
plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Glacier contribution to annual\ncatchment-wide runoff (%)',fontsize=14)
plt.ylim(82.5,90)
plt.text(1980,89,s='a)',fontsize=15,fontweight='bold')
plt.grid()






# =============================================================================
# total Kaskawulsh glacier contribution to total catchment runoff:
# =============================================================================   
yrly_snowmelt = []
yrly_icemelt = []
yrly_SImelt = []
yrly_rain = []
for year in years:
    i = year-years[0]
    # get total of each component
    yrly_snowmelt.append(np.nansum(netsnowmelt_dist[i]))
    yrly_icemelt.append(np.nansum(gl_icemelt_dist[i]))
    yrly_SImelt.append(np.nansum(superimp_icemelt_dist[i]))
    yrly_rain.append(np.nansum(rain_runoff_dist[i]))

total_runoff_krh = np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt) + np.array(yrly_icemelt)  

yrly_snowmelt = []
yrly_icemelt = []
yrly_SImelt = []
yrly_rain = []
for year in years:
    i = year-years[0]
    # get total of each component
    yrly_snowmelt.append(np.nansum(netsnowmelt_dist[i][np.isfinite(KRH_tributaries)]))
    yrly_icemelt.append(np.nansum(gl_icemelt_dist[i][np.isfinite(KRH_tributaries)]))
    yrly_SImelt.append(np.nansum(superimp_icemelt_dist[i][np.isfinite(KRH_tributaries)]))
    yrly_rain.append(np.nansum(rain_runoff_dist[i][np.isfinite(KRH_tributaries)]))

total_runoff_gl = np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt) + np.array(yrly_icemelt)  

gl_fraction = total_runoff_gl/total_runoff_krh*100

m, b = np.polyfit(years,gl_fraction,deg=1)

plt.figure()
plt.plot(years,gl_fraction)
plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Kaskawulsh contribution to annual\ncatchment-wide runoff (%)',fontsize=14)
plt.ylim(82.5,90)
plt.text(1980,89,s='b)',fontsize=15,fontweight='bold')
plt.grid()


