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

# Import functions for model:
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import write_config_file, save_to_netcdf, shiftedColorMap
from Model_functions_ver4 import get_meanSP, Calculate_Pmean, MassBalance, max_superimposed_ice

sys.path.insert(2,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\ProcessOutputs')
from MBM_plottingfunctions import load_hydrologic_year, calculate_mb_components_timeseries, \
calculate_mb_components_distributed, runoff_timeseries_12years, runoff_piecharts_12years, \
runoff_timeseries_average, distributed_average_mass_balance, massbalance_timeseries_12years, \
massbalance_timeseries_average, snowrunoff_timeseries_12years, snowrunoff_timeseries_average

MODEL_OUTPUTS = 'D:/TuningOutputs/Passed_StageI_All_runs_averaged'
NARR_INPUTS = 'D:/BiasCorrected_files/KRH' 
sim = 9999

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
massbal_kw, totalsnowmelt_kw, refreezing_kw, netsnowmelt_kw, superimp_icemelt_kw, gl_icemelt_kw, rain_runoff_kw, refrozen_rain_kw, rain_kw, accumulation_kw \
= calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

dist_massbal, dist_totalsnowmelt, dist_refrozenmelt, dist_netsnowmelt, dist_superimp_icemelt, dist_gl_icemelt, dist_rain_runoff, dist_refrozen_rain, dist_rain, dist_accumulation \
= calculate_mb_components_distributed(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)

# =============================================================================
#  Get timeseries (catchment average) 
# =============================================================================
massbal_krh, totalsnowmelt_krh, refreezing_krh, netsnowmelt_krh, superimp_icemelt_krh, gl_icemelt_krh, rain_runoff_krh, refrozen_rain_krh, rain_krh, accumulation_krh \
= calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)


# =============================================================================
#  Get timeseries (all glacierized area) 
# =============================================================================
massbal_allgl, totalsnowmelt_allgl, refreezing_allgl, netsnowmelt_allgl, superimp_icemelt_allgl, gl_icemelt_allgl, rain_runoff_allgl, refrozen_rain_allgl, rain_allgl, accumulation_allgl \
= calculate_mb_components_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)



# =============================================================================
#  Get timeseries in the ablation zone
# =============================================================================
massbal_abl, totalsnowmelt_abl, refreezing_abl, netsnowmelt_abl, superimp_icemelt_abl, gl_icemelt_abl, rain_runoff_abl, refrozen_rain_abl, rain_abl, accumulation_abl \
= calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,57,216,PointScale=True,KWaverage=False,Catchmentaverage=False)


# =============================================================================
#  Get timeseries in the acc zone
# =============================================================================
massbal_acc, totalsnowmelt_acc, refreezing_acc, netsnowmelt_acc, superimp_icemelt_acc, gl_icemelt_acc, rain_runoff_acc, refrozen_rain_acc, rain_acc, accumulation_acc \
= calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,210,231,PointScale=True,KWaverage=False,Catchmentaverage=False)


# =============================================================================
#  Get timeseries near the ela
# =============================================================================
massbal_ela, totalsnowmelt_ela, refreezing_ela, netsnowmelt_ela, superimp_icemelt_ela, gl_icemelt_ela, rain_runoff_ela, refrozen_rain_ela, rain_ela, accumulation_ela \
= calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,87,67,PointScale=True,KWaverage=False,Catchmentaverage=False)


# =============================================================================
#  KW wide average plots: 
# =============================================================================
runoff_timeseries_12years('Glacier-wide average runoff',2007,years,0.03,2.2,gl_icemelt_kw,netsnowmelt_kw, rain_runoff_kw, superimp_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

runoff_timeseries_average('Glacier-wide average runoff: ',np.arange(2007,2018+1),years,0.01,1.6,gl_icemelt_kw, netsnowmelt_kw, rain_runoff_kw, superimp_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_average_runoff_2007-2018.png'),bbox_inches='tight')

runoff_piecharts_12years('Glacier-wide average',2007,years,gl_icemelt_kw,netsnowmelt_kw, rain_runoff_kw, superimp_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_12years('Glacier-wide average mass balance',2007, years, 0.055, 1.5, accumulation_kw, refrozen_rain_kw ,netsnowmelt_kw, superimp_icemelt_kw, gl_icemelt_kw)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('Glacier-wide average mass balance: ',np.arange(2007,2018+1),years,0.025,0.8, accumulation_kw, refrozen_rain_kw, netsnowmelt_kw, superimp_icemelt_kw, gl_icemelt_kw)

snowrunoff_timeseries_12years('Glacier-wide average snow melt runoff',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_kw, refreezing_kw)

snowrunoff_timeseries_average('Glacier-wide average snow melt runoff: ',np.arange(2007,2018+1),years,-0.005, 0.01, -0.233333, 0.7, totalsnowmelt_kw, refreezing_kw)


# =============================================================================
#  Catchment wide average plots: 
# =============================================================================
runoff_timeseries_12years('Catchment-wide average runoff',2007,years,0.03,2,gl_icemelt_krh,netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

runoff_timeseries_average('Catchment-wide average runoff',np.arange(2007,2018+1),years,0.01,1.4,gl_icemelt_krh, netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_average_runoff_2007-2018.png'),bbox_inches='tight')

runoff_piecharts_12years('Catchment-wide average',2007,years,gl_icemelt_krh,netsnowmelt_krh, rain_runoff_krh, superimp_icemelt_krh)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_12years('Catchment-wide average mass balance',2007, years, 0.055, 1.5, accumulation_krh, refrozen_rain_krh ,netsnowmelt_krh, superimp_icemelt_krh, gl_icemelt_krh)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('Catchment-wide average mass balance: ',np.arange(2007,2018+1),years,0.025,0.8, accumulation_krh, refrozen_rain_krh, netsnowmelt_krh, superimp_icemelt_krh, gl_icemelt_krh)

snowrunoff_timeseries_12years('Catchment-wide average snow melt runoff',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_krh, refreezing_krh)

snowrunoff_timeseries_average('Catchment-wide average snow melt runoff: ',np.arange(2007,2018+1),years,-0.005, 0.015, -0.233333, 0.7, totalsnowmelt_krh, refreezing_krh)


# =============================================================================
#  Ablation zone plots: 
# =============================================================================
runoff_timeseries_12years('KW_M',2007,years,0.03,2,gl_icemelt_abl,netsnowmelt_abl, rain_runoff_abl, superimp_icemelt_abl)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

runoff_timeseries_average('KW_M: ',np.arange(2007,2018+1),years,0.06,4,gl_icemelt_abl, netsnowmelt_abl, rain_runoff_abl, superimp_icemelt_abl)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_average_runoff_2007-2018.png'),bbox_inches='tight')

runoff_piecharts_12years('KW_M',2007,years,gl_icemelt_abl,netsnowmelt_abl, rain_runoff_abl, superimp_icemelt_abl)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_12years('KW_M',2007, years, 0.055, 1.5, accumulation_abl, refrozen_rain_abl ,netsnowmelt_abl, superimp_icemelt_abl, gl_icemelt_abl)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('KW_M: ',np.arange(2007,2018+1),years,0.025,0.8, accumulation_abl, refrozen_rain_abl, netsnowmelt_abl, superimp_icemelt_abl, gl_icemelt_abl)

snowrunoff_timeseries_12years('KW_M',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_abl, refreezing_abl)

snowrunoff_timeseries_average('KW_M: ',np.arange(2007,2018+1),years,-0.005, 0.015, -0.233333, 0.7, totalsnowmelt_abl, refreezing_abl)



# =============================================================================
#  Accumulation zone plots: 
# =============================================================================
runoff_timeseries_12years('SA2670',2007,years,0.03,2,gl_icemelt_acc,netsnowmelt_acc, rain_runoff_acc, superimp_icemelt_acc)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

runoff_timeseries_average('SA2670',np.arange(2007,2018+1),years,0.06,6,gl_icemelt_acc, netsnowmelt_acc, rain_runoff_acc, superimp_icemelt_acc)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_average_runoff_2007-2018.png'),bbox_inches='tight')

runoff_piecharts_12years('SA2670',2007,years,gl_icemelt_acc,netsnowmelt_acc, rain_runoff_acc, superimp_icemelt_acc)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_12years('SA2670',2007, years, 0.055, 1.5, accumulation_acc, refrozen_rain_acc ,netsnowmelt_acc, superimp_icemelt_acc, gl_icemelt_acc)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('SA2670: ',np.arange(2007,2018+1),years,0.025,0.8, accumulation_acc, refrozen_rain_acc, netsnowmelt_acc, superimp_icemelt_acc, gl_icemelt_acc)

snowrunoff_timeseries_12years('SA2670',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_acc, refreezing_acc)

snowrunoff_timeseries_average('SA2670: ',np.arange(2007,2018+1),years,-0.005, 0.015, -0.233333, 0.7, totalsnowmelt_acc, refreezing_acc)

# =============================================================================
#  ELA zone plots: 
# =============================================================================
runoff_timeseries_12years('CA2400',2007,years,0.03,2,gl_icemelt_ela,netsnowmelt_ela, rain_runoff_ela, superimp_icemelt_ela)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_runoff_2007-2018.png'),bbox_inches='tight')

runoff_timeseries_average('CA2400: ',np.arange(2007,2018+1),years,0.06,4,gl_icemelt_ela, netsnowmelt_ela, rain_runoff_ela, superimp_icemelt_ela)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_average_runoff_2007-2018.png'),bbox_inches='tight')

runoff_piecharts_12years('CA2400',2007,years,gl_icemelt_ela,netsnowmelt_ela, rain_runoff_ela, superimp_icemelt_ela)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_piecharts_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_12years('CA2400',2007, years, 0.055, 1.5, accumulation_ela, refrozen_rain_ela ,netsnowmelt_ela, superimp_icemelt_ela, gl_icemelt_ela)
#fig.savefig(os.path.join(MODEL_OUTPUTS,'KW_massbalance_2007-2018.png'),bbox_inches='tight')

massbalance_timeseries_average('CA2400: ',np.arange(2007,2018+1),years,0.025,0.8, accumulation_ela, refrozen_rain_ela, netsnowmelt_ela, superimp_icemelt_ela, gl_icemelt_ela)

snowrunoff_timeseries_12years('CA2400',2007, years, -0.015, 0.035, -0.3428571428571429, 0.8, totalsnowmelt_ela, refreezing_ela)

snowrunoff_timeseries_average('CA2400: ',np.arange(2007,2018+1),years,-0.005, 0.015, -0.233333, 0.7, totalsnowmelt_ela, refreezing_ela)





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
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.7734806629834254,stop=1,name='massbal')
distributed_average_mass_balance(np.arange(2007,2018+1),years,dist_massbal,np.linspace(-7,2.05,17),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
plt.scatter(Xgrid[57,216],Ygrid[57,216],s=100,c='k')
plt.scatter(Xgrid[210,231],Ygrid[210,231],s=100,c='k')
plt.scatter(Xgrid[87,67],Ygrid[87,67],s=100,c='k')
#plt.savefig(os.path.join(MODEL_OUTPUTS,'2007-2018_massbalance.png'),bbox_inches='tight')



# =============================================================================
# Annual total runoff, stacked bar plot 
# =============================================================================

# get total runoff from each component for all years:
yrly_snowmelt = []
yrly_icemelt = []
yrly_SImelt = []
yrly_rain = []
for year in years:
    i = year-years[0]
    # get total of each component
    yrly_snowmelt.append(np.sum(netsnowmelt_kw[i]))
    yrly_icemelt.append(np.sum(gl_icemelt_kw[i]))
    yrly_SImelt.append(np.sum(superimp_icemelt_kw[i]))
    yrly_rain.append(np.sum(rain_runoff_kw[i]))

total_runoff = np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt) + np.array(yrly_icemelt)
m, b = np.polyfit(years,total_runoff,deg=1)

width = 0.75       # the width of the bars: can also be len(x) sequence

fig, ax = plt.subplots(figsize=(12,5))
ax.bar(years, yrly_rain, width, label='Rain',color='deeppink')
ax.bar(years, yrly_SImelt, width, bottom=yrly_rain,label='Superimposed ice melt',color='darkorange')
ax.bar(years, yrly_snowmelt, width, bottom=(np.array(yrly_rain) + np.array(yrly_SImelt)),label='Snow melt',color='royalblue')
ax.bar(years, yrly_icemelt, width, bottom=(np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt)),label='Glacier ice melt',color='turquoise')
plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2)

ax.set_ylabel('Runoff (m w.e. yr$^{-1}$)',fontsize=14)
plt.xticks(np.arange(1980,2022,5),fontsize=14,rotation=45)
plt.yticks(fontsize=14)
#plt.ylim(0,2.5)
#ax.set_title('Monthly Runoff (2007-2018) \n Debris Case',fontsize = 15)
ax.legend(fontsize=14)
plt.title('Glacier-wide total annual runoff',fontsize=14)

plt.tight_layout()    


# transition date
transition_dates = []
transition_DOYs = []
for year in years:
    dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
    massbal = (accumulation_kw[year-years[0]] + refrozen_rain_kw[year-years[0]]) - (netsnowmelt_kw[year-years[0]] + superimp_icemelt_kw[year-years[0]] + gl_icemelt_kw[year-years[0]])
    
    d_tr = [np.where(np.cumsum(massbal)[50:] <=0)][0][0]
    if len(d_tr) == 0:
        transition_date = 'N/A'
        transition_DOY = np.nan
    else:
        transition_date = (dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]])
        transition_DOY = np.where(np.cumsum(massbal)[50:] <=0)[0][0]
    
    transition_dates.append(transition_date)
    transition_DOYs.append(transition_DOY+50)
    
monthday = []
for i in np.arange(min(transition_DOYs),max(transition_DOYs)+1,10):
    monthday.append( str(dates[i])[5:10])
    
    
m, b = np.polyfit(years[np.isfinite(transition_DOYs)],np.array(transition_DOYs)[np.isfinite(transition_DOYs)],deg=1)

plt.figure(figsize=(8,4))
plt.title('Date of zero balance',fontsize=14)
plt.plot(years,transition_DOYs,c='darkblue')
plt.yticks(np.arange(min(transition_DOYs),max(transition_DOYs)+1,10),monthday,fontsize=14)
plt.xticks(fontsize=14)
plt.ylabel('date of $\dot{B}$ = 0',fontsize=14)
plt.grid()
plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2)

# annual non renewable component of runoff (exceeds previous years accumulation)
yrly_snowmelt = []
yrly_icemelt = []
yrly_SImelt = []
yrly_rain = []
for year in years:
    i = year-years[0]
    # get total of each component
    yrly_snowmelt.append(np.sum(netsnowmelt_kw[i]))
    yrly_icemelt.append(np.sum(gl_icemelt_kw[i]))
    yrly_SImelt.append(np.sum(superimp_icemelt_kw[i]))
    yrly_rain.append(np.sum(rain_runoff_kw[i]))


# total annual glacier contribution to total catchment runoff:
