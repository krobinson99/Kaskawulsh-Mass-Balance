# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 17:49:20 2023

Compared REF_MODEL to ROUNCE DEBRIS model outputs.

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
distributed_average_mass_balance, massbalance_timeseries_12years, \
massbalance_timeseries_average, snowrunoff_timeseries_12years, snowrunoff_timeseries_average, \
annualrunoff_stackedbar, date_of_zero_balance, runoff_timeseries_average_discharge_withstddev, calculate_stddev_timeseries, \
runoff_timeseries_average_SLRformat, distributed_runoff, distributed_glaciericemelt, distributed_snowrunoff, distributed_rainrunoff, \
distributed_SImelt

# Functions to plot the differences between two models:
from MBM_plottingfunctions import distributed_mass_balance_difference

model_name = 'ROUNCE_DEBRIS'
MODEL_OUTPUTS_REFMODEL = 'D:/Model Runs/REF_MODEL'
MODEL_OUTPUTS_ROUNCEDEBRIS = 'D:/Model Runs/ROUNCE_DEBRIS'

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
snowdepth_kw, Ptau_kw, SI_kw, temp_kw = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

massbal_kw_std, totalsnowmelt_kw_std, refreezing_kw_std, netsnowmelt_kw_std, gl_icemelt_kw_std, superimp_icemelt_kw_std, rain_kw_std, refrozen_rain_kw_std, rain_runoff_kw_std, accumulation_kw_std, \
snowdepth_kw_std, Ptau_kw_std, SI_kw_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

# =============================================================================
#  REF MODEL: KW AVERAGE
# =============================================================================
massbal_kw_ref, totalsnowmelt_kw_ref, refreezing_kw_ref, netsnowmelt_kw_ref, gl_icemelt_kw_ref, superimp_icemelt_kw_ref, rain_kw_ref, refrozen_rain_kw_ref, rain_runoff_kw_ref, accumulation_kw_ref, \
snowdepth_kw_ref, Ptau_kw_ref, SI_kw_ref, temp_kw_ref = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

massbal_kw_ref_std, totalsnowmelt_kw_ref_std, refreezing_kw_ref_std, netsnowmelt_kw_ref_std, gl_icemelt_kw_ref_std, superimp_icemelt_kw_ref_std, rain_kw_ref_std, refrozen_rain_kw_ref_std, rain_runoff_kw_ref_std, accumulation_kw_ref_std, \
snowdepth_kw_ref_std, Ptau_kw_ref_std, SI_kw_ref_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)



# =============================================================================
#  Get timeseries (catchment average) 
# =============================================================================
massbal_krh, totalsnowmelt_krh, refreezing_krh, netsnowmelt_krh, gl_icemelt_krh, superimp_icemelt_krh, rain_krh, refrozen_rain_krh, rain_runoff_krh, accumulation_krh, \
snowdepth_krh, Ptau_krh, SI_krh, temp_krh = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

massbal_krh_std, totalsnowmelt_krh_std, refreezing_krh_std, netsnowmelt_krh_std, gl_icemelt_krh_std, superimp_icemelt_krh_std, rain_krh_std, refrozen_rain_krh_std, rain_runoff_krh_std, accumulation_krh_std, \
snowdepth_krh_std, Ptau_krh_std, SI_krh_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

# =============================================================================
#  Get timeseries (all glacierized area) 
# =============================================================================
massbal_allgl, totalsnowmelt_allgl, refreezing_allgl, netsnowmelt_allgl, gl_icemelt_allgl, superimp_icemelt_allgl, rain_allgl, refrozen_rain_allgl, rain_runoff_allgl, accumulation_allgl, \
snowdepth_allgl, Ptau_allgl, SI_allgl, temp_allgl = calculate_mb_components_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

massbal_allgl_std, totalsnowmelt_allgl_std, refreezing_allgl_std, netsnowmelt_allgl_std, gl_icemelt_allgl_std, superimp_icemelt_allgl_std, rain_allgl_std, refrozen_rain_allgl_std, rain_runoff_allgl_std, accumulation_allgl_std, \
snowdepth_allgl_std, Ptau_allgl_std, SI_allgl_std = calculate_stddev_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

# =============================================================================
#  Distributed totals
# =============================================================================
massbal_dist, totalsnowmelt_dist, refreezing_dist, netsnowmelt_dist, gl_icemelt_dist, superimp_icemelt_dist, rain_dist, refrozen_rain_dist, rain_runoff_dist, accumulation_dist, \
snowdepth_dist, Ptau_dist, SI_dist, temp_dist = calculate_mb_components_distributed(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID)






# =============================================================================
#  REF_MODEL Get timeseries (catchment average) 
# =============================================================================
massbal_krh_ref, totalsnowmelt_krh_ref, refreezing_krh_ref, netsnowmelt_krh_ref, gl_icemelt_krh_ref, superimp_icemelt_krh_ref, rain_krh_ref, refrozen_rain_krh_ref, rain_runoff_krh_ref, accumulation_krh_ref, \
snowdepth_krh_ref, Ptau_krh_ref, SI_krh_ref, temp_krh_ref = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

massbal_krh_ref_std, totalsnowmelt_krh_ref_std, refreezing_krh_ref_std, netsnowmelt_krh_ref_std, gl_icemelt_krh_ref_std, superimp_icemelt_krh_ref_std, rain_krh_ref_std, refrozen_rain_krh_ref_std, rain_runoff_krh_ref_std, accumulation_krh_ref_std, \
snowdepth_krh_ref_std, Ptau_krh_ref_std, SI_krh_ref_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

# =============================================================================
#   REF_MODEL Get timeseries (all glacierized area) 
# =============================================================================
massbal_allgl_ref, totalsnowmelt_allgl_ref, refreezing_allgl_ref, netsnowmelt_allgl_ref, gl_icemelt_allgl_ref, superimp_icemelt_allgl_ref, rain_allgl_ref, refrozen_rain_allgl_ref, rain_runoff_allgl_ref, accumulation_allgl_ref, \
snowdepth_allgl_ref, Ptau_allgl_ref, SI_allgl_ref, temp_allgl_ref = calculate_mb_components_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

massbal_allgl_ref_std, totalsnowmelt_allgl_ref_std, refreezing_allgl_ref_std, netsnowmelt_allgl_ref_std, gl_icemelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, rain_allgl_ref_std, refrozen_rain_allgl_ref_std, rain_runoff_allgl_ref_std, accumulation_allgl_ref_std, \
snowdepth_allgl_ref_std, Ptau_allgl_ref_std, SI_allgl_ref_std = calculate_stddev_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

# =============================================================================
#  REF_MODEL Distributed totals
# =============================================================================
massbal_dist_ref, totalsnowmelt_dist_ref, refreezing_dist_ref, netsnowmelt_dist_ref, gl_icemelt_dist_ref, superimp_icemelt_dist_ref, rain_dist_ref, refrozen_rain_dist_ref, rain_runoff_dist_ref, accumulation_dist_ref, \
snowdepth_dist_ref, Ptau_dist_ref, SI_dist_ref, temp_dist_ref = calculate_mb_components_distributed(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID)


# =============================================================================
# DISTRIBUTED PLOTS
# =============================================================================

shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.888283378746594,stop=1,name='massbal')
distributed_average_mass_balance(model_name,np.arange(2007,2017+1),years,massbal_dist,np.linspace(-16.3,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
#plt.savefig(os.path.join(MODEL_OUTPUTS,'2007-2018_massbalance.png'),bbox_inches='tight')

shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.21052631578947367,stop=1,name='massbal_diff')
distributed_mass_balance_difference(model_name,np.arange(1980,2021+1),years,massbal_dist_ref,massbal_dist,np.linspace(-2,7.5,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_runoff(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

def distributed_runoff_difference(model_name,avg_years,all_years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,netsnowmelt_dist_alt,gl_icemelt_dist_alt,superimp_icemelt_dist_alt,rain_runoff_dist_alt,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY_REF = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_REF[year-avg_years[0],:,:] = netsnowmelt_dist_ref[year-all_years[0]] + gl_icemelt_dist_ref[year-all_years[0]] + superimp_icemelt_dist_ref[year-all_years[0]] + rain_runoff_dist_ref[year-all_years[0]]
    MB_REF = np.nanmean(MB_ARRAY_REF,axis=0)
        
    MB_ARRAY_ALT = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_ALT[year-avg_years[0],:,:] = netsnowmelt_dist_alt[year-all_years[0]] + gl_icemelt_dist_alt[year-all_years[0]] + superimp_icemelt_dist_alt[year-all_years[0]] + rain_runoff_dist_alt[year-all_years[0]]
    MB_ALT = np.nanmean(MB_ARRAY_ALT,axis=0)
    
    DIFF = np.array(MB_REF-MB_ALT)
    DIFF[np.where(DIFF==0)] = np.nan
        
    print(np.nanmin(DIFF),np.nanmax(DIFF))    
    
    plt.figure(figsize=(9,5))
    plt.title('Total Runoff ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='runoff_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,20))
    legend.ax.set_ylabel('Difference (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.tight_layout()

#shiftedColorMap(matplotlib.cm.BrBG,start=0,midpoint=0.7894736842105263,stop=1,name='runoff_diff')
shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.79,stop=1,name='runoff_diff')
distributed_runoff_difference(model_name,np.arange(1980,2021+1),years,netsnowmelt_dist_ref,\
                              gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,\
                              netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,\
                              np.linspace(-7.5,2,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)




distributed_glaciericemelt(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_snowrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,0.7,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_rainrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,0.3,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_SImelt(np.arange(1980,2021+1),years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,np.linspace(0,0.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)



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
