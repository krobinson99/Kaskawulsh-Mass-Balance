# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 15:09:24 2023

@author: katierobinson
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import os


# Import functions for model:
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import shiftedColorMap

sys.path.insert(2,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\ProcessOutputs')
from MBM_plottingfunctions import load_hydrologic_year, calculate_mb_components_timeseries, \
calculate_mb_components_distributed, runoff_timeseries_12years, runoff_piecharts_12years, \
distributed_average_mass_balance, massbalance_timeseries_12years, cumulative_massbalance, \
massbalance_timeseries_average, snowrunoff_timeseries_12years, snowrunoff_timeseries_average, \
annualrunoff_stackedbar, date_of_zero_balance, runoff_timeseries_average_discharge_withstddev, calculate_stddev_timeseries, \
runoff_timeseries_average_SLRformat, distributed_runoff, distributed_glaciericemelt, distributed_snowrunoff, distributed_rainrunoff, \
distributed_SImelt, compare_hydrographs_differences, annual_massbalance_barchart, massbalance_timeseries_comparison, \
massbalance_timeseries_difference, compare_hydrographs, compare_date_of_zero_balance, distributed_mass_balance_difference, \
distributed_runoff_difference, distributed_glaciermelt_difference, distributed_rainrunoff_difference, \
distributed_snowrunoff_difference, distributed_SImelt_difference

REF_MODEL_PATH = 'D:/Model Runs/REF_MODEL/Sim_99999_v2'
ROUNCE_DEBRIS_PATH = 'D:/Model Runs/ROUNCE_DEBRIS'
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

varnames = ['massbal','totalsnowmelt','refreezing','netsnowmelt','gl_icemelt','superimp_icemelt','rain','refrozen_rain','rain_runoff','accumulation','snowdepth','Ptau','SI']

# Make a function to load the daily timeseries for each model output saved as .text files:
def load_daily_timeseries(MODEL_PATH,modelname,varnames,domain,years):
    outputs = []
    
    for var in varnames:
        print(var)
        outputs.append(np.loadtxt(os.path.join(MODEL_PATH, modelname + '_' + var + '_' + domain + '_' + str(years[0]) + '-' + str(years[-1]) + '.txt')))
    
    return outputs

def load_distributed_fields(MODEL_PATH,modelname,varnames,years):
    outputs = []
    
    for var in varnames:
        print(var)
        varlist = []
        for year in years:
            print(year)
            directory = os.path.join(MODEL_PATH,var + '_distributed')
            varlist.append(np.loadtxt(os.path.join(directory,modelname + '_' + var + '_distributed_' + str(year) + '.txt')))
        
        outputs.append(varlist)
            
    return outputs

# =============================================================================
#  REF_MODEL: Catchment-wide average
# =============================================================================
massbal_krh_ref, totalsnowmelt_krh_ref, refreezing_krh_ref, netsnowmelt_krh_ref, gl_icemelt_krh_ref, superimp_icemelt_krh_ref, rain_krh_ref, refrozen_rain_krh_ref, rain_runoff_krh_ref, accumulation_krh_ref, \
snowdepth_krh_ref, Ptau_krh_ref, SI_krh_ref = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'krh',years)

massbal_krh_ref_std, totalsnowmelt_krh_ref_std, refreezing_krh_ref_std, netsnowmelt_krh_ref_std, gl_icemelt_krh_ref_std, superimp_icemelt_krh_ref_std, rain_krh_ref_std, refrozen_rain_krh_ref_std, rain_runoff_krh_ref_std, accumulation_krh_ref_std, \
snowdepth_krh_ref_std, Ptau_krh_ref_std, SI_krh_ref_std = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'krh_std',years)

# =============================================================================
#  REF_MODEL: All glacierized area
# =============================================================================
massbal_allgl_ref, totalsnowmelt_allgl_ref, refreezing_allgl_ref, netsnowmelt_allgl_ref, gl_icemelt_allgl_ref, superimp_icemelt_allgl_ref, rain_allgl_ref, refrozen_rain_allgl_ref, rain_runoff_allgl_ref, accumulation_allgl_ref, \
snowdepth_allgl_ref, Ptau_allgl_ref, SI_allgl_ref = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'allgl',years)

massbal_allgl_ref_std, totalsnowmelt_allgl_ref_std, refreezing_allgl_ref_std, netsnowmelt_allgl_ref_std, gl_icemelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, rain_allgl_ref_std, refrozen_rain_allgl_ref_std, rain_runoff_allgl_ref_std, accumulation_allgl_ref_std, \
snowdepth_allgl_ref_std, Ptau_allgl_ref_std, SI_allgl_ref_std = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'allgl_std',years)

# =============================================================================
#  REF MODEL: Kaskawulsh
# =============================================================================
massbal_kw_ref, totalsnowmelt_kw_ref, refreezing_kw_ref, netsnowmelt_kw_ref, gl_icemelt_kw_ref, superimp_icemelt_kw_ref, rain_kw_ref, refrozen_rain_kw_ref, rain_runoff_kw_ref, accumulation_kw_ref, \
snowdepth_kw_ref, Ptau_kw_ref, SI_kw_ref = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'kw',years)

massbal_kw_ref_std, totalsnowmelt_kw_ref_std, refreezing_kw_ref_std, netsnowmelt_kw_ref_std, gl_icemelt_kw_ref_std, superimp_icemelt_kw_ref_std, rain_kw_ref_std, refrozen_rain_kw_ref_std, rain_runoff_kw_ref_std, accumulation_kw_ref_std, \
snowdepth_kw_ref_std, Ptau_kw_ref_std, SI_kw_ref_std = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'kw_std',years)

# =============================================================================
#  REF MODEL: Debris-covered cell at the terminus
# =============================================================================
massbal_deb_ref, totalsnowmelt_deb_ref, refreezing_deb_ref, netsnowmelt_deb_ref, gl_icemelt_deb_ref, superimp_icemelt_deb_ref, rain_deb_ref, refrozen_rain_deb_ref, rain_runoff_deb_ref, accumulation_deb_ref, \
snowdepth_deb_ref, Ptau_deb_ref, SI_deb_ref = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'deb',years)

massbal_deb_ref_std, totalsnowmelt_deb_ref_std, refreezing_deb_ref_std, netsnowmelt_deb_ref_std, gl_icemelt_deb_ref_std, superimp_icemelt_deb_ref_std, rain_deb_ref_std, refrozen_rain_deb_ref_std, rain_runoff_deb_ref_std, accumulation_deb_ref_std, \
snowdepth_deb_ref_std, Ptau_deb_ref_std, SI_deb_ref_std = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'deb_std',years)

# =============================================================================
#  REF_MODEL: Distributed totals
# =============================================================================
massbal_dist_ref, totalsnowmelt_dist_ref, refreezing_dist_ref, netsnowmelt_dist_ref, gl_icemelt_dist_ref, superimp_icemelt_dist_ref, rain_dist_ref, refrozen_rain_dist_ref, rain_runoff_dist_ref, accumulation_dist_ref, \
snowdepth_dist_ref, Ptau_dist_ref, SI_dist_ref = load_distributed_fields(REF_MODEL_PATH,'REF_MODEL',varnames,years)









# =============================================================================
#  ROUNCE_DEBRIS: Catchment-wide average
# =============================================================================
massbal_krh_rounce, totalsnowmelt_krh_rounce, refreezing_krh_rounce, netsnowmelt_krh_rounce, gl_icemelt_krh_rounce, superimp_icemelt_krh_rounce, rain_krh_rounce, refrozen_rain_krh_rounce, rain_runoff_krh_rounce, accumulation_krh_rounce, \
snowdepth_krh_rounce, Ptau_krh_rounce, SI_krh_rounce = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'krh',years)

massbal_krh_rounce_std, totalsnowmelt_krh_rounce_std, refreezing_krh_rounce_std, netsnowmelt_krh_rounce_std, gl_icemelt_krh_rounce_std, superimp_icemelt_krh_rounce_std, rain_krh_rounce_std, refrozen_rain_krh_rounce_std, rain_runoff_krh_rounce_std, accumulation_krh_rounce_std, \
snowdepth_krh_rounce_std, Ptau_krh_rounce_std, SI_krh_rounce_std = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'krh_std',years)

# =============================================================================
#  ROUNCE_DEBRIS: All glacierized area
# =============================================================================
massbal_allgl_rounce, totalsnowmelt_allgl_rounce, refreezing_allgl_rounce, netsnowmelt_allgl_rounce, gl_icemelt_allgl_rounce, superimp_icemelt_allgl_rounce, rain_allgl_rounce, refrozen_rain_allgl_rounce, rain_runoff_allgl_rounce, accumulation_allgl_rounce, \
snowdepth_allgl_rounce, Ptau_allgl_rounce, SI_allgl_rounce = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'allgl',years)

massbal_allgl_rounce_std, totalsnowmelt_allgl_rounce_std, refreezing_allgl_rounce_std, netsnowmelt_allgl_rounce_std, gl_icemelt_allgl_rounce_std, superimp_icemelt_allgl_rounce_std, rain_allgl_rounce_std, refrozen_rain_allgl_rounce_std, rain_runoff_allgl_rounce_std, accumulation_allgl_rounce_std, \
snowdepth_allgl_rounce_std, Ptau_allgl_rounce_std, SI_allgl_rounce_std = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'allgl_std',years)

# =============================================================================
#  ROUNCE_DEBRIS: Kaskawulsh
# =============================================================================
massbal_kw_rounce, totalsnowmelt_kw_rounce, refreezing_kw_rounce, netsnowmelt_kw_rounce, gl_icemelt_kw_rounce, superimp_icemelt_kw_rounce, rain_kw_rounce, refrozen_rain_kw_rounce, rain_runoff_kw_rounce, accumulation_kw_rounce, \
snowdepth_kw_rounce, Ptau_kw_rounce, SI_kw_rounce = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'kw',years)

massbal_kw_rounce_std, totalsnowmelt_kw_rounce_std, refreezing_kw_rounce_std, netsnowmelt_kw_rounce_std, gl_icemelt_kw_rounce_std, superimp_icemelt_kw_rounce_std, rain_kw_rounce_std, refrozen_rain_kw_rounce_std, rain_runoff_kw_rounce_std, accumulation_kw_rounce_std, \
snowdepth_kw_rounce_std, Ptau_kw_rounce_std, SI_kw_rounce_std = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'kw_std',years)

# =============================================================================
#  ROUNCE_DEBRIS: Debris-covered cell at the terminus
# =============================================================================
massbal_deb_rounce, totalsnowmelt_deb_rounce, refreezing_deb_rounce, netsnowmelt_deb_rounce, gl_icemelt_deb_rounce, superimp_icemelt_deb_rounce, rain_deb_rounce, refrozen_rain_deb_rounce, rain_runoff_deb_rounce, accumulation_deb_rounce, \
snowdepth_deb_rounce, Ptau_deb_rounce, SI_deb_rounce = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'deb',years)

massbal_deb_rounce_std, totalsnowmelt_deb_rounce_std, refreezing_deb_rounce_std, netsnowmelt_deb_rounce_std, gl_icemelt_deb_rounce_std, superimp_icemelt_deb_rounce_std, rain_deb_rounce_std, refrozen_rain_deb_rounce_std, rain_runoff_deb_rounce_std, accumulation_deb_rounce_std, \
snowdepth_deb_rounce_std, Ptau_deb_rounce_std, SI_deb_rounce_std = load_daily_timeseries(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,'deb_std',years)

# =============================================================================
#  ROUNCE_DEBRIS: Distributed totals
# =============================================================================
massbal_dist_rounce, totalsnowmelt_dist_rounce, refreezing_dist_rounce, netsnowmelt_dist_rounce, gl_icemelt_dist_rounce, superimp_icemelt_dist_rounce, rain_dist_rounce, refrozen_rain_dist_rounce, rain_runoff_dist_rounce, accumulation_dist_rounce, \
snowdepth_dist_rounce, Ptau_dist_rounce, SI_dist_rounce = load_distributed_fields(ROUNCE_DEBRIS_PATH,'ROUNCE_DEBRIS',varnames,years)







# =============================================================================
# Reference model plots
# =============================================================================
# Average mass balance
massbalance_timeseries_average('Catchment-wide average mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_krh_ref, refrozen_rain_krh_ref, netsnowmelt_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref, accumulation_krh_ref_std, refrozen_rain_krh_ref_std, netsnowmelt_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_ref_std)
massbalance_timeseries_average('Glacierized area mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_allgl_ref, refrozen_rain_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref, accumulation_allgl_ref_std, refrozen_rain_allgl_ref_std, netsnowmelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_ref_std)
massbalance_timeseries_average('Kaskawulsh mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std)
massbalance_timeseries_average('Terminus mass balance: ',np.arange(2007,2017+1),years,0.15,10, accumulation_deb_ref, refrozen_rain_deb_ref, netsnowmelt_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref, accumulation_deb_ref_std, refrozen_rain_deb_ref_std, netsnowmelt_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_ref_std)
# 12 year mass balance
massbalance_timeseries_12years('Catchment-wide average mass balance',2007, years, 0.055, 1.5, accumulation_krh_ref, refrozen_rain_krh_ref ,netsnowmelt_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref)
massbalance_timeseries_12years('Glacierized area mass balance',1980, years, 0.055, 1.5, accumulation_allgl_ref, refrozen_rain_allgl_ref ,netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref)
massbalance_timeseries_12years('Kaskawulsh mass balance',1980, years, 0.055, 1.5, accumulation_kw_ref, refrozen_rain_kw_ref ,netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref)
massbalance_timeseries_12years('Terminus mass balance',2007, years, 0.15, 10, accumulation_deb_ref, refrozen_rain_deb_ref ,netsnowmelt_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref)
# Average runoff
runoff_timeseries_average_SLRformat('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
runoff_timeseries_average_SLRformat('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref,gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std,All_glacierized_area)
runoff_timeseries_average_SLRformat('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_kw_ref, netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref,gl_icemelt_kw_ref_std, netsnowmelt_kw_ref_std, rain_runoff_kw_ref_std, superimp_icemelt_kw_ref_std,KRH_tributaries)
runoff_timeseries_average_SLRformat('m3','Runoff at the terminus: ',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref,gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std,np.ones((1,1)))
# 12 year runoff (eventually put individual pie charts on these plots)
runoff_timeseries_12years('m3','Catchment-wide average runoff',2007,years,450,3,gl_icemelt_krh_ref,netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,Sfc)
runoff_timeseries_12years('m3','Average runoff from the glacierized area',2007,years,450,3,gl_icemelt_allgl_ref,netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref,Sfc)
runoff_timeseries_12years('m3','Kaskawulsh runoff',2007,years,450,3,gl_icemelt_kw_ref,netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref,Sfc)
# 12 year pie charts
runoff_piecharts_12years('Catchment-wide average runoff',2007,years,gl_icemelt_krh_ref,netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref)
runoff_piecharts_12years('Average runoff from the glacierized area',2007,years,gl_icemelt_allgl_ref,netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref)
runoff_piecharts_12years('Kaskawulsh runoff',2007,years,gl_icemelt_kw_ref,netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref)
# annual runoff stacked bar
annualrunoff_stackedbar('Kaskawulsh River Headwaters: total annual runoff',years,netsnowmelt_krh_ref,gl_icemelt_krh_ref,superimp_icemelt_krh_ref,rain_runoff_krh_ref,KRH_tributaries)
annualrunoff_stackedbar('Glacierized area: total annual runoff',years,netsnowmelt_allgl_ref,gl_icemelt_allgl_ref,superimp_icemelt_allgl_ref,rain_runoff_allgl_ref,KRH_tributaries)
annualrunoff_stackedbar('Kaskawulsh: total annual runoff',years,netsnowmelt_kw_ref,gl_icemelt_kw_ref,superimp_icemelt_kw_ref,rain_runoff_kw_ref,KRH_tributaries)
# date of zero balance
date_of_zero_balance('Catchment',years,massbal_krh_ref)
date_of_zero_balance('All glacierized area',years,massbal_allgl_ref)
date_of_zero_balance('Kaskawulsh',years,massbal_kw_ref)
# annual glacier contribution to catchment runoff

# cumulative mass balance
cumulative_massbalance('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref)
cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref)
# annual mass balance bar chart
annual_massbalance_barchart('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref)
annual_massbalance_barchart('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref)

# distributed plots
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.888283378746594,stop=1,name='massbal')
distributed_average_mass_balance('REF_MODEL',np.arange(2007,2017+1),years,massbal_dist_ref,np.linspace(-16.3,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_runoff(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_glaciericemelt(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_snowrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,0.7,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_rainrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,0.3,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_SImelt(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,0.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)





# =============================================================================
# Rounce debris model plots
# =============================================================================
# Average mass balance
massbalance_timeseries_average('Catchment-wide average mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_krh_rounce, refrozen_rain_krh_rounce, netsnowmelt_krh_rounce, superimp_icemelt_krh_rounce, gl_icemelt_krh_rounce, accumulation_krh_rounce_std, refrozen_rain_krh_rounce_std, netsnowmelt_krh_rounce_std, superimp_icemelt_krh_rounce_std, gl_icemelt_krh_rounce_std)
massbalance_timeseries_average('Glacierized area mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_allgl_rounce, refrozen_rain_allgl_rounce, netsnowmelt_allgl_rounce, superimp_icemelt_allgl_rounce, gl_icemelt_allgl_rounce, accumulation_allgl_rounce_std, refrozen_rain_allgl_rounce_std, netsnowmelt_allgl_rounce_std, superimp_icemelt_allgl_rounce_std, gl_icemelt_allgl_rounce_std)
massbalance_timeseries_average('Kaskawulsh mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_kw_rounce, refrozen_rain_kw_rounce, netsnowmelt_kw_rounce, superimp_icemelt_kw_rounce, gl_icemelt_kw_rounce, accumulation_kw_rounce_std, refrozen_rain_kw_rounce_std, netsnowmelt_kw_rounce_std, superimp_icemelt_kw_rounce_std, gl_icemelt_kw_rounce_std)
massbalance_timeseries_average('Terminus mass balance: ',np.arange(2007,2017+1),years,0.3,20, accumulation_deb_rounce, refrozen_rain_deb_rounce, netsnowmelt_deb_rounce, superimp_icemelt_deb_rounce, gl_icemelt_deb_rounce, accumulation_deb_rounce_std, refrozen_rain_deb_rounce_std, netsnowmelt_deb_rounce_std, superimp_icemelt_deb_rounce_std, gl_icemelt_deb_rounce_std)
# 12 year mass balance
massbalance_timeseries_12years('Catchment-wide average mass balance',2007, years, 0.055, 1.5, accumulation_krh_rounce, refrozen_rain_krh_rounce ,netsnowmelt_krh_rounce, superimp_icemelt_krh_rounce, gl_icemelt_krh_rounce)
massbalance_timeseries_12years('Glacierized area mass balance',2007, years, 0.055, 1.5, accumulation_allgl_rounce, refrozen_rain_allgl_rounce ,netsnowmelt_allgl_rounce, superimp_icemelt_allgl_rounce, gl_icemelt_allgl_rounce)
massbalance_timeseries_12years('Kaskawulsh mass balance',2007, years, 0.055, 1.5, accumulation_kw_rounce, refrozen_rain_kw_rounce ,netsnowmelt_kw_rounce, superimp_icemelt_kw_rounce, gl_icemelt_kw_rounce)
massbalance_timeseries_12years('Terminus mass balance',2007, years, 0.15, 10, accumulation_deb_rounce, refrozen_rain_deb_rounce ,netsnowmelt_deb_rounce, superimp_icemelt_deb_rounce, gl_icemelt_deb_rounce)
# Average runoff
runoff_timeseries_average_SLRformat('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_krh_rounce, netsnowmelt_krh_rounce, rain_runoff_krh_rounce, superimp_icemelt_krh_rounce,gl_icemelt_krh_rounce_std, netsnowmelt_krh_rounce_std, rain_runoff_krh_rounce_std, superimp_icemelt_krh_rounce_std,Sfc)
runoff_timeseries_average_SLRformat('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_allgl_rounce, netsnowmelt_allgl_rounce, rain_runoff_allgl_rounce, superimp_icemelt_allgl_rounce,gl_icemelt_allgl_rounce_std, netsnowmelt_allgl_rounce_std, rain_runoff_allgl_rounce_std, superimp_icemelt_allgl_rounce_std,All_glacierized_area)
runoff_timeseries_average_SLRformat('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_kw_rounce, netsnowmelt_kw_rounce, rain_runoff_kw_rounce, superimp_icemelt_kw_rounce,gl_icemelt_kw_rounce_std, netsnowmelt_kw_rounce_std, rain_runoff_kw_rounce_std, superimp_icemelt_kw_rounce_std,KRH_tributaries)
runoff_timeseries_average_SLRformat('m3','Runoff at the terminus: ',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_rounce, netsnowmelt_deb_rounce, rain_runoff_deb_rounce, superimp_icemelt_deb_rounce,gl_icemelt_deb_rounce_std, netsnowmelt_deb_rounce_std, rain_runoff_deb_rounce_std, superimp_icemelt_deb_rounce_std,np.ones((1,1)))
# 12 year runoff (eventually put individual pie charts on these plots)
runoff_timeseries_12years('m3','Catchment-wide average runoff',2007,years,450,3,gl_icemelt_krh_rounce,netsnowmelt_krh_rounce, rain_runoff_krh_rounce, superimp_icemelt_krh_rounce,Sfc)
runoff_timeseries_12years('m3','Average runoff from the glacierized area',2007,years,450,3,gl_icemelt_allgl_rounce,netsnowmelt_allgl_rounce, rain_runoff_allgl_rounce, superimp_icemelt_allgl_rounce,Sfc)
runoff_timeseries_12years('m3','Kaskawulsh runoff',2007,years,450,3,gl_icemelt_kw_rounce,netsnowmelt_kw_rounce, rain_runoff_kw_rounce, superimp_icemelt_kw_rounce,Sfc)
# 12 year pie charts
runoff_piecharts_12years('Catchment-wide average runoff',2007,years,gl_icemelt_krh_rounce,netsnowmelt_krh_rounce, rain_runoff_krh_rounce, superimp_icemelt_krh_rounce)
runoff_piecharts_12years('Average runoff from the glacierized area',2007,years,gl_icemelt_allgl_rounce,netsnowmelt_allgl_rounce, rain_runoff_allgl_rounce, superimp_icemelt_allgl_rounce)
runoff_piecharts_12years('Kaskawulsh runoff',2007,years,gl_icemelt_kw_rounce,netsnowmelt_kw_rounce, rain_runoff_kw_rounce, superimp_icemelt_kw_rounce)
# annual runoff stacked bar
annualrunoff_stackedbar('Kaskawulsh River Headwaters: total annual runoff',years,netsnowmelt_krh_rounce,gl_icemelt_krh_rounce,superimp_icemelt_krh_rounce,rain_runoff_krh_rounce,KRH_tributaries)
annualrunoff_stackedbar('Glacierized area: total annual runoff',years,netsnowmelt_allgl_rounce,gl_icemelt_allgl_rounce,superimp_icemelt_allgl_rounce,rain_runoff_allgl_rounce,KRH_tributaries)
annualrunoff_stackedbar('Kaskawulsh: total annual runoff',years,netsnowmelt_kw_rounce,gl_icemelt_kw_rounce,superimp_icemelt_kw_rounce,rain_runoff_kw_rounce,KRH_tributaries)
# date of zero balance
date_of_zero_balance('Catchment',years,massbal_krh_rounce)
date_of_zero_balance('All glacierized area',years,massbal_allgl_rounce)
date_of_zero_balance('Kaskawulsh',years,massbal_kw_rounce)
# annual glacier contribution to catchment runoff

# cumulative mass balance
cumulative_massbalance('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_kw_rounce, refrozen_rain_kw_rounce, gl_icemelt_kw_rounce, netsnowmelt_kw_rounce, superimp_icemelt_kw_rounce)
cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_allgl_rounce, refrozen_rain_allgl_rounce, gl_icemelt_allgl_rounce, netsnowmelt_allgl_rounce, superimp_icemelt_allgl_rounce)
# annual mass balance bar chart
annual_massbalance_barchart('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_kw_rounce, refrozen_rain_kw_rounce, gl_icemelt_kw_rounce, netsnowmelt_kw_rounce, superimp_icemelt_kw_rounce)
annual_massbalance_barchart('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_allgl_rounce, refrozen_rain_allgl_rounce, gl_icemelt_allgl_rounce, netsnowmelt_allgl_rounce, superimp_icemelt_allgl_rounce)

# distributed plots
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.888283378746594,stop=1,name='massbal')
distributed_average_mass_balance('REF_MODEL',np.arange(2007,2017+1),years,massbal_dist_rounce,np.linspace(-16.3,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_runoff(np.arange(1980,2021+1),years,netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_glaciericemelt(np.arange(1980,2021+1),years,netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_snowrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce,np.linspace(0,0.7,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_rainrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce,np.linspace(0,0.3,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_SImelt(np.arange(1980,2021+1),years,netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce,np.linspace(0,0.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)








# =============================================================================
# REF_MODEL vs ROUNCE DEBRIS
# =============================================================================
# runoff differences
compare_hydrographs_differences('m3','Catchment-wide runoff: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_rounce, netsnowmelt_krh_rounce, rain_runoff_krh_rounce, superimp_icemelt_krh_rounce, gl_icemelt_krh_rounce_std, netsnowmelt_krh_rounce_std, rain_runoff_krh_rounce_std, superimp_icemelt_krh_rounce_std,Sfc)
compare_hydrographs_differences('m3','Runoff from the glacierized area: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_rounce, netsnowmelt_allgl_rounce, rain_runoff_allgl_rounce, superimp_icemelt_allgl_rounce, gl_icemelt_allgl_rounce_std, netsnowmelt_allgl_rounce_std, rain_runoff_allgl_rounce_std, superimp_icemelt_allgl_rounce_std,All_glacierized_area)
compare_hydrographs_differences('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_kw_ref, netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref_std, netsnowmelt_kw_ref_std, rain_runoff_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_rounce, netsnowmelt_kw_rounce, rain_runoff_kw_rounce, superimp_icemelt_kw_rounce, gl_icemelt_kw_rounce_std, netsnowmelt_kw_rounce_std, rain_runoff_kw_rounce_std, superimp_icemelt_kw_rounce_std,KRH_tributaries)
compare_hydrographs_differences('m3','Runoff from a debris-covered cell at the terminus: ',np.arange(1980,2021+1),years,-0.035,0.03, gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_rounce, netsnowmelt_deb_rounce, rain_runoff_deb_rounce, superimp_icemelt_deb_rounce, gl_icemelt_deb_rounce_std, netsnowmelt_deb_rounce_std, rain_runoff_deb_rounce_std, superimp_icemelt_deb_rounce_std,np.ones((1,1)))
# hydrograph comparison
compare_hydrographs('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_rounce, netsnowmelt_allgl_rounce, rain_runoff_allgl_rounce, superimp_icemelt_allgl_rounce, gl_icemelt_allgl_rounce_std, netsnowmelt_allgl_rounce_std, rain_runoff_allgl_rounce_std, superimp_icemelt_allgl_rounce_std,All_glacierized_area)
compare_hydrographs('m3','Runoff from a debris-covered cell at the terminus:\n',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_rounce, netsnowmelt_deb_rounce, rain_runoff_deb_rounce, superimp_icemelt_deb_rounce, gl_icemelt_deb_rounce_std, netsnowmelt_deb_rounce_std, rain_runoff_deb_rounce_std, superimp_icemelt_deb_rounce_std,np.ones((1,1)))
# mass balance difference:
massbalance_timeseries_difference('Catchment-wide mass balance: ',np.arange(2007,2017+1),years,0.001,0.03, accumulation_krh_ref, refrozen_rain_krh_ref, netsnowmelt_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref, accumulation_krh_ref_std, refrozen_rain_krh_ref_std, netsnowmelt_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_ref_std, accumulation_krh_rounce, refrozen_rain_krh_rounce, netsnowmelt_krh_rounce, superimp_icemelt_krh_rounce, gl_icemelt_krh_rounce, accumulation_krh_rounce_std, refrozen_rain_krh_rounce_std, netsnowmelt_krh_rounce_std, superimp_icemelt_krh_rounce_std, gl_icemelt_krh_rounce_std)
massbalance_timeseries_difference('Glacierized area mass balance: ',np.arange(2007,2017+1),years,0.001,0.03, accumulation_allgl_ref, refrozen_rain_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref, accumulation_allgl_ref_std, refrozen_rain_allgl_ref_std, netsnowmelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_ref_std, accumulation_allgl_rounce, refrozen_rain_allgl_rounce, netsnowmelt_allgl_rounce, superimp_icemelt_allgl_rounce, gl_icemelt_allgl_rounce, accumulation_allgl_rounce_std, refrozen_rain_allgl_rounce_std, netsnowmelt_allgl_rounce_std, superimp_icemelt_allgl_rounce_std, gl_icemelt_allgl_rounce_std)
massbalance_timeseries_difference('Kaskawulsh Glacier mass balance: ',np.arange(2007,2017+1),years,0.001,0.03, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std, accumulation_kw_rounce, refrozen_rain_kw_rounce, netsnowmelt_kw_rounce, superimp_icemelt_kw_rounce, gl_icemelt_kw_rounce, accumulation_kw_rounce_std, refrozen_rain_kw_rounce_std, netsnowmelt_kw_rounce_std, superimp_icemelt_kw_rounce_std, gl_icemelt_kw_rounce_std)
# mass balance comparison
massbalance_timeseries_comparison('Kaskawulsh mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std, accumulation_kw_rounce, refrozen_rain_kw_rounce, netsnowmelt_kw_rounce, superimp_icemelt_kw_rounce, gl_icemelt_kw_rounce, accumulation_kw_rounce_std, refrozen_rain_kw_rounce_std, netsnowmelt_kw_rounce_std, superimp_icemelt_kw_rounce_std, gl_icemelt_kw_rounce_std)
# date of zero balance 
compare_date_of_zero_balance(years,massbal_allgl_ref,massbal_allgl_rounce)
# Distributed differences
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.21052631578947367,stop=1,name='massbal_diff')
distributed_mass_balance_difference('ROUNCE_DEBRIS',np.arange(1980,2021+1),years,massbal_dist_ref,massbal_dist_rounce,np.linspace(-2,7.5,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.79,stop=1,name='runoff_diff')
distributed_runoff_difference('ROUNCE_DEBRIS',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce,np.linspace(-7.5,2,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_glaciermelt_difference('ROUNCE_DEBRIS',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce, np.linspace(-7.5,2,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.9649122807017544,stop=1,name='snowmelt_diff')
distributed_snowrunoff_difference('ROUNCE_DEBRIS',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce, np.linspace(-0.055,0.001,40),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=1,stop=0.7,name='rainrunoff_diff')
distributed_rainrunoff_difference('ROUNCE_DEBRIS',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce, np.linspace(-0.003,0,10),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.5,stop=1,name='SImelt_diff')
distributed_SImelt_difference('ROUNCE_DEBRIS',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_rounce,gl_icemelt_dist_rounce,superimp_icemelt_dist_rounce,rain_runoff_dist_rounce, np.linspace(-0.02,0.02,21),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)






















































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
    yrly_snowmelt.append(np.nansum(netsnowmelt_dist_ref[i]))
    yrly_icemelt.append(np.nansum(gl_icemelt_dist_ref[i]))
    yrly_SImelt.append(np.nansum(superimp_icemelt_dist_ref[i]))
    yrly_rain.append(np.nansum(rain_runoff_dist_ref[i]))

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


