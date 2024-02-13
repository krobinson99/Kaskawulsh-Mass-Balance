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
import cmocean
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable



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
distributed_snowrunoff_difference, distributed_SImelt_difference, surface_elevation_profiles, \
cumulative_and_annual_massbal, Bmod_vs_Bcal, distributed_runoff_components, \
runoff_contribution_timeseries, runoff_timeseries_with_zerophaseshift, distributed_runoff_allcomponents

REF_MODEL_PATH = 'D:/Model Runs/REF_MODEL/Sim_99999_v2'
ROUNCE_DEBRIS_PATH = 'D:/Model Runs/ROUNCE_DEBRIS'
UNCORRECTED_ACC_PATH = 'D:/Model Runs/UNCORRECTED_ACC'
NO_DEBRIS_PATH = 'D:/Model Runs/NO_DEBRIS'
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
Zgrid = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Zgrids/DEM_KRH_2022.txt')
Sfc = np.loadtxt(Sfc_grid)

Catchmentoutline = np.array(Sfc)
Catchmentoutline[np.where(np.isfinite(Sfc))] = 0
Catchmentoutline[np.where(np.isnan(Sfc))] = 1

KRH_tributaries = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Tributaries.txt')

KW_fluxgates = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Fluxgates.txt')

All_glacierized_area = np.array(Sfc)
All_glacierized_area[np.where(Sfc==1)] = np.nan

Otherice = np.array((Sfc))
Otherice[np.isfinite(KRH_tributaries)] = np.nan
Otherice[Otherice==1] = np.nan

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
#  REF_MODEL: Non KW Ice
# =============================================================================
massbal_nonKWice_ref, totalsnowmelt_nonKWice_ref, refreezing_nonKWice_ref, netsnowmelt_nonKWice_ref, gl_icemelt_nonKWice_ref, superimp_icemelt_nonKWice_ref, rain_nonKWice_ref, refrozen_rain_nonKWice_ref, rain_runoff_nonKWice_ref, accumulation_nonKWice_ref, \
snowdepth_nonKWice_ref, Ptau_nonKWice_ref, SI_nonKWice_ref = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'nonKWice',years)

massbal_nonKWice_ref_std, totalsnowmelt_nonKWice_ref_std, refreezing_nonKWice_ref_std, netsnowmelt_nonKWice_ref_std, gl_icemelt_nonKWice_ref_std, superimp_icemelt_nonKWice_ref_std, rain_nonKWice_ref_std, refrozen_rain_nonKWice_ref_std, rain_runoff_nonKWice_ref_std, accumulation_nonKWice_ref_std, \
snowdepth_nonKWice_ref_std, Ptau_nonKWice_ref_std, SI_nonKWice_ref_std = load_daily_timeseries(REF_MODEL_PATH,'REF_MODEL',varnames,'nonKWice_std',years)

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
#  NO_DEBRIS: Catchment-wide average
# =============================================================================
massbal_krh_nodeb, totalsnowmelt_krh_nodeb, refreezing_krh_nodeb, netsnowmelt_krh_nodeb, gl_icemelt_krh_nodeb, superimp_icemelt_krh_nodeb, rain_krh_nodeb, refrozen_rain_krh_nodeb, rain_runoff_krh_nodeb, accumulation_krh_nodeb, \
snowdepth_krh_nodeb, Ptau_krh_nodeb, SI_krh_nodeb = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'krh',years)

massbal_krh_nodeb_std, totalsnowmelt_krh_nodeb_std, refreezing_krh_nodeb_std, netsnowmelt_krh_nodeb_std, gl_icemelt_krh_nodeb_std, superimp_icemelt_krh_nodeb_std, rain_krh_nodeb_std, refrozen_rain_krh_nodeb_std, rain_runoff_krh_nodeb_std, accumulation_krh_nodeb_std, \
snowdepth_krh_nodeb_std, Ptau_krh_nodeb_std, SI_krh_nodeb_std = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'krh_std',years)

# =============================================================================
#  NO_DEBRIS: All glacierized area
# =============================================================================
massbal_allgl_nodeb, totalsnowmelt_allgl_nodeb, refreezing_allgl_nodeb, netsnowmelt_allgl_nodeb, gl_icemelt_allgl_nodeb, superimp_icemelt_allgl_nodeb, rain_allgl_nodeb, refrozen_rain_allgl_nodeb, rain_runoff_allgl_nodeb, accumulation_allgl_nodeb, \
snowdepth_allgl_nodeb, Ptau_allgl_nodeb, SI_allgl_nodeb = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'allgl',years)

massbal_allgl_nodeb_std, totalsnowmelt_allgl_nodeb_std, refreezing_allgl_nodeb_std, netsnowmelt_allgl_nodeb_std, gl_icemelt_allgl_nodeb_std, superimp_icemelt_allgl_nodeb_std, rain_allgl_nodeb_std, refrozen_rain_allgl_nodeb_std, rain_runoff_allgl_nodeb_std, accumulation_allgl_nodeb_std, \
snowdepth_allgl_nodeb_std, Ptau_allgl_nodeb_std, SI_allgl_nodeb_std = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'allgl_std',years)

# =============================================================================
#  NO_DEBRIS: Kaskawulsh
# =============================================================================
massbal_kw_nodeb, totalsnowmelt_kw_nodeb, refreezing_kw_nodeb, netsnowmelt_kw_nodeb, gl_icemelt_kw_nodeb, superimp_icemelt_kw_nodeb, rain_kw_nodeb, refrozen_rain_kw_nodeb, rain_runoff_kw_nodeb, accumulation_kw_nodeb, \
snowdepth_kw_nodeb, Ptau_kw_nodeb, SI_kw_nodeb = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'kw',years)

massbal_kw_nodeb_std, totalsnowmelt_kw_nodeb_std, refreezing_kw_nodeb_std, netsnowmelt_kw_nodeb_std, gl_icemelt_kw_nodeb_std, superimp_icemelt_kw_nodeb_std, rain_kw_nodeb_std, refrozen_rain_kw_nodeb_std, rain_runoff_kw_nodeb_std, accumulation_kw_nodeb_std, \
snowdepth_kw_nodeb_std, Ptau_kw_nodeb_std, SI_kw_nodeb_std = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'kw_std',years)

# =============================================================================
#  NO_DEBRIS: Debris-covered cell at the terminus: NOT GENERATED YET
# =============================================================================
massbal_deb_nodeb, totalsnowmelt_deb_nodeb, refreezing_deb_nodeb, netsnowmelt_deb_nodeb, gl_icemelt_deb_nodeb, superimp_icemelt_deb_nodeb, rain_deb_nodeb, refrozen_rain_deb_nodeb, rain_runoff_deb_nodeb, accumulation_deb_nodeb, \
snowdepth_deb_nodeb, Ptau_deb_nodeb, SI_deb_nodeb = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'deb',years)

massbal_deb_nodeb_std, totalsnowmelt_deb_nodeb_std, refreezing_deb_nodeb_std, netsnowmelt_deb_nodeb_std, gl_icemelt_deb_nodeb_std, superimp_icemelt_deb_nodeb_std, rain_deb_nodeb_std, refrozen_rain_deb_nodeb_std, rain_runoff_deb_nodeb_std, accumulation_deb_nodeb_std, \
snowdepth_deb_nodeb_std, Ptau_deb_nodeb_std, SI_deb_nodeb_std = load_daily_timeseries(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,'deb_std',years)

# =============================================================================
#  NO_DEBRIS: Distributed totals
# =============================================================================
massbal_dist_nodeb, totalsnowmelt_dist_nodeb, refreezing_dist_nodeb, netsnowmelt_dist_nodeb, gl_icemelt_dist_nodeb, superimp_icemelt_dist_nodeb, rain_dist_nodeb, refrozen_rain_dist_nodeb, rain_runoff_dist_nodeb, accumulation_dist_nodeb, \
snowdepth_dist_nodeb, Ptau_dist_nodeb, SI_dist_nodeb = load_distributed_fields(NO_DEBRIS_PATH,'NO_DEBRIS',varnames,years)












# =============================================================================
#  UNCORRECTED_ACC: Catchment-wide average
# =============================================================================
massbal_krh_uncorracc, totalsnowmelt_krh_uncorracc, refreezing_krh_uncorracc, netsnowmelt_krh_uncorracc, gl_icemelt_krh_uncorracc, superimp_icemelt_krh_uncorracc, rain_krh_uncorracc, refrozen_rain_krh_uncorracc, rain_runoff_krh_uncorracc, accumulation_krh_uncorracc, \
snowdepth_krh_uncorracc, Ptau_krh_uncorracc, SI_krh_uncorracc = load_daily_timeseries(UNCORRECTED_ACC_PATH,'UNCORRECTED_ACC',varnames,'krh',years)

massbal_krh_uncorracc_std, totalsnowmelt_krh_uncorracc_std, refreezing_krh_uncorracc_std, netsnowmelt_krh_uncorracc_std, gl_icemelt_krh_uncorracc_std, superimp_icemelt_krh_uncorracc_std, rain_krh_uncorracc_std, refrozen_rain_krh_uncorracc_std, rain_runoff_krh_uncorracc_std, accumulation_krh_uncorracc_std, \
snowdepth_krh_uncorracc_std, Ptau_krh_uncorracc_std, SI_krh_uncorracc_std = load_daily_timeseries(UNCORRECTED_ACC_PATH,'UNCORRECTED_ACC',varnames,'krh_std',years)

# =============================================================================
#  UNCORRECTED_ACC: All glacierized area
# =============================================================================
massbal_allgl_uncorracc, totalsnowmelt_allgl_uncorracc, refreezing_allgl_uncorracc, netsnowmelt_allgl_uncorracc, gl_icemelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc, rain_allgl_uncorracc, refrozen_rain_allgl_uncorracc, rain_runoff_allgl_uncorracc, accumulation_allgl_uncorracc, \
snowdepth_allgl_uncorracc, Ptau_allgl_uncorracc, SI_allgl_uncorracc = load_daily_timeseries(UNCORRECTED_ACC_PATH,'UNCORRECTED_ACC',varnames,'allgl',years)

massbal_allgl_uncorracc_std, totalsnowmelt_allgl_uncorracc_std, refreezing_allgl_uncorracc_std, netsnowmelt_allgl_uncorracc_std, gl_icemelt_allgl_uncorracc_std, superimp_icemelt_allgl_uncorracc_std, rain_allgl_uncorracc_std, refrozen_rain_allgl_uncorracc_std, rain_runoff_allgl_uncorracc_std, accumulation_allgl_uncorracc_std, \
snowdepth_allgl_uncorracc_std, Ptau_allgl_uncorracc_std, SI_allgl_uncorracc_std = load_daily_timeseries(UNCORRECTED_ACC_PATH,'UNCORRECTED_ACC',varnames,'allgl_std',years)

# =============================================================================
#  UNCORRECTED_ACC: Kaskawulsh
# =============================================================================
massbal_kw_uncorracc, totalsnowmelt_kw_uncorracc, refreezing_kw_uncorracc, netsnowmelt_kw_uncorracc, gl_icemelt_kw_uncorracc, superimp_icemelt_kw_uncorracc, rain_kw_uncorracc, refrozen_rain_kw_uncorracc, rain_runoff_kw_uncorracc, accumulation_kw_uncorracc, \
snowdepth_kw_uncorracc, Ptau_kw_uncorracc, SI_kw_uncorracc = load_daily_timeseries(UNCORRECTED_ACC_PATH,'UNCORRECTED_ACC',varnames,'kw',years)

massbal_kw_uncorracc_std, totalsnowmelt_kw_uncorracc_std, refreezing_kw_uncorracc_std, netsnowmelt_kw_uncorracc_std, gl_icemelt_kw_uncorracc_std, superimp_icemelt_kw_uncorracc_std, rain_kw_uncorracc_std, refrozen_rain_kw_uncorracc_std, rain_runoff_kw_uncorracc_std, accumulation_kw_uncorracc_std, \
snowdepth_kw_uncorracc_std, Ptau_kw_uncorracc_std, SI_kw_uncorracc_std = load_daily_timeseries(UNCORRECTED_ACC_PATH,'UNCORRECTED_ACC',varnames,'kw_std',years)

# =============================================================================
#  UNCORRECTED_ACC: Distributed totals
# =============================================================================
massbal_dist_uncorracc, totalsnowmelt_dist_uncorracc, refreezing_dist_uncorracc, netsnowmelt_dist_uncorracc, gl_icemelt_dist_uncorracc, superimp_icemelt_dist_uncorracc, rain_dist_uncorracc, refrozen_rain_dist_uncorracc, rain_runoff_dist_uncorracc, accumulation_dist_uncorracc, \
snowdepth_dist_uncorracc, Ptau_dist_uncorracc, SI_dist_uncorracc = load_distributed_fields(UNCORRECTED_ACC_PATH,'UNCORRECTED_ACC',varnames,years)







# =============================================================================
# Reference model plots
# =============================================================================
# Average mass balance
massbalance_timeseries_average('Catchment-wide average mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_krh_ref, refrozen_rain_krh_ref, netsnowmelt_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref, accumulation_krh_ref_std, refrozen_rain_krh_ref_std, netsnowmelt_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_ref_std)
massbalance_timeseries_average('Glacierized area mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_allgl_ref, refrozen_rain_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref, accumulation_allgl_ref_std, refrozen_rain_allgl_ref_std, netsnowmelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_ref_std,All_glacierized_area)
massbalance_timeseries_average('Kaskawulsh mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std,KRH_tributaries)
massbalance_timeseries_average('Terminus mass balance: ',np.arange(2007,2017+1),years,0.15,10, accumulation_deb_ref, refrozen_rain_deb_ref, netsnowmelt_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref, accumulation_deb_ref_std, refrozen_rain_deb_ref_std, netsnowmelt_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_ref_std)
# 12 year mass balance
massbalance_timeseries_12years('Catchment-wide average mass balance',2007, years, 0.055, 1.5, accumulation_krh_ref, refrozen_rain_krh_ref ,netsnowmelt_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref)
massbalance_timeseries_12years('Glacierized area mass balance',1980, years, 0.055, 1.5, accumulation_allgl_ref, refrozen_rain_allgl_ref ,netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref)
massbalance_timeseries_12years('Kaskawulsh mass balance',1980, years, 0.055, 1.5, accumulation_kw_ref, refrozen_rain_kw_ref ,netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref)
massbalance_timeseries_12years('Terminus mass balance',2007, years, 0.15, 10, accumulation_deb_ref, refrozen_rain_deb_ref ,netsnowmelt_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref)

# Average runoff
runoff_timeseries_average_SLRformat('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,450,2.7,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc) #plt.savefig('D:/Model Runs/REF_MODEL/Plots/Hydrograph_Catchmentwide_1980-2022_REFMODEL.pdf',bbox_inches='tight')
runoff_timeseries_average_SLRformat('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref,gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std,All_glacierized_area)
runoff_timeseries_average_SLRformat('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_kw_ref, netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref,gl_icemelt_kw_ref_std, netsnowmelt_kw_ref_std, rain_runoff_kw_ref_std, superimp_icemelt_kw_ref_std,KRH_tributaries)
runoff_timeseries_average_SLRformat('m3','Runoff at the terminus: ',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref,gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std,np.ones((1,1)))

runoff_percentile_curves('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)


t, tr1980s, gl1980s, sn1980s, ra1980s, si1980s = runoff_timeseries_with_zerophaseshift(ax1,51,3,'m3','Catchment-wide average runoff: ',np.arange(1980,1990),years,370,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
t, tr1990s, gl1990s, sn1990s, ra1990s, si1990s = runoff_timeseries_with_zerophaseshift(ax1,51,3,'m3','Catchment-wide average runoff: ',np.arange(1990,2000),years,370,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
t, tr2000s, gl2000s, sn2000s, ra2000s, si2000s = runoff_timeseries_with_zerophaseshift(ax1,51,3,'m3','Catchment-wide average runoff: ',np.arange(2000,2010),years,370,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
t, tr2010s, gl2010s, sn2010s, ra2010s, si2010s = runoff_timeseries_with_zerophaseshift(ax1,51,3,'m3','Catchment-wide average runoff: ',np.arange(2010,2020),years,370,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)


fig = plt.figure(figsize=(10,6),constrained_layout=True)
gs = GridSpec(4, 2, figure=fig)
ax1 = fig.add_subplot(gs[0:2, 0])
ax1.plot(t[182:],tr1980s[182:],label='1980-1990',color='khaki',linewidth=2)
ax1.plot(t[182:],tr1990s[182:],label='1990-2000',color='mediumaquamarine',linewidth=2)
ax1.plot(t[182:],tr2000s[182:],label='2000-2010',color='steelblue',linewidth=2)
ax1.plot(t[182:],tr2010s[182:],label='2010-2020',color='midnightblue',linewidth=2)
ax1.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
ax1.set_xticklabels(labels=['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],fontsize=14)
ax1.grid()
ax1.tick_params(axis='both',labelsize=14)
ax1.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=14)
ax1.set_xlim(182,365)
ax1.set_ylim(0,325)
ax1.text(185,304,'a) Total runoff',fontsize=14,weight='bold')
ax1.text(185,251,'Max Q (m$^3$ s$^{-1}$)',fontsize=14)
ax1.text(185,225,str(int(np.round(np.max(tr1980s)))),color='khaki',fontsize=14,weight='bold')
ax1.text(185,200,str(int(np.round(np.max(tr1990s)))),color='mediumaquamarine',fontsize=14,weight='bold')
ax1.text(185,175,str(int(np.round(np.max(tr2000s)))),color='steelblue',fontsize=14,weight='bold')
ax1.text(185,150,str(int(np.round(np.max(tr2010s)))),color='midnightblue',fontsize=14,weight='bold')
#ax1.vlines(t[np.where(tr1980s==np.max(tr1980s))],0,np.max(tr1980s),color='khaki',linewidth=3,linestyle='--')
#ax1.vlines(t[np.where(tr1990s==np.max(tr1990s))],0,np.max(tr1990s),color='mediumaquamarine',linewidth=3,linestyle='--')
#ax1.vlines(t[np.where(tr2000s==np.max(tr2000s))],0,np.max(tr2000s),color='steelblue',linewidth=3,linestyle='--')
#ax1.vlines(t[np.where(tr2010s==np.max(tr2010s))],0,np.max(tr2010s),color='midnightblue',linewidth=3,linestyle='--')

# Print max discharge rate and date on 

ax2 = fig.add_subplot(gs[2:, 0])
ax2.plot(t[182:],gl1980s[182:],label='1980-1990',color='khaki',linewidth=2)
ax2.plot(t[182:],gl1990s[182:],label='1990-2000',color='mediumaquamarine',linewidth=2)
ax2.plot(t[182:],gl2000s[182:],label='2000-2010',color='steelblue',linewidth=2)
ax2.plot(t[182:],gl2010s[182:],label='2010-2020',color='midnightblue',linewidth=2)
ax2.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
ax2.set_xticklabels(labels=['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],fontsize=14)
ax2.grid()
ax2.tick_params(axis='both',labelsize=14)
ax2.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=14)
ax2.set_xlim(182,365)
ax2.set_ylim(0,225)
ax2.text(185,203,'b) Glacier ice runoff',fontsize=14,weight='bold')
ax2.text(185,175,'Max Q (m$^3$ s$^{-1}$)',fontsize=14)
ax2.text(185,155,str(int(np.round(np.max(gl1980s)))),color='khaki',fontsize=14,weight='bold')
ax2.text(185,135,str(int(np.round(np.max(gl1990s)))),color='mediumaquamarine',fontsize=14,weight='bold')
ax2.text(185,115,str(int(np.round(np.max(gl2000s)))),color='steelblue',fontsize=14,weight='bold')
ax2.text(185,95,str(int(np.round(np.max(gl2010s)))),color='midnightblue',fontsize=14,weight='bold')
#ax2.vlines(t[np.where(gl1980s==np.max(gl1980s))],0,np.max(gl1980s),color='khaki',linewidth=3,linestyle='--')
#ax2.vlines(t[np.where(gl1990s==np.max(gl1990s))],0,np.max(gl1990s),color='mediumaquamarine',linewidth=3,linestyle='--')
#ax2.vlines(t[np.where(gl2000s==np.max(gl2000s))],0,np.max(gl2000s),color='steelblue',linewidth=3,linestyle='--')
#ax2.vlines(t[np.where(gl2010s==np.max(gl2010s))],0,np.max(gl2010s),color='midnightblue',linewidth=3,linestyle='--')

ax3 = fig.add_subplot(gs[0:2,1])
ax3.plot(t[182:],sn1980s[182:],label='1980-1990',color='khaki',linewidth=2)
ax3.plot(t[182:],sn1990s[182:],label='1990-2000',color='mediumaquamarine',linewidth=2)
ax3.plot(t[182:],sn2000s[182:],label='2000-2010',color='steelblue',linewidth=2)
ax3.plot(t[182:],sn2010s[182:],label='2010-2020',color='midnightblue',linewidth=2)
ax3.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
ax3.set_xticklabels(labels=['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],fontsize=14)
ax3.grid()
ax3.tick_params(axis='both',labelsize=14)
ax3.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=14)
ax3.set_xlim(182,365)
ax3.set_ylim(0,100)
ax3.text(185,92,'c) Snow runoff',fontsize=14,weight='bold')
ax3.text(185,80,'Max Q (m$^3$ s$^{-1}$)',fontsize=14)
ax3.text(185,72,str(int(np.round(np.max(sn1980s)))),color='khaki',fontsize=14,weight='bold')
ax3.text(185,64,str(int(np.round(np.max(sn1990s)))),color='mediumaquamarine',fontsize=14,weight='bold')
ax3.text(185,56,str(int(np.round(np.max(sn2000s)))),color='steelblue',fontsize=14,weight='bold')
ax3.text(185,48,str(int(np.round(np.max(sn2010s)))),color='midnightblue',fontsize=14,weight='bold')
#ax3.vlines(t[np.where(sn1980s==np.max(sn1980s))],0,np.max(sn1980s),color='khaki',linewidth=3,linestyle='--')
#ax3.vlines(t[np.where(sn1990s==np.max(sn1990s))],0,np.max(sn1990s),color='mediumaquamarine',linewidth=3,linestyle='--')
#ax3.vlines(t[np.where(sn2000s==np.max(sn2000s))],0,np.max(sn2000s),color='steelblue',linewidth=3,linestyle='--')
#ax3.vlines(t[np.where(sn2010s==np.max(sn2010s))],0,np.max(sn2010s),color='midnightblue',linewidth=3,linestyle='--')


ax4 = fig.add_subplot(gs[2, 1])
ax4.plot(t[182:],ra1980s[182:],label='1980-1990',color='khaki',linewidth=2)
ax4.plot(t[182:],ra1990s[182:],label='1990-2000',color='mediumaquamarine',linewidth=2)
ax4.plot(t[182:],ra2000s[182:],label='2000-2010',color='steelblue',linewidth=2)
ax4.plot(t[182:],ra2010s[182:],label='2010-2020',color='midnightblue',linewidth=2)
ax4.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
ax4.set_xticklabels(labels=['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],fontsize=14)
ax4.grid()
ax4.tick_params(axis='both',labelsize=14)
#ax4.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=12)
ax4.set_xlim(182,365)
ax4.set_ylim(0,20)
ax4.set_yticks(ticks=[0,5,10,15,20])
ax4.text(185,16,'d) Rain',fontsize=14,weight='bold')
ax4.text(185,11,'Max Q (m$^3$ s$^{-1}$)',fontsize=14)
ax4.text(185,8,str(int(np.round(np.max(ra1980s)))),color='khaki',fontsize=14,weight='bold')
ax4.text(185,5,str(int(np.round(np.max(ra1990s)))),color='mediumaquamarine',fontsize=14,weight='bold')
ax4.text(215,8,str(int(np.round(np.max(ra2000s)))),color='steelblue',fontsize=14,weight='bold')
ax4.text(215,5,str(int(np.round(np.max(ra2010s)))),color='midnightblue',fontsize=14,weight='bold')
#ax4.vlines(t[np.where(ra1980s==np.max(ra1980s))],0,np.max(ra1980s),color='khaki',linewidth=3,linestyle='--')
#ax4.vlines(t[np.where(ra1990s==np.max(ra1990s))],0,np.max(ra1990s),color='mediumaquamarine',linewidth=3,linestyle='--')
#ax4.vlines(t[np.where(ra2000s==np.max(ra2000s))],0,np.max(ra2000s),color='steelblue',linewidth=3,linestyle='--')
#ax4.vlines(t[np.where(ra2010s==np.max(ra2010s))],0,np.max(ra2010s),color='midnightblue',linewidth=3,linestyle='--')

ax5 = fig.add_subplot(gs[3, 1])
ax5.plot(t[182:],si1980s[182:],label='1980-1990',color='khaki',linewidth=2)
ax5.plot(t[182:],si1990s[182:],label='1990-2000',color='mediumaquamarine',linewidth=2)
ax5.plot(t[182:],si2000s[182:],label='2000-2010',color='steelblue',linewidth=2)
ax5.plot(t[182:],si2010s[182:],label='2010-2020',color='midnightblue',linewidth=2)
ax5.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
ax5.set_xticklabels(labels=['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],fontsize=14)
ax5.grid()
ax5.tick_params(axis='both',labelsize=14)
#ax5.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=12)
ax5.set_xlim(182,365)
ax5.set_ylim(0,14)
ax5.set_yticks(ticks=np.arange(0,15,4))
ax5.text(185,11,'e) Refrozen ice runoff',fontsize=14,weight='bold')
ax5.text(185,8,'Max Q (m$^3$ s$^{-1}$)',fontsize=14)
ax5.text(185,5,str(int(np.round(np.max(si1980s)))),color='khaki',fontsize=14,weight='bold')
ax5.text(185,2,str(int(np.round(np.max(si1990s)))),color='mediumaquamarine',fontsize=14,weight='bold')
ax5.text(205,5,str(int(np.round(np.max(si2000s)))),color='steelblue',fontsize=14,weight='bold')
ax5.text(205,2,str(int(np.round(np.max(si2010s)))),color='midnightblue',fontsize=14,weight='bold')
#ax5.vlines(t[np.where(si1980s==np.max(si1980s))],0,np.max(si1980s),color='khaki',linewidth=3,linestyle='--')
#ax5.vlines(t[np.where(si1990s==np.max(si1990s))],0,np.max(si1990s),color='mediumaquamarine',linewidth=3,linestyle='--')
#ax5.vlines(t[np.where(si2000s==np.max(si2000s))],0,np.max(si2000s),color='steelblue',linewidth=3,linestyle='--')
#ax5.vlines(t[np.where(si2010s==np.max(si2010s))],0,np.max(si2010s),color='midnightblue',linewidth=3,linestyle='--')

fig.text(0.51,0.26,'Discharge (m$^3$ s$^{-1}$)',fontsize=14,rotation=90,verticalalignment='center')

handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(0.999,1.07),fontsize=14, ncol=4, borderaxespad=0.19)
#fig.savefig('D:/Model Runs/REF_MODEL/Plots/Refmodel_smoothed_runoff_components.pdf',bbox_inches='tight')




# 12 year runoff (eventually put individual pie charts on these plots)
runoff_timeseries_12years('m3','Catchment-wide average runoff',2007,years,450,3,gl_icemelt_krh_ref,netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,Sfc)
runoff_timeseries_12years('m3','Average runoff from the glacierized area',2007,years,450,3,gl_icemelt_allgl_ref,netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref,Sfc)
runoff_timeseries_12years('m3','Kaskawulsh runoff',2007,years,450,3,gl_icemelt_kw_ref,netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref,Sfc)

# 12 year pie charts
runoff_piecharts_12years('Catchment-wide average runoff',2007,years,gl_icemelt_krh_ref,netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref)
runoff_piecharts_12years('Average runoff from the glacierized area',2007,years,gl_icemelt_allgl_ref,netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref)
runoff_piecharts_12years('Kaskawulsh runoff',2007,years,gl_icemelt_kw_ref,netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref)

# runoff contribution timeseries
runoff_contribution_timeseries(years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,11,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
# annual runoff stacked bar
annualrunoff_stackedbar('Kaskawulsh River Headwaters: total annual runoff',np.arange(1980,2022),years,netsnowmelt_krh_ref,gl_icemelt_krh_ref,superimp_icemelt_krh_ref,rain_runoff_krh_ref,Sfc) #plt.savefig('D:\Model Runs\REF_MODEL\Plots\KRH_annualrunoff_stacked.pdf',bbox_inches='tight')
annualrunoff_stackedbar('Glacierized area: total annual runoff',np.arange(1980,2022),years,netsnowmelt_allgl_ref,gl_icemelt_allgl_ref,superimp_icemelt_allgl_ref,rain_runoff_allgl_ref,All_glacierized_area)
annualrunoff_stackedbar('Kaskawulsh: total annual runoff',np.arange(1980,2022),years,netsnowmelt_kw_ref,gl_icemelt_kw_ref,superimp_icemelt_kw_ref,rain_runoff_kw_ref,KRH_tributaries)

# date of zero balance
date_of_zero_balance('Catchment',years,massbal_krh_ref)
date_of_zero_balance(years,massbal_allgl_ref) #plt.savefig('D:\Model Runs\REF_MODEL\Plots\Refmodel_B0_date.pdf',bbox_inches='tight')
date_of_zero_balance('Kaskawulsh',years,massbal_kw_ref)

date_of_zero_balance_with_mbcurve(years,massbal_allgl_ref,np.linspace(-10,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries,years,np.arange(1980,2022),0.027,1, accumulation_allgl_ref, refrozen_rain_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref, accumulation_allgl_ref_std, refrozen_rain_allgl_ref_std, netsnowmelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_ref_std,All_glacierized_area) #fig.savefig('D:/Model Runs/REF_MODEL/Plots/Refmodel_B0date_MBcurve.pdf',bbox_inches='tight')

# Cumulative mass balance
cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area)
Calculate_mass_loss_rates() # Calculate the cumulative mass loss from each decade
Calculate_KRH_contribution_to_Alaska_massloss()
cumulative_and_annual_massbal('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,22,2.5,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area) #plt.savefig('D:\Model Runs\REF_MODEL\Plots\Refmodel_MB_timeseries_1979-2022.pdf',bbox_inches='tight')

mb_vs_runoff(years,netsnowmelt_allgl_ref,gl_icemelt_allgl_ref,superimp_icemelt_allgl_ref,rain_runoff_allgl_ref,All_glacierized_area,accumulation_allgl_ref,refrozen_rain_allgl_ref)

# distributed plots
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.8298755186721992,stop=1,name='massbal')
distributed_average_mass_balance('REF_MODEL',np.arange(1980,2022),years,massbal_dist_ref,np.linspace(-10,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries) #plt.savefig('D:\Model Runs\REF_MODEL\Plots\Refmodel_MB_1979-2022.pdf',bbox_inches='tight')

distributed_runoff(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,11,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_glaciericemelt(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,10,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_snowrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,0.7,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_rainrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,0.3,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_SImelt(np.arange(1980,2021+1),years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,0.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_runoff_allcomponents(years,Sfc,Xgrid,Ygrid,Zgrid,Catchmentoutline,gl_icemelt_dist_ref,netsnowmelt_dist_ref,rain_runoff_dist_ref,superimp_icemelt_dist_ref)

distributed_snowpack_decadal(years,Sfc,Xgrid,Ygrid,Zgrid,Catchmentoutline,snowdepth_dist_ref)


AAR_over_time(years,years,massbal_dist_ref,massbal_dist_uncorracc)
Bmod_vs_Bcal(np.arange(2007,2018),KRH_tributaries,KW_fluxgates,massbal_dist_ref)

glacier_contribution_to_annualrunoff(years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,11,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
t, g, s, r, si = runoff_contribution_timeseries2(years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,np.linspace(0,11,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

# =============================================================================
# Plot refreezing vs Ptau over time 
# =============================================================================
Rmelt = []
Rrain = []
Ptau = []
for i in range(0,43):
    Ptau.append(np.mean(Ptau_dist_ref[i][np.isfinite(KRH_tributaries)]))
    Rmelt.append(np.mean(refreezing_dist_ref[i][np.isfinite(KRH_tributaries)]))
    Rrain.append(np.mean(refrozen_rain_dist_ref[i][np.isfinite(KRH_tributaries)]))

Rmelt = []
Rrain = []
Ptau = []
for i in range(0,43):
    Ptau.append(np.mean(Ptau_dist_ref[i][KRH_tributaries == 1]))
    Rmelt.append(np.mean(refreezing_dist_ref[i][KRH_tributaries == 1]))
    Rrain.append(np.mean(refrozen_rain_dist_ref[i][KRH_tributaries == 1]))

plt.plot(np.array(Rmelt)+np.array(Rrain));plt.plot(Ptau)

Ptau_area = []
Dry_snow_zone = []
for i in range(0,43):
    #Ptau_area.append(np.where(Ptau_dist_ref[i] >np.nanmax(totalsnowmelt_dist_ref) )[0].shape[0]*(0.2*0.2))
    Ptau_area.append(np.where(Ptau_dist_ref[i] >1)[0].shape[0]*(0.2*0.2))
    #Dry_snow_zone.append(np.where(refreezing_dist_ref[i]==0)[0].shape[0]*(0.2*0.2))
    Dry_snow_zone.append(np.where(totalsnowmelt_dist_ref[i]==0)[0].shape[0]*(0.2*0.2))



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
distributed_average_mass_balance('ROUNCE_DEBRIS',np.arange(2007,2017+1),years,massbal_dist_rounce,np.linspace(-16.3,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
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
massbalance_timeseries_difference('Kaskawulsh Glacier mass balance: ','ROUNCE_DEBRIS',np.arange(2007,2017+1),years,0.008,0.4, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std, accumulation_kw_rounce, refrozen_rain_kw_rounce, netsnowmelt_kw_rounce, superimp_icemelt_kw_rounce, gl_icemelt_kw_rounce, accumulation_kw_rounce_std, refrozen_rain_kw_rounce_std, netsnowmelt_kw_rounce_std, superimp_icemelt_kw_rounce_std, gl_icemelt_kw_rounce_std)
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

surface_elevation_profiles()


# =============================================================================
# NO_DEBRIS model plots
# =============================================================================
# Average mass balance
massbalance_timeseries_average('Catchment-wide average mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_krh_nodeb, refrozen_rain_krh_nodeb, netsnowmelt_krh_nodeb, superimp_icemelt_krh_nodeb, gl_icemelt_krh_nodeb, accumulation_krh_nodeb_std, refrozen_rain_krh_nodeb_std, netsnowmelt_krh_nodeb_std, superimp_icemelt_krh_nodeb_std, gl_icemelt_krh_nodeb_std,Sfc)
massbalance_timeseries_average('Glacierized area mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_allgl_nodeb, refrozen_rain_allgl_nodeb, netsnowmelt_allgl_nodeb, superimp_icemelt_allgl_nodeb, gl_icemelt_allgl_nodeb, accumulation_allgl_nodeb_std, refrozen_rain_allgl_nodeb_std, netsnowmelt_allgl_nodeb_std, superimp_icemelt_allgl_nodeb_std, gl_icemelt_allgl_nodeb_std,All_glacierized_area)
massbalance_timeseries_average('Kaskawulsh mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_kw_nodeb, refrozen_rain_kw_nodeb, netsnowmelt_kw_nodeb, superimp_icemelt_kw_nodeb, gl_icemelt_kw_nodeb, accumulation_kw_nodeb_std, refrozen_rain_kw_nodeb_std, netsnowmelt_kw_nodeb_std, superimp_icemelt_kw_nodeb_std, gl_icemelt_kw_nodeb_std,KRH_tributaries)
massbalance_timeseries_average('Terminus mass balance: ',np.arange(2007,2017+1),years,0.3,20, accumulation_deb_nodeb, refrozen_rain_deb_nodeb, netsnowmelt_deb_nodeb, superimp_icemelt_deb_nodeb, gl_icemelt_deb_nodeb, accumulation_deb_nodeb_std, refrozen_rain_deb_nodeb_std, netsnowmelt_deb_nodeb_std, superimp_icemelt_deb_nodeb_std, gl_icemelt_deb_nodeb_std,np.ones((1,1)))
# 12 year mass balance
massbalance_timeseries_12years('Catchment-wide average mass balance',2007, years, 0.055, 1.5, accumulation_krh_nodeb, refrozen_rain_krh_nodeb ,netsnowmelt_krh_nodeb, superimp_icemelt_krh_nodeb, gl_icemelt_krh_nodeb)
massbalance_timeseries_12years('Glacierized area mass balance',2007, years, 0.055, 1.5, accumulation_allgl_nodeb, refrozen_rain_allgl_nodeb ,netsnowmelt_allgl_nodeb, superimp_icemelt_allgl_nodeb, gl_icemelt_allgl_nodeb)
massbalance_timeseries_12years('Kaskawulsh mass balance',2007, years, 0.055, 1.5, accumulation_kw_nodeb, refrozen_rain_kw_nodeb ,netsnowmelt_kw_nodeb, superimp_icemelt_kw_nodeb, gl_icemelt_kw_nodeb)
massbalance_timeseries_12years('Terminus mass balance',2007, years, 0.15, 10, accumulation_deb_nodeb, refrozen_rain_deb_nodeb ,netsnowmelt_deb_nodeb, superimp_icemelt_deb_nodeb, gl_icemelt_deb_nodeb)
# Average runoff
runoff_timeseries_average_SLRformat('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_krh_nodeb, netsnowmelt_krh_nodeb, rain_runoff_krh_nodeb, superimp_icemelt_krh_nodeb,gl_icemelt_krh_nodeb_std, netsnowmelt_krh_nodeb_std, rain_runoff_krh_nodeb_std, superimp_icemelt_krh_nodeb_std,Sfc)
runoff_timeseries_average_SLRformat('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_allgl_nodeb, netsnowmelt_allgl_nodeb, rain_runoff_allgl_nodeb, superimp_icemelt_allgl_nodeb,gl_icemelt_allgl_nodeb_std, netsnowmelt_allgl_nodeb_std, rain_runoff_allgl_nodeb_std, superimp_icemelt_allgl_nodeb_std,All_glacierized_area)
runoff_timeseries_average_SLRformat('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_kw_nodeb, netsnowmelt_kw_nodeb, rain_runoff_kw_nodeb, superimp_icemelt_kw_nodeb,gl_icemelt_kw_nodeb_std, netsnowmelt_kw_nodeb_std, rain_runoff_kw_nodeb_std, superimp_icemelt_kw_nodeb_std,KRH_tributaries)
runoff_timeseries_average_SLRformat('m3','Runoff at the terminus: ',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_nodeb, netsnowmelt_deb_nodeb, rain_runoff_deb_nodeb, superimp_icemelt_deb_nodeb,gl_icemelt_deb_nodeb_std, netsnowmelt_deb_nodeb_std, rain_runoff_deb_nodeb_std, superimp_icemelt_deb_nodeb_std,np.ones((1,1)))
# 12 year runoff (eventually put individual pie charts on these plots)
runoff_timeseries_12years('m3','Catchment-wide average runoff',2007,years,450,3,gl_icemelt_krh_nodeb,netsnowmelt_krh_nodeb, rain_runoff_krh_nodeb, superimp_icemelt_krh_nodeb,Sfc)
runoff_timeseries_12years('m3','Average runoff from the glacierized area',2007,years,450,3,gl_icemelt_allgl_nodeb,netsnowmelt_allgl_nodeb, rain_runoff_allgl_nodeb, superimp_icemelt_allgl_nodeb,Sfc)
runoff_timeseries_12years('m3','Kaskawulsh runoff',2007,years,450,3,gl_icemelt_kw_nodeb,netsnowmelt_kw_nodeb, rain_runoff_kw_nodeb, superimp_icemelt_kw_nodeb,Sfc)
# 12 year pie charts
runoff_piecharts_12years('Catchment-wide average runoff',2007,years,gl_icemelt_krh_nodeb,netsnowmelt_krh_nodeb, rain_runoff_krh_nodeb, superimp_icemelt_krh_nodeb)
runoff_piecharts_12years('Average runoff from the glacierized area',2007,years,gl_icemelt_allgl_nodeb,netsnowmelt_allgl_nodeb, rain_runoff_allgl_nodeb, superimp_icemelt_allgl_nodeb)
runoff_piecharts_12years('Kaskawulsh runoff',2007,years,gl_icemelt_kw_nodeb,netsnowmelt_kw_nodeb, rain_runoff_kw_nodeb, superimp_icemelt_kw_nodeb)
# annual runoff stacked bar
annualrunoff_stackedbar('Kaskawulsh River Headwaters: total annual runoff',years,netsnowmelt_krh_nodeb,gl_icemelt_krh_nodeb,superimp_icemelt_krh_nodeb,rain_runoff_krh_nodeb,KRH_tributaries)
annualrunoff_stackedbar('Glacierized area: total annual runoff',years,netsnowmelt_allgl_nodeb,gl_icemelt_allgl_nodeb,superimp_icemelt_allgl_nodeb,rain_runoff_allgl_nodeb,KRH_tributaries)
annualrunoff_stackedbar('Kaskawulsh: total annual runoff',years,netsnowmelt_kw_nodeb,gl_icemelt_kw_nodeb,superimp_icemelt_kw_nodeb,rain_runoff_kw_nodeb,KRH_tributaries)
# date of zero balance
date_of_zero_balance('Catchment',years,massbal_krh_nodeb)
date_of_zero_balance('All glacierized area',years,massbal_allgl_nodeb)
date_of_zero_balance('Kaskawulsh',years,massbal_kw_nodeb)
# annual glacier contribution to catchment runoff

# cumulative mass balance
cumulative_massbalance('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_kw_nodeb, refrozen_rain_kw_nodeb, gl_icemelt_kw_nodeb, netsnowmelt_kw_nodeb, superimp_icemelt_kw_nodeb)
cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_allgl_nodeb, refrozen_rain_allgl_nodeb, gl_icemelt_allgl_nodeb, netsnowmelt_allgl_nodeb, superimp_icemelt_allgl_nodeb)
# annual mass balance bar chart
annual_massbalance_barchart('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_kw_nodeb, refrozen_rain_kw_nodeb, gl_icemelt_kw_nodeb, netsnowmelt_kw_nodeb, superimp_icemelt_kw_nodeb)
annual_massbalance_barchart('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_allgl_nodeb, refrozen_rain_allgl_nodeb, gl_icemelt_allgl_nodeb, netsnowmelt_allgl_nodeb, superimp_icemelt_allgl_nodeb)

# distributed plots
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.888283378746594,stop=1,name='massbal')
distributed_average_mass_balance('ROUNCE_DEBRIS',np.arange(2007,2017+1),years,massbal_dist_nodeb,np.linspace(-16.3,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_runoff(np.arange(1980,2021+1),years,netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_glaciericemelt(np.arange(1980,2021+1),years,netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb,np.linspace(0,17,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_snowrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb,np.linspace(0,0.7,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_rainrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb,np.linspace(0,0.3,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_SImelt(np.arange(1980,2021+1),years,netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb,np.linspace(0,0.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)


# =============================================================================
# REF_MODEL vs ROUNCE DEBRIS
# =============================================================================
# runoff differences
compare_hydrographs_differences('m3','Catchment-wide runoff: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_nodeb, netsnowmelt_krh_nodeb, rain_runoff_krh_nodeb, superimp_icemelt_krh_nodeb, gl_icemelt_krh_nodeb_std, netsnowmelt_krh_nodeb_std, rain_runoff_krh_nodeb_std, superimp_icemelt_krh_nodeb_std,Sfc)
compare_hydrographs_differences('m3','Runoff from the glacierized area: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_nodeb, netsnowmelt_allgl_nodeb, rain_runoff_allgl_nodeb, superimp_icemelt_allgl_nodeb, gl_icemelt_allgl_nodeb_std, netsnowmelt_allgl_nodeb_std, rain_runoff_allgl_nodeb_std, superimp_icemelt_allgl_nodeb_std,All_glacierized_area)
compare_hydrographs_differences('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_kw_ref, netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref_std, netsnowmelt_kw_ref_std, rain_runoff_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_nodeb, netsnowmelt_kw_nodeb, rain_runoff_kw_nodeb, superimp_icemelt_kw_nodeb, gl_icemelt_kw_nodeb_std, netsnowmelt_kw_nodeb_std, rain_runoff_kw_nodeb_std, superimp_icemelt_kw_nodeb_std,KRH_tributaries)
compare_hydrographs_differences('m3','Runoff from a debris-covered cell at the terminus: ',np.arange(1980,2021+1),years,-0.035,0.03, gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_nodeb, netsnowmelt_deb_nodeb, rain_runoff_deb_nodeb, superimp_icemelt_deb_nodeb, gl_icemelt_deb_nodeb_std, netsnowmelt_deb_nodeb_std, rain_runoff_deb_nodeb_std, superimp_icemelt_deb_nodeb_std,np.ones((1,1)))
# hydrograph comparison
compare_hydrographs('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_nodeb, netsnowmelt_allgl_nodeb, rain_runoff_allgl_nodeb, superimp_icemelt_allgl_nodeb, gl_icemelt_allgl_nodeb_std, netsnowmelt_allgl_nodeb_std, rain_runoff_allgl_nodeb_std, superimp_icemelt_allgl_nodeb_std,All_glacierized_area)
compare_hydrographs('m3','Runoff from a debris-covered cell at the terminus:\n',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_nodeb, netsnowmelt_deb_nodeb, rain_runoff_deb_nodeb, superimp_icemelt_deb_nodeb, gl_icemelt_deb_nodeb_std, netsnowmelt_deb_nodeb_std, rain_runoff_deb_nodeb_std, superimp_icemelt_deb_nodeb_std,np.ones((1,1)))
# mass balance difference:
massbalance_timeseries_difference('Catchment-wide mass balance: ',np.arange(2007,2017+1),years,0.001,0.03, accumulation_krh_ref, refrozen_rain_krh_ref, netsnowmelt_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref, accumulation_krh_ref_std, refrozen_rain_krh_ref_std, netsnowmelt_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_ref_std, accumulation_krh_nodeb, refrozen_rain_krh_nodeb, netsnowmelt_krh_nodeb, superimp_icemelt_krh_nodeb, gl_icemelt_krh_nodeb, accumulation_krh_nodeb_std, refrozen_rain_krh_nodeb_std, netsnowmelt_krh_nodeb_std, superimp_icemelt_krh_nodeb_std, gl_icemelt_krh_nodeb_std)
massbalance_timeseries_difference('Glacierized area mass balance: ',np.arange(2007,2017+1),years,0.001,0.03, accumulation_allgl_ref, refrozen_rain_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref, accumulation_allgl_ref_std, refrozen_rain_allgl_ref_std, netsnowmelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_ref_std, accumulation_allgl_nodeb, refrozen_rain_allgl_nodeb, netsnowmelt_allgl_nodeb, superimp_icemelt_allgl_nodeb, gl_icemelt_allgl_nodeb, accumulation_allgl_nodeb_std, refrozen_rain_allgl_nodeb_std, netsnowmelt_allgl_nodeb_std, superimp_icemelt_allgl_nodeb_std, gl_icemelt_allgl_nodeb_std)
massbalance_timeseries_difference('Kaskawulsh Glacier mass balance: ','ROUNCE_DEBRIS',np.arange(2007,2017+1),years,0.008,0.4, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std, accumulation_kw_nodeb, refrozen_rain_kw_nodeb, netsnowmelt_kw_nodeb, superimp_icemelt_kw_nodeb, gl_icemelt_kw_nodeb, accumulation_kw_nodeb_std, refrozen_rain_kw_nodeb_std, netsnowmelt_kw_nodeb_std, superimp_icemelt_kw_nodeb_std, gl_icemelt_kw_nodeb_std)
# mass balance comparison
massbalance_timeseries_comparison('Kaskawulsh mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std, accumulation_kw_nodeb, refrozen_rain_kw_nodeb, netsnowmelt_kw_nodeb, superimp_icemelt_kw_nodeb, gl_icemelt_kw_nodeb, accumulation_kw_nodeb_std, refrozen_rain_kw_nodeb_std, netsnowmelt_kw_nodeb_std, superimp_icemelt_kw_nodeb_std, gl_icemelt_kw_nodeb_std)
# date of zero balance 
compare_date_of_zero_balance(years,massbal_allgl_ref,massbal_allgl_nodeb)
# Distributed differences
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.21052631578947367,stop=1,name='massbal_diff')
distributed_mass_balance_difference('DEBRIS_FREE',np.arange(1980,2021+1),years,massbal_dist_ref,massbal_dist_nodeb,np.linspace(-2,7.5,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.79,stop=1,name='runoff_diff')
distributed_runoff_difference('DEBRIS_FREE',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb,np.linspace(-7.5,2,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_glaciermelt_difference('DEBRIS_FREE',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb, np.linspace(-7.5,2,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.9649122807017544,stop=1,name='snowmelt_diff')
distributed_snowrunoff_difference('DEBRIS_FREE',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb, np.linspace(-0.055,0.001,40),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=1,stop=0.7,name='rainrunoff_diff')
distributed_rainrunoff_difference('DEBRIS_FREE',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb, np.linspace(-0.003,0,10),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.5,stop=1,name='SImelt_diff')
distributed_SImelt_difference('DEBRIS_FREE',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_nodeb,gl_icemelt_dist_nodeb,superimp_icemelt_dist_nodeb,rain_runoff_dist_nodeb, np.linspace(-0.02,0.02,21),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

surface_elevation_profiles()





# =============================================================================
# Uncorrected accumulation model plots
# =============================================================================
# Average mass balance
massbalance_timeseries_average('Catchment-wide average mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_krh_uncorracc, refrozen_rain_krh_uncorracc, netsnowmelt_krh_uncorracc, superimp_icemelt_krh_uncorracc, gl_icemelt_krh_uncorracc, accumulation_krh_uncorracc_std, refrozen_rain_krh_uncorracc_std, netsnowmelt_krh_uncorracc_std, superimp_icemelt_krh_uncorracc_std, gl_icemelt_krh_uncorracc_std,Sfc)
massbalance_timeseries_average('Glacierized area mass balance: ',np.arange(2007,2018+1),years,0.035,0.8, accumulation_allgl_uncorracc, refrozen_rain_allgl_uncorracc, netsnowmelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc, gl_icemelt_allgl_uncorracc, accumulation_allgl_uncorracc_std, refrozen_rain_allgl_uncorracc_std, netsnowmelt_allgl_uncorracc_std, superimp_icemelt_allgl_uncorracc_std, gl_icemelt_allgl_uncorracc_std,All_glacierized_area)
massbalance_timeseries_average('Kaskawulsh mass balance: ',np.arange(2013,2013+1),years,0.035,1.8, accumulation_kw_uncorracc, refrozen_rain_kw_uncorracc, netsnowmelt_kw_uncorracc, superimp_icemelt_kw_uncorracc, gl_icemelt_kw_uncorracc, accumulation_kw_uncorracc_std, refrozen_rain_kw_uncorracc_std, netsnowmelt_kw_uncorracc_std, superimp_icemelt_kw_uncorracc_std, gl_icemelt_kw_uncorracc_std)
# 12 year mass balance
massbalance_timeseries_12years('Catchment-wide average mass balance',2007, years, 0.055, 1.5, accumulation_krh_uncorracc, refrozen_rain_krh_uncorracc ,netsnowmelt_krh_uncorracc, superimp_icemelt_krh_uncorracc, gl_icemelt_krh_uncorracc)
massbalance_timeseries_12years('Glacierized area mass balance',2010, years, 0.055, 1.5, accumulation_allgl_uncorracc, refrozen_rain_allgl_uncorracc ,netsnowmelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc, gl_icemelt_allgl_uncorracc)
massbalance_timeseries_12years('Kaskawulsh mass balance',2007, years, 0.055, 1.5, accumulation_kw_uncorracc, refrozen_rain_kw_uncorracc ,netsnowmelt_kw_uncorracc, superimp_icemelt_kw_uncorracc, gl_icemelt_kw_uncorracc)
massbalance_timeseries_12years('Terminus mass balance',2007, years, 0.15, 10, accumulation_deb_uncorracc, refrozen_rain_deb_uncorracc ,netsnowmelt_deb_uncorracc, superimp_icemelt_deb_uncorracc, gl_icemelt_deb_uncorracc)
# Average runoff
runoff_timeseries_average_SLRformat('m3','Catchment-wide average runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_krh_uncorracc, netsnowmelt_krh_uncorracc, rain_runoff_krh_uncorracc, superimp_icemelt_krh_uncorracc,gl_icemelt_krh_uncorracc_std, netsnowmelt_krh_uncorracc_std, rain_runoff_krh_uncorracc_std, superimp_icemelt_krh_uncorracc_std,Sfc)
runoff_timeseries_average_SLRformat('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_allgl_uncorracc, netsnowmelt_allgl_uncorracc, rain_runoff_allgl_uncorracc, superimp_icemelt_allgl_uncorracc,gl_icemelt_allgl_uncorracc_std, netsnowmelt_allgl_uncorracc_std, rain_runoff_allgl_uncorracc_std, superimp_icemelt_allgl_uncorracc_std,All_glacierized_area)
runoff_timeseries_average_SLRformat('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,275,2.75,gl_icemelt_kw_uncorracc, netsnowmelt_kw_uncorracc, rain_runoff_kw_uncorracc, superimp_icemelt_kw_uncorracc,gl_icemelt_kw_uncorracc_std, netsnowmelt_kw_uncorracc_std, rain_runoff_kw_uncorracc_std, superimp_icemelt_kw_uncorracc_std,KRH_tributaries)
runoff_timeseries_average_SLRformat('m3','Runoff at the terminus: ',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_uncorracc, netsnowmelt_deb_uncorracc, rain_runoff_deb_uncorracc, superimp_icemelt_deb_uncorracc,gl_icemelt_deb_uncorracc_std, netsnowmelt_deb_uncorracc_std, rain_runoff_deb_uncorracc_std, superimp_icemelt_deb_uncorracc_std,np.ones((1,1)))
# 12 year runoff (eventually put individual pie charts on these plots)
runoff_timeseries_12years('m3','Catchment-wide average runoff',2007,years,450,3,gl_icemelt_krh_uncorracc,netsnowmelt_krh_uncorracc, rain_runoff_krh_uncorracc, superimp_icemelt_krh_uncorracc,Sfc)
runoff_timeseries_12years('m3','Average runoff from the glacierized area',2007,years,450,3,gl_icemelt_allgl_uncorracc,netsnowmelt_allgl_uncorracc, rain_runoff_allgl_uncorracc, superimp_icemelt_allgl_uncorracc,Sfc)
runoff_timeseries_12years('m3','Kaskawulsh runoff',2007,years,450,3,gl_icemelt_kw_uncorracc,netsnowmelt_kw_uncorracc, rain_runoff_kw_uncorracc, superimp_icemelt_kw_uncorracc,Sfc)
# 12 year pie charts
runoff_piecharts_12years('Catchment-wide average runoff',2007,years,gl_icemelt_krh_uncorracc,netsnowmelt_krh_uncorracc, rain_runoff_krh_uncorracc, superimp_icemelt_krh_uncorracc)
runoff_piecharts_12years('Average runoff from the glacierized area',2007,years,gl_icemelt_allgl_uncorracc,netsnowmelt_allgl_uncorracc, rain_runoff_allgl_uncorracc, superimp_icemelt_allgl_uncorracc)
runoff_piecharts_12years('Kaskawulsh runoff',2007,years,gl_icemelt_kw_uncorracc,netsnowmelt_kw_uncorracc, rain_runoff_kw_uncorracc, superimp_icemelt_kw_uncorracc)
# annual runoff stacked bar
annualrunoff_stackedbar('Kaskawulsh River Headwaters: total annual runoff',years,netsnowmelt_krh_uncorracc,gl_icemelt_krh_uncorracc,superimp_icemelt_krh_uncorracc,rain_runoff_krh_uncorracc,KRH_tributaries)
annualrunoff_stackedbar('Glacierized area: total annual runoff',years,netsnowmelt_allgl_uncorracc,gl_icemelt_allgl_uncorracc,superimp_icemelt_allgl_uncorracc,rain_runoff_allgl_uncorracc,KRH_tributaries)
annualrunoff_stackedbar('Kaskawulsh: total annual runoff',years,netsnowmelt_kw_uncorracc,gl_icemelt_kw_uncorracc,superimp_icemelt_kw_uncorracc,rain_runoff_kw_uncorracc,KRH_tributaries)
# date of zero balance
date_of_zero_balance('Catchment',years,massbal_krh_uncorracc)
date_of_zero_balance('All glacierized area',years,massbal_allgl_uncorracc)
date_of_zero_balance('Kaskawulsh',years,massbal_kw_uncorracc)
# annual glacier contribution to catchment runoff

# cumulative mass balance
cumulative_massbalance('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_kw_uncorracc, refrozen_rain_kw_uncorracc, gl_icemelt_kw_uncorracc, netsnowmelt_kw_uncorracc, superimp_icemelt_kw_uncorracc)
cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_allgl_uncorracc, refrozen_rain_allgl_uncorracc, gl_icemelt_allgl_uncorracc, netsnowmelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc)
# annual mass balance bar chart
annual_massbalance_barchart('Kaskawulsh Glacier Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_kw_uncorracc, refrozen_rain_kw_uncorracc, gl_icemelt_kw_uncorracc, netsnowmelt_kw_uncorracc, superimp_icemelt_kw_uncorracc)
annual_massbalance_barchart('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),2.5,accumulation_allgl_uncorracc, refrozen_rain_allgl_uncorracc, gl_icemelt_allgl_uncorracc, netsnowmelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc)

annual_massbalance_barchart_5yraverages('Glacierized Area Mass Balance\n',years,np.arange(1980,2020),1.5,accumulation_allgl_uncorracc, refrozen_rain_allgl_uncorracc, gl_icemelt_allgl_uncorracc, netsnowmelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc)
annual_massbalance_barchart_10yraverages('Glacierized Area Mass Balance\n',years,np.arange(1980,2020),1.5,accumulation_allgl_uncorracc, refrozen_rain_allgl_uncorracc, gl_icemelt_allgl_uncorracc, netsnowmelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc)


# distributed plots
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.8298755186721992,stop=1,name='massbal')
distributed_average_mass_balance('UNCORRECTED_ACC',np.arange(2007,2017+1),years,massbal_dist_uncorracc,np.linspace(-10,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_runoff(np.arange(1980,2021+1),years,netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc,np.linspace(0,11,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_glaciericemelt(np.arange(1980,2021+1),years,netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc,np.linspace(0,10,18),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_snowrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc,np.linspace(0,0.7,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_rainrunoff(np.arange(1980,2021+1),years,netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc,np.linspace(0,0.3,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)
distributed_SImelt(np.arange(1980,2021+1),years,netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc,np.linspace(0,0.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)



# =============================================================================
# REF_MODEL vs UNCORRECTED_ACC
# =============================================================================
# runoff differences
compare_hydrographs_differences('m3','Catchment-wide runoff: ',np.arange(1980,2021+1),years,-125,55, gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_uncorracc, netsnowmelt_krh_uncorracc, rain_runoff_krh_uncorracc, superimp_icemelt_krh_uncorracc, gl_icemelt_krh_uncorracc_std, netsnowmelt_krh_uncorracc_std, rain_runoff_krh_uncorracc_std, superimp_icemelt_krh_uncorracc_std,Sfc)
compare_hydrographs_differences('m3','Runoff from the glacierized area: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_uncorracc, netsnowmelt_allgl_uncorracc, rain_runoff_allgl_uncorracc, superimp_icemelt_allgl_uncorracc, gl_icemelt_allgl_uncorracc_std, netsnowmelt_allgl_uncorracc_std, rain_runoff_allgl_uncorracc_std, superimp_icemelt_allgl_uncorracc_std,All_glacierized_area)
compare_hydrographs_differences('m3','Kaskawulsh runoff: ',np.arange(1980,2021+1),years,-7,10, gl_icemelt_kw_ref, netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref_std, netsnowmelt_kw_ref_std, rain_runoff_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_uncorracc, netsnowmelt_kw_uncorracc, rain_runoff_kw_uncorracc, superimp_icemelt_kw_uncorracc, gl_icemelt_kw_uncorracc_std, netsnowmelt_kw_uncorracc_std, rain_runoff_kw_uncorracc_std, superimp_icemelt_kw_uncorracc_std,KRH_tributaries)
compare_hydrographs_differences('m3','Runoff from a debris-covered cell at the terminus: ',np.arange(1980,2021+1),years,-0.035,0.03, gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_uncorracc, netsnowmelt_deb_uncorracc, rain_runoff_deb_uncorracc, superimp_icemelt_deb_uncorracc, gl_icemelt_deb_uncorracc_std, netsnowmelt_deb_uncorracc_std, rain_runoff_deb_uncorracc_std, superimp_icemelt_deb_uncorracc_std,np.ones((1,1)))
# hydrograph comparison
compare_hydrographs('m3','Catchment-wide average runoff: ','UNCORRECTED_ACC',np.arange(1980,2021+1),years,275,2.75, gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_uncorracc, netsnowmelt_krh_uncorracc, rain_runoff_krh_uncorracc, superimp_icemelt_krh_uncorracc, gl_icemelt_krh_uncorracc_std, netsnowmelt_krh_uncorracc_std, rain_runoff_krh_uncorracc_std, superimp_icemelt_krh_uncorracc_std,Sfc)
compare_hydrographs('m3','Average runoff from the glacierized area: ',np.arange(1980,2021+1),years,275,2.75, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, rain_runoff_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref_std, netsnowmelt_allgl_ref_std, rain_runoff_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_uncorracc, netsnowmelt_allgl_uncorracc, rain_runoff_allgl_uncorracc, superimp_icemelt_allgl_uncorracc, gl_icemelt_allgl_uncorracc_std, netsnowmelt_allgl_uncorracc_std, rain_runoff_allgl_uncorracc_std, superimp_icemelt_allgl_uncorracc_std,All_glacierized_area)
compare_hydrographs('m3','Runoff from a debris-covered cell at the terminus:\n',np.arange(1980,2021+1),years,0.1,0.0008,gl_icemelt_deb_ref, netsnowmelt_deb_ref, rain_runoff_deb_ref, superimp_icemelt_deb_ref, gl_icemelt_deb_ref_std, netsnowmelt_deb_ref_std, rain_runoff_deb_ref_std, superimp_icemelt_deb_ref_std, gl_icemelt_deb_uncorracc, netsnowmelt_deb_uncorracc, rain_runoff_deb_uncorracc, superimp_icemelt_deb_uncorracc, gl_icemelt_deb_uncorracc_std, netsnowmelt_deb_uncorracc_std, rain_runoff_deb_uncorracc_std, superimp_icemelt_deb_uncorracc_std,np.ones((1,1)))
# mass balance difference:
massbalance_timeseries_difference('Catchment-wide mass balance: ',np.arange(2007,2017+1),years,0.005,0.4, accumulation_krh_ref, refrozen_rain_krh_ref, netsnowmelt_krh_ref, superimp_icemelt_krh_ref, gl_icemelt_krh_ref, accumulation_krh_ref_std, refrozen_rain_krh_ref_std, netsnowmelt_krh_ref_std, superimp_icemelt_krh_ref_std, gl_icemelt_krh_ref_std, accumulation_krh_uncorracc, refrozen_rain_krh_uncorracc, netsnowmelt_krh_uncorracc, superimp_icemelt_krh_uncorracc, gl_icemelt_krh_uncorracc, accumulation_krh_uncorracc_std, refrozen_rain_krh_uncorracc_std, netsnowmelt_krh_uncorracc_std, superimp_icemelt_krh_uncorracc_std, gl_icemelt_krh_uncorracc_std)
massbalance_timeseries_difference('Glacierized area mass balance: ',np.arange(2007,2017+1),years,0.001,0.03, accumulation_allgl_ref, refrozen_rain_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref, accumulation_allgl_ref_std, refrozen_rain_allgl_ref_std, netsnowmelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_ref_std, accumulation_allgl_uncorracc, refrozen_rain_allgl_uncorracc, netsnowmelt_allgl_uncorracc, superimp_icemelt_allgl_uncorracc, gl_icemelt_allgl_uncorracc, accumulation_allgl_uncorracc_std, refrozen_rain_allgl_uncorracc_std, netsnowmelt_allgl_uncorracc_std, superimp_icemelt_allgl_uncorracc_std, gl_icemelt_allgl_uncorracc_std)
massbalance_timeseries_difference('Kaskawulsh Glacier mass balance: ','UNCORRECTED_ACC',np.arange(2007,2017+1),years,0.008,0.4, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std, accumulation_kw_uncorracc, refrozen_rain_kw_uncorracc, netsnowmelt_kw_uncorracc, superimp_icemelt_kw_uncorracc, gl_icemelt_kw_uncorracc, accumulation_kw_uncorracc_std, refrozen_rain_kw_uncorracc_std, netsnowmelt_kw_uncorracc_std, superimp_icemelt_kw_uncorracc_std, gl_icemelt_kw_uncorracc_std)
# mass balance comparison
massbalance_timeseries_comparison('Kaskawulsh mass balance: ',np.arange(2007,2017+1),years,0.035,0.8, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std, accumulation_kw_uncorracc, refrozen_rain_kw_uncorracc, netsnowmelt_kw_uncorracc, superimp_icemelt_kw_uncorracc, gl_icemelt_kw_uncorracc, accumulation_kw_uncorracc_std, refrozen_rain_kw_uncorracc_std, netsnowmelt_kw_uncorracc_std, superimp_icemelt_kw_uncorracc_std, gl_icemelt_kw_uncorracc_std)
# date of zero balance 
compare_date_of_zero_balance(years,massbal_allgl_ref,massbal_allgl_uncorracc)
# Distributed differences
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.65625,stop=1,name='massbal_diff')
distributed_mass_balance_difference('UNCORRECTED_ACC',np.arange(1980,2021+1),years,massbal_dist_ref,massbal_dist_uncorracc,np.linspace(-2.1,1.1,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.2142857142857143,stop=1,name='runoff_diff')
distributed_runoff_difference('UNCORRECTED_ACC',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc,np.linspace(-0.6,2.2,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

distributed_glaciermelt_difference('UNCORRECTED_ACC',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc, np.linspace(-0.6,2.2,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.09090909090909083,stop=1,name='snowmelt_diff')
distributed_snowrunoff_difference('UNCORRECTED_ACC',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc, np.linspace(-0.03,0.3,40),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.9375,stop=0.7,name='rainrunoff_diff')
distributed_rainrunoff_difference('UNCORRECTED_ACC',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc, np.linspace(-0.06,0.004,20),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)

shiftedColorMap(matplotlib.cm.twilight_shifted_r,start=0.1,midpoint=0.8,stop=1,name='SImelt_diff')
distributed_SImelt_difference('UNCORRECTED_ACC',np.arange(1980,2021+1),years,netsnowmelt_dist_ref, gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref, netsnowmelt_dist_uncorracc,gl_icemelt_dist_uncorracc,superimp_icemelt_dist_uncorracc,rain_runoff_dist_uncorracc, np.linspace(-0.08,0.02,41),Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries)












































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


