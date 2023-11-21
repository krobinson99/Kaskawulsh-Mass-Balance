# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:13:38 2023

Process outputs from the mass balance model
Save daily timeseries as text files to make plotting more efficient.

@author: katierobinson
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(2,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\ProcessOutputs')
from MBM_plottingfunctions import calculate_mb_components_timeseries, \
calculate_stddev_timeseries, calculate_mb_components_distributed, full_save


MODEL_OUTPUTS_REFMODEL = 'D:/Model Runs/REF_MODEL/Sim_99999_v2'
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

varnames = ['massbal','totalsnowmelt','refreezing','netsnowmelt','gl_icemelt','superimp_icemelt','rain','refrozen_rain','rain_runoff','accumulation','snowdepth','Ptau','SI','temp']


# =============================================================================
#  REF_MODEL: Catchment-wide average
# =============================================================================
massbal_krh_ref, totalsnowmelt_krh_ref, refreezing_krh_ref, netsnowmelt_krh_ref, gl_icemelt_krh_ref, superimp_icemelt_krh_ref, rain_krh_ref, refrozen_rain_krh_ref, rain_runoff_krh_ref, accumulation_krh_ref, \
snowdepth_krh_ref, Ptau_krh_ref, SI_krh_ref, temp_krh_ref = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

full_save([massbal_krh_ref, totalsnowmelt_krh_ref, refreezing_krh_ref, netsnowmelt_krh_ref, gl_icemelt_krh_ref, superimp_icemelt_krh_ref, rain_krh_ref, refrozen_rain_krh_ref, rain_runoff_krh_ref, accumulation_krh_ref, \
snowdepth_krh_ref, Ptau_krh_ref, SI_krh_ref, temp_krh_ref],years,MODEL_OUTPUTS_REFMODEL,'REF_MODEL',varnames,'krh')

massbal_krh_ref_std, totalsnowmelt_krh_ref_std, refreezing_krh_ref_std, netsnowmelt_krh_ref_std, gl_icemelt_krh_ref_std, superimp_icemelt_krh_ref_std, rain_krh_ref_std, refrozen_rain_krh_ref_std, rain_runoff_krh_ref_std, accumulation_krh_ref_std, \
snowdepth_krh_ref_std, Ptau_krh_ref_std, SI_krh_ref_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

full_save([massbal_krh_ref_std, totalsnowmelt_krh_ref_std, refreezing_krh_ref_std, netsnowmelt_krh_ref_std, gl_icemelt_krh_ref_std, superimp_icemelt_krh_ref_std, rain_krh_ref_std, refrozen_rain_krh_ref_std, rain_runoff_krh_ref_std, accumulation_krh_ref_std, \
snowdepth_krh_ref_std, Ptau_krh_ref_std, SI_krh_ref_std],years,MODEL_OUTPUTS_REFMODEL,'REF_MODEL',varnames,'krh_std')
# =============================================================================
#  REF_MODEL: All glacierized area
# =============================================================================
massbal_allgl_ref, totalsnowmelt_allgl_ref, refreezing_allgl_ref, netsnowmelt_allgl_ref, gl_icemelt_allgl_ref, superimp_icemelt_allgl_ref, rain_allgl_ref, refrozen_rain_allgl_ref, rain_runoff_allgl_ref, accumulation_allgl_ref, \
snowdepth_allgl_ref, Ptau_allgl_ref, SI_allgl_ref, temp_allgl_ref = calculate_mb_components_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_allgl_ref, totalsnowmelt_allgl_ref, refreezing_allgl_ref, netsnowmelt_allgl_ref, gl_icemelt_allgl_ref, superimp_icemelt_allgl_ref, rain_allgl_ref, refrozen_rain_allgl_ref, rain_runoff_allgl_ref, accumulation_allgl_ref, \
snowdepth_allgl_ref, Ptau_allgl_ref, SI_allgl_ref, temp_allgl_ref],years,MODEL_OUTPUTS_REFMODEL,'REF_MODEL',varnames,'allgl')

massbal_allgl_ref_std, totalsnowmelt_allgl_ref_std, refreezing_allgl_ref_std, netsnowmelt_allgl_ref_std, gl_icemelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, rain_allgl_ref_std, refrozen_rain_allgl_ref_std, rain_runoff_allgl_ref_std, accumulation_allgl_ref_std, \
snowdepth_allgl_ref_std, Ptau_allgl_ref_std, SI_allgl_ref_std = calculate_stddev_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_allgl_ref_std, totalsnowmelt_allgl_ref_std, refreezing_allgl_ref_std, netsnowmelt_allgl_ref_std, gl_icemelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, rain_allgl_ref_std, refrozen_rain_allgl_ref_std, rain_runoff_allgl_ref_std, accumulation_allgl_ref_std, \
snowdepth_allgl_ref_std, Ptau_allgl_ref_std, SI_allgl_ref_std],years,MODEL_OUTPUTS_REFMODEL,'REF_MODEL',varnames,'allgl_std')
# =============================================================================
#  REF MODEL: Kaskawulsh
# =============================================================================
massbal_kw_ref, totalsnowmelt_kw_ref, refreezing_kw_ref, netsnowmelt_kw_ref, gl_icemelt_kw_ref, superimp_icemelt_kw_ref, rain_kw_ref, refrozen_rain_kw_ref, rain_runoff_kw_ref, accumulation_kw_ref, \
snowdepth_kw_ref, Ptau_kw_ref, SI_kw_ref, temp_kw_ref = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_kw_ref, totalsnowmelt_kw_ref, refreezing_kw_ref, netsnowmelt_kw_ref, gl_icemelt_kw_ref, superimp_icemelt_kw_ref, rain_kw_ref, refrozen_rain_kw_ref, rain_runoff_kw_ref, accumulation_kw_ref, \
snowdepth_kw_ref, Ptau_kw_ref, SI_kw_ref, temp_kw_ref],years,MODEL_OUTPUTS_REFMODEL,'REF_MODEL',varnames,'kw')

massbal_kw_ref_std, totalsnowmelt_kw_ref_std, refreezing_kw_ref_std, netsnowmelt_kw_ref_std, gl_icemelt_kw_ref_std, superimp_icemelt_kw_ref_std, rain_kw_ref_std, refrozen_rain_kw_ref_std, rain_runoff_kw_ref_std, accumulation_kw_ref_std, \
snowdepth_kw_ref_std, Ptau_kw_ref_std, SI_kw_ref_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_kw_ref_std, totalsnowmelt_kw_ref_std, refreezing_kw_ref_std, netsnowmelt_kw_ref_std, gl_icemelt_kw_ref_std, superimp_icemelt_kw_ref_std, rain_kw_ref_std, refrozen_rain_kw_ref_std, rain_runoff_kw_ref_std, accumulation_kw_ref_std, \
snowdepth_kw_ref_std, Ptau_kw_ref_std, SI_kw_ref_std],years,MODEL_OUTPUTS_REFMODEL,'REF_MODEL',varnames,'kw_std')

# =============================================================================
#  REF MODEL: Debris-covered gridcell at the terminus
# =============================================================================
massbal_deb_ref, totalsnowmelt_deb_ref, refreezing_deb_ref, netsnowmelt_deb_ref, gl_icemelt_deb_ref, superimp_icemelt_deb_ref, rain_deb_ref, refrozen_rain_deb_ref, rain_runoff_deb_ref, accumulation_deb_ref, \
snowdepth_deb_ref, Ptau_deb_ref, SI_deb_ref, temp_deb_ref = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,13,314,PointScale=False,KWaverage=True,Catchmentaverage=False)

massbal_deb_ref_std, totalsnowmelt_deb_ref_std, refreezing_deb_ref_std, netsnowmelt_deb_ref_std, gl_icemelt_deb_ref_std, superimp_icemelt_deb_ref_std, rain_deb_ref_std, refrozen_rain_deb_ref_std, rain_runoff_deb_ref_std, accumulation_deb_ref_std, \
snowdepth_deb_ref_std, Ptau_deb_ref_std, SI_deb_ref_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID,13,314,PointScale=True,KWaverage=False,Catchmentaverage=False)

# =============================================================================
#  REF_MODEL: Distributed totals
# =============================================================================
massbal_dist_ref, totalsnowmelt_dist_ref, refreezing_dist_ref, netsnowmelt_dist_ref, gl_icemelt_dist_ref, superimp_icemelt_dist_ref, rain_dist_ref, refrozen_rain_dist_ref, rain_runoff_dist_ref, accumulation_dist_ref, \
snowdepth_dist_ref, Ptau_dist_ref, SI_dist_ref, temp_dist_ref = calculate_mb_components_distributed(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_REFMODEL,Glacier_ID)










# =============================================================================
#  ROUNCE_DEBRIS: Catchment-wide average
# =============================================================================
massbal_krh, totalsnowmelt_krh, refreezing_krh, netsnowmelt_krh, gl_icemelt_krh, superimp_icemelt_krh, rain_krh, refrozen_rain_krh, rain_runoff_krh, accumulation_krh, \
snowdepth_krh, Ptau_krh, SI_krh, temp_krh = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

full_save([massbal_krh, totalsnowmelt_krh, refreezing_krh, netsnowmelt_krh, gl_icemelt_krh, superimp_icemelt_krh, rain_krh, refrozen_rain_krh, rain_runoff_krh, accumulation_krh, \
snowdepth_krh, Ptau_krh, SI_krh, temp_krh],years,MODEL_OUTPUTS_ROUNCEDEBRIS,'ROUNCE_DEBRIS',varnames,'krh')

massbal_krh_std, totalsnowmelt_krh_std, refreezing_krh_std, netsnowmelt_krh_std, gl_icemelt_krh_std, superimp_icemelt_krh_std, rain_krh_std, refrozen_rain_krh_std, rain_runoff_krh_std, accumulation_krh_std, \
snowdepth_krh_std, Ptau_krh_std, SI_krh_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=False,Catchmentaverage=True)

full_save([massbal_krh_std, totalsnowmelt_krh_std, refreezing_krh_std, netsnowmelt_krh_std, gl_icemelt_krh_std, superimp_icemelt_krh_std, rain_krh_std, refrozen_rain_krh_std, rain_runoff_krh_std, accumulation_krh_std, \
snowdepth_krh_std, Ptau_krh_std, SI_krh_std],years,MODEL_OUTPUTS_ROUNCEDEBRIS,'ROUNCE_DEBRIS',varnames[:-1],'krh_std')

# =============================================================================
#  ROUNCE_DEBRIS: All glacierized area
# =============================================================================
massbal_allgl, totalsnowmelt_allgl, refreezing_allgl, netsnowmelt_allgl, gl_icemelt_allgl, superimp_icemelt_allgl, rain_allgl, refrozen_rain_allgl, rain_runoff_allgl, accumulation_allgl, \
snowdepth_allgl, Ptau_allgl, SI_allgl, temp_allgl = calculate_mb_components_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_allgl, totalsnowmelt_allgl, refreezing_allgl, netsnowmelt_allgl, gl_icemelt_allgl, superimp_icemelt_allgl, rain_allgl, refrozen_rain_allgl, rain_runoff_allgl, accumulation_allgl, \
snowdepth_allgl, Ptau_allgl, SI_allgl, temp_allgl],years,MODEL_OUTPUTS_ROUNCEDEBRIS,'ROUNCE_DEBRIS',varnames,'allgl')

massbal_allgl_std, totalsnowmelt_allgl_std, refreezing_allgl_std, netsnowmelt_allgl_std, gl_icemelt_allgl_std, superimp_icemelt_allgl_std, rain_allgl_std, refrozen_rain_allgl_std, rain_runoff_allgl_std, accumulation_allgl_std, \
snowdepth_allgl_std, Ptau_allgl_std, SI_allgl_std = calculate_stddev_timeseries(sim,years,R2S,All_glacierized_area,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_allgl_std, totalsnowmelt_allgl_std, refreezing_allgl_std, netsnowmelt_allgl_std, gl_icemelt_allgl_std, superimp_icemelt_allgl_std, rain_allgl_std, refrozen_rain_allgl_std, rain_runoff_allgl_std, accumulation_allgl_std, \
snowdepth_allgl_std, Ptau_allgl_std, SI_allgl_std],years,MODEL_OUTPUTS_ROUNCEDEBRIS,'ROUNCE_DEBRIS',varnames[:-1],'allgl_std')

# =============================================================================
#  ROUNCE_DEBRIS: Kaskawulsh
# =============================================================================
massbal_kw, totalsnowmelt_kw, refreezing_kw, netsnowmelt_kw, gl_icemelt_kw, superimp_icemelt_kw, rain_kw, refrozen_rain_kw, rain_runoff_kw, accumulation_kw, \
snowdepth_kw, Ptau_kw, SI_kw, temp_kw = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_kw, totalsnowmelt_kw, refreezing_kw, netsnowmelt_kw, gl_icemelt_kw, superimp_icemelt_kw, rain_kw, refrozen_rain_kw, rain_runoff_kw, accumulation_kw, \
snowdepth_kw, Ptau_kw, SI_kw, temp_kw],years,MODEL_OUTPUTS_ROUNCEDEBRIS,'ROUNCE_DEBRIS',varnames,'kw')

massbal_kw_std, totalsnowmelt_kw_std, refreezing_kw_std, netsnowmelt_kw_std, gl_icemelt_kw_std, superimp_icemelt_kw_std, rain_kw_std, refrozen_rain_kw_std, rain_runoff_kw_std, accumulation_kw_std, \
snowdepth_kw_std, Ptau_kw_std, SI_kw_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,0,0,PointScale=False,KWaverage=True,Catchmentaverage=False)

full_save([massbal_kw_std, totalsnowmelt_kw_std, refreezing_kw_std, netsnowmelt_kw_std, gl_icemelt_kw_std, superimp_icemelt_kw_std, rain_kw_std, refrozen_rain_kw_std, rain_runoff_kw_std, accumulation_kw_std, \
snowdepth_kw_std, Ptau_kw_std, SI_kw_std],years,MODEL_OUTPUTS_ROUNCEDEBRIS,'ROUNCE_DEBRIS',varnames[:-1],'kw_std')

# =============================================================================
# ROUNCE_DEBRIS Debris-covered gridcell at the terminus
# =============================================================================
massbal_deb, totalsnowmelt_deb, refreezing_deb, netsnowmelt_deb, gl_icemelt_deb, superimp_icemelt_deb, rain_deb, refrozen_rain_deb, rain_runoff_deb, accumulation_deb, \
snowdepth_deb, Ptau_deb, SI_deb, temp_deb = calculate_mb_components_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,13,314,PointScale=True,KWaverage=False,Catchmentaverage=False)

massbal_deb_std, totalsnowmelt_deb_std, refreezing_deb_std, netsnowmelt_deb_std, gl_icemelt_deb_std, superimp_icemelt_deb_std, rain_deb_std, refrozen_rain_deb_std, rain_runoff_deb_std, accumulation_deb_std, \
snowdepth_deb_std, Ptau_deb_std, SI_deb_std = calculate_stddev_timeseries(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID,13,314,PointScale=True,KWaverage=False,Catchmentaverage=False)

# =============================================================================
#  ROUNCE_DEBRIS: Distributed totals
# =============================================================================
massbal_dist, totalsnowmelt_dist, refreezing_dist, netsnowmelt_dist, gl_icemelt_dist, superimp_icemelt_dist, rain_dist, refrozen_rain_dist, rain_runoff_dist, accumulation_dist, \
snowdepth_dist, Ptau_dist, SI_dist, temp_dist = calculate_mb_components_distributed(sim,years,R2S,KRH_tributaries,NARR_INPUTS,MODEL_OUTPUTS_ROUNCEDEBRIS,Glacier_ID)




