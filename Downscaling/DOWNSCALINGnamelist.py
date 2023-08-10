# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""
import numpy as np

#time_control
start_year = 1979
end_year = 1979
time_step = 3 #timestep in hours)

#domain
Glacier_ID = 'KRH'  #for naming output files - dont change this!!

# DOWNSCALING parameters:
UTM = '7N'

NARR_subregions_original = [8,9,10,14,15,16,20,21,22,26,27,28]
NARR_subregions_full = np.arange(0,36)
NARR_subregions_inner = [7,8,9,10,13,14,15,16,19,20,21,22,25,26,27,28]
NARR_subregions_kw = [14,15,19,20,21]
NARR_subregions_orographic = [9,10,14,15,16,19,20,21,22,25,26,27,28]

NARR_subregions = NARR_subregions_orographic #Indices of subregions used for precip downscaling, picked manually to omit points on opposite side of divide

# INPUTS
Climate_inputs = 'F:/Mass Balance Model/CoarseNARR_KW'
Coarse_DEM_input = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/kaskCE.nc' # From Datasets/NARR/time_invariant/hgt.sfc.nc
Easting_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Xgrid.txt'
Northing_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Ygrid.txt'
Elev_inputs = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Zgrids'

# OUTPUTS
OUTPUT_PATH = 'D:/Downscaled_files/Catchment/downscaling_v2_test/VectorizedDS_test1'