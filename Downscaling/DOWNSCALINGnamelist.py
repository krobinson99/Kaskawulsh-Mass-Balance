# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

#time_control
start_year = 2016
end_year = 2016
time_step = 3 # hours

#domain
Glacier_ID = 'KRH' 

# DOWNSCALING parameters:
UTM = 7 # UTM zone

NARR_subregions = [9,10,14,15,16,19,20,21,22,25,26,27,28] #Indices of subregions used for precip downscaling, picked manually to omit points on opposite side of divide

# INPUTS
Model_functions = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel'
Climate_inputs = 'F:/Mass Balance Model/CoarseNARR_KW'
Coarse_DEM_input = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/kaskCE.nc' # From Datasets/NARR/time_invariant/hgt.sfc.nc
Easting_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Xgrid.txt'
Northing_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Ygrid.txt'
Elev_inputs = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Zgrids'

# OUTPUTS
OUTPUT_PATH = 'D:/Downscaled_files/Catchment/downscaling_v2_test/VectorizedDS_test1'