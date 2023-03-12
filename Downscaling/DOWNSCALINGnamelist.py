# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

#time_control
start_year = 2010
end_year = 2010
delta_t = 3 #timestep in hours)

#domain
glacier_id = 'kaskonly'  #for naming output files - dont change this!!
glacier_outline = 'kask_catchment' #change this to import a different outline file
catchment = True #are we looking at the kaskawulsh only (False) or the whole catchment? (True)
static_surface = True

# DOWNSCALING parameters:
Topo_param = False #False = on ice grid, True = off ice grid #eventually remove this param also - doesnt do anything 
Debris = False #false = debris IS considered, true = debris NOT considered (Leave false for downscaling - eventually remove the option)
snow_start = False #True to set snow initial condition to 0 in first year, False to carry over snow from previous year
save_snow = False #Followed by boolean for whether to generate snow.txt output (will overwrite previous files, leave FALSE to reuse snow from previous downscaling run)
UTM = '7N'
#ELA = 1900
NARR_subregions = [8,9,10,14,15,16,20,21,22,26,27,28] #Indices of subregions used for precip downscaling, picked manually to omit points on opposite side of divide

# Bias correction
D_T = False # LEAVE OFF - BIAS CORRECT SEPERATELY
D_P = False

#physics
cfactor = 1.
snowfactor = 1. #(rain to snow threshold)

# PATHS for inputs/outputs
solar_in = 'F:\Mass Balance Model\Solar_Inputs'
rawnarr_inputs = 'F:/Mass Balance Model/CoarseNARR_KW'
downscaled_outputs = 'D:/Downscaled_files/DynamicSurface/downscaling_v2_test'