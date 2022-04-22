# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

#time_control
start_year = 2011
end_year = 2011      
delta_t = 3 #timestep in hours)

#domain
glacier_id = 'kaskonly'       #points to a specific glacier outline - dont change this!!
Glacier_outline = 'kaskonly' #change this to import a different outline file
considering_kaskonly = True #are we looking at the kaskawulsh only (True) or the whole catchment? (False)

# DOWNSCALING parameters:
Topo_param = False #False = on ice grid, True = off ice grid
Debris = False #false = debris IS considered, true = debris NOT considered (?)
SNOW_START = True #usually FALSE
gen_snow = True
UTM = '7N'
ELA = 1900
NARR_subregions = [8,9,10,14,15,16,20,21,22,26,27,28] #Indices of subregions used for precip downscaling, picked manually to omit points on opposite side of divide

# Bias correction
D_T = False
D_P = False

#physics
cfactor = 1.
snowfactor = 1.

# PATHS for inputs/outputs
downscaled_outputs = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'