# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:32:38 2023

This script takes all the outputs from a set of model runs and 
returns the average of all runs for each variable.

@author: katierobinson
"""

import numpy as np
import os
import sys
import pandas as pd
from netCDF4 import Dataset

from MBMnamelist import Easting_grid, Northing_grid, Sfc_grid

# Import functions for model:
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import save_to_netcdf

# =============================================================================
# User specified options:
# =============================================================================

INPUT_FOLDER = 'D:/TuningOutputs/Tuning_Sept26'
OUTPUT_FOLDER = 'D:/TuningOutputs/Tuning_Sept26'
OUTPUT_ID = 99999

sims = np.arange(0,100)

# Get sim param from the job array
# =============================================================================
year = 1979
#year = int(sys.argv[1])
#print('this is year',year)
#print(type(year))

Glacier_ID = 'KRH'

# =============================================================================
# Load input geometry:
# =============================================================================

Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
Sfc = np.loadtxt(Sfc_grid)

# Function to concatenate any variable:
def concatenate_model_ouputs(year, sims, varname, var):
    '''
    Inputs: 
        year (1979--2022)
        sims (array or list of integers corresponding to model runs to be averaged)
        varname (string, e.g. 'Icemelt')
        var: netcdf variable (string, e.g. 'Ice melt')
    Returns:
        Average and standard deviation of all model runs for a given year and sim number. 
    '''
    
    print('year:',year)
    sys.stdout.flush()
    dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')  
    
# =============================================================================
#   Calculate mean of all sims:
# =============================================================================
    all_sims = np.empty((len(sims),len(dates),Xgrid.shape[0],Xgrid.shape[1]))
    for sim in sims:
        print('Calculating mean and std dev. for sim:',sim)
        sys.stdout.flush()
        
        inMB = Dataset(os.path.join(INPUT_FOLDER,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
        mb = inMB.variables[var][:]
        sys.stdout.flush()
        inMB.close()
        
        # Combine data from every sim into one array
        all_sims[sim,:,:,:] = mb
        
    # Get the mean and std over axis 0 (the 100 sims), outputs should have shape (time,x,y)
    mean_mb = np.nanmean(all_sims,axis=0)
    std = np.std(all_sims,axis=0)
        
    save_to_netcdf(mean_mb, var, os.path.join(OUTPUT_FOLDER,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(OUTPUT_ID) + '.nc'), year, Xgrid, Ygrid) 
    save_to_netcdf(std, var, os.path.join(OUTPUT_FOLDER,varname + '_' + 'std' + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(OUTPUT_ID) + '.nc'), year, Xgrid, Ygrid) 
    
    
#Call function for each of the model outputs:    
# =============================================================================
concatenate_model_ouputs(year,sims,'Snowmelt', 'Snowmelt')   
concatenate_model_ouputs(year,sims,'Refrozenmelt', 'Refrozenmelt')   
concatenate_model_ouputs(year,sims,'Netsnowmelt', 'Netsnowmelt')  
concatenate_model_ouputs(year,sims,'Glaciericemelt', 'Glaciericemelt')   
concatenate_model_ouputs(year,sims,'Superimposedicemelt', 'Superimposedicemelt')  
concatenate_model_ouputs(year,sims,'Rain', 'Rain')    
concatenate_model_ouputs(year,sims,'Refrozenrain', 'Refrozenrain')  
concatenate_model_ouputs(year,sims,'Rainrunoff', 'Rainrunoff') 
concatenate_model_ouputs(year,sims,'Accumulation', 'Accumulation') 
concatenate_model_ouputs(year,sims,'Snowdepth', 'Snowdepth') 
concatenate_model_ouputs(year,sims,'Ptau', 'Ptau')
concatenate_model_ouputs(year,sims,'Superimposedice', 'Superimposedice')  
concatenate_model_ouputs(year,sims,'Netbalance', 'Netbalance')
# =============================================================================

  