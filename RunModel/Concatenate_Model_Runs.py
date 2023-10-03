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
# Import functions for model:
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import save_to_netcdf

INPUT_FOLDER = 'D:/TuningOutputs/Tuning_Sept26'
OUTPUT_FOLDER = 'D:/TuningOutputs/Tuning_Sept26'
OUTPUT_ID = 9999

years = np.arange(1979,2022+1)
sims = np.arange(0,2)

Glacier_ID = 'KRH'
R2S = 1

# Input geometry
Easting_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Xgrid.txt' # Paths to text files defining Easting/Northing coords of every model gridcell
Northing_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Ygrid.txt'
Sfc_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_SfcType.txt'      # Path to text file where 1 = off-glacier, 0 = on-glacier, NaN = not in the domain.

Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
Sfc = np.loadtxt(Sfc_grid)

# Function to concatenate any variable:
def concatenate_model_ouputs(varname, var):
    
    for year in years:
        print('year:',year)
        
        dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq=str(3)+'H')
        MB_sum = np.zeros((len(dates),Xgrid.shape[0],Xgrid.shape[1]))
        
        for sim in sims:
            print('sim:',sim)
            
            inMB = Dataset(os.path.join(INPUT_FOLDER,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
            mb_array = inMB.variables[var][:]
            sys.stdout.flush()
            inMB.close()
        
            MB_sum += np.array(mb_array)
            
        Average_MB = np.array(MB_sum)/len(sims)
        save_to_netcdf(Average_MB, var, os.path.join(OUTPUT_FOLDER,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(OUTPUT_ID) + '.nc'), year, Xgrid, Ygrid) 
        
concatenate_model_ouputs('Icemelt', 'Ice melt')   
concatenate_model_ouputs('Snowmelt', 'Snow melt')   
concatenate_model_ouputs('Refrozenmelt', 'Refreezing melt')  
concatenate_model_ouputs('Refrozenrain', 'Refreezing rain')   
concatenate_model_ouputs('Netbalance', 'Net balance')  
concatenate_model_ouputs('Superimposedice', 'Superimposed ice')    
concatenate_model_ouputs('Snowdepth', 'Snow depth')  
concatenate_model_ouputs('Ptau', 'Ptau')    

  