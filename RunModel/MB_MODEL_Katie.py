# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:02:43 2021

Revised version of glacier mass balance model (originally by Young et al. 2021)
Calculates distributed glacier mass balance and runoff.

Set configuration options in MBMnamelist.py

Inputs:
    Melt parameters (MF, asnow, aice)
    Downscaled/bias-corrected NARR 3-hourly surface air temperature
    Downscaled/bias-corrected NARR 3-hourly surface precipitation
    3 hourly potential direct incoming solar radiation
    Coordinates for model grid (X,Y)
    Surface Type Grid (1 = off-glacier, 0 = on-glacier, NaN = not in the domain.)
    Debris cover grid
Outputs:
    Distributed snow melt, ice melt, refreezing, and mass balance as NetCDF files.

@author: katierobinson
"""

# Import packages:
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import sys
import os

# Import parameters from config file
from MBMnamelist import start_year, end_year, timestep, params_file, sim, R2S, Glacier_ID, \
Model_functions, Precip_inputs, Temp_inputs, Solar_inputs, Easting_grid, Northing_grid, \
Sfc_grid, debris_parameterization, debris_map, cleanice_melt, peak_melt, peak_melt_thickness, \
transition_thickness, b0, k, OUTPUT_PATH, SaveMBonly

# Import functions for model:
sys.path.insert(1,Model_functions)
from Model_functions_ver4 import write_config_file, save_to_netcdf
from Model_functions_ver4 import debris, get_meanSP, cold_content, MassBalance

# Save configuration file for this run to output directory:
write_config_file(OUTPUT_PATH,"MBMnamelist.py")

# Load model grid:
# =============================================================================
years = np.arange(start_year,end_year+1)
Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
Sfc = np.loadtxt(Sfc_grid)
print('Model coordinates loaded.')
sys.stdout.flush()
# =============================================================================

# Load melt parameters for this run:
# =============================================================================
# KR_note: params file should have 3 columns (aice, asnow, MF) with one row per simulation    
params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice, asnow, MF = params[sim]
print('Parameters for this run:\naice = ',aice,'\nasnow = ',asnow,'\nMF = ',MF)
sys.stdout.flush()
# =============================================================================
  
# Set up the debris parameterization:
# =============================================================================
debris_m = debris(debris_parameterization,debris_map,Sfc,cleanice_melt,peak_melt,peak_melt_thickness,transition_thickness,b0,k)
print('Debris parameterization loaded with option:',debris_parameterization)
sys.stdout.flush()
# =============================================================================

# Calculate mean annual snowpack (for cold content/refreezing):
# =============================================================================
Cmean = get_meanSP(years,Glacier_ID,R2S,Precip_inputs,Temp_inputs)
print('Mean annual snowpack calculated') 
sys.stdout.flush()
# =============================================================================

for year in years:
    print('Running model for',year)
    sys.stdout.flush()
    dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq=str(timestep)+'H')
    
    # Load inputs for model (Temperature, Precip, Solar Radiation) 
    # =========================================================================
    
    # Temperature:
    inT = Dataset(os.path.join(Temp_inputs,'Temperature_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    T_array = inT.variables['Temperature'][:]
    sys.stdout.flush()
    
    # Precipitation:
    inP = Dataset(os.path.join(Precip_inputs,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    P_array = inP.variables['Precipitation'][:]
    sys.stdout.flush()
    
    # Solar radiation:
    inS = Dataset(os.path.join(Solar_inputs,'Solar_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    S_array = inS.variables['SolarRadiation'][:]
    sys.stdout.flush()
    
    print(year,'inputs (temperature, precipitation, solar) loaded.') 
    sys.stdout.flush()
    # =========================================================================     

    # Get cold content for year:
    # =========================================================================
    if year == years[-1]:
        pass
    else:
        CC = cold_content(year, P_array,T_array, Glacier_ID, Cmean, R2S, Precip_inputs, Temp_inputs)

    # Set up output files:
    # =========================================================================
    TotalMelt = np.empty(T_array.shape)
    IceMelt = np.empty(T_array.shape)
    SnowMelt = np.empty(T_array.shape)
    RefrozenMelt = np.empty(T_array.shape)
    MassBal = np.empty(T_array.shape)
    
    if year == years[0]:
        # Trackers should have +1 extra timestep to account for carry over values for following year:
        Snowpack_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        CC_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
    else:
        # Load snowpack and CC from previous year here.
        Snowpack_tracker[0] = Snowpack_tracker[-1] # Carry over snowpack data from end of previous year to beginning of new year
        Snowpack_tracker[1:] = 0 # Reset snowpack for rest of the year to zero

        CC_tracker[0] = CC_tracker[-1]
        CC_tracker[1:] = 0

    # Set up timestepping (loop through every timestep in year)
    # =========================================================================
    for timestamp in range(0,len(dates)):
        print(dates[timestamp])
        sys.stdout.flush()
        
        # Update cold content and snowpack trackers
        # =====================================================================
        # (from original model): add cold content in late october to snow pack to prevent winter melt events
        if ((dates[timestamp] >= pd.Timestamp(str(year)+'-10-01T00')) and (dates[timestamp] < pd.Timestamp(str(year)+'-10-02T00'))):
            CC_tracker[timestamp] += CC/8
        else:
            pass
        
        New_snowfall =  P_array[timestamp,:,:]
        New_snowfall[np.where(T_array[timestamp,:,:] > R2S)] = 0
        Snowpack_tracker[timestamp,:,:] += New_snowfall
        
        # Calculate Melt:
        # =====================================================================        
        Msnow, Mice, Refreezing, CC_out, SP_out = MassBalance(MF,asnow,aice,T_array[timestamp,:,:],S_array[timestamp,:,:],Snowpack_tracker[timestamp,:,:],CC_tracker[timestamp,:,:],debris_m,debris_parameterization,Sfc)
        # KR_note: check that this is working as expected
        
        # Update output arrays for this timestep:
        # Total Melt = Snow Melt + Ice Melt
        # Net ablation = total melt - refreezing
        # Net balance = Accumulation - net ablation
        # ===================================================================== 
        IceMelt[timestamp,:,:] = Mice
        SnowMelt[timestamp,:,:] = Msnow
        RefrozenMelt[timestamp,:,:] = Refreezing
        MassBal[timestamp,:,:] = New_snowfall - ((Msnow + Mice) - Refreezing)  

        # Update snowpack and CC trackers for next timestep:
        # ===================================================================== 
        Snowpack_tracker[timestamp+1,:,:] = SP_out
        CC_tracker[timestamp+1,:,:] = CC_out   
        
    # Save outputs for the year before starting next year:
    # =========================================================================
    print('Saving model outputs for',year)
    if SaveMBonly == False:
        save_to_netcdf(IceMelt, 'Ice melt', os.path.join(OUTPUT_PATH,'Icemelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 
        save_to_netcdf(SnowMelt, 'Snow melt', os.path.join(OUTPUT_PATH,'Snowmelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 
        save_to_netcdf(RefrozenMelt, 'Refreezing', os.path.join(OUTPUT_PATH,'Refreezing_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 
    else:
        pass
    
    save_to_netcdf(MassBal, 'Net balance', os.path.join(OUTPUT_PATH,'Netbalance_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 

    inT.close()
    inP.close()
    inS.close()

print('RUN COMPLETE!')