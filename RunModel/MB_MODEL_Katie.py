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
from Model_functions_ver4 import debris, Calculate_Pmean, max_superimposed_ice,  \
max_superimposed_ice_finalyear,MassBalance

# Save configuration file for this run to output directory:
write_config_file(OUTPUT_PATH,"MBMnamelist.py")

# Get sim param from the job array
# =============================================================================
#sim = int(sys.argv[1])
#print('this is sim',sim)
#print(type(sim))

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
params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice, asnow, MF = params[sim][1:]
# Save params for this run in a text file just in case:
np.savetxt(os.path.join(OUTPUT_PATH,'sim' + str(sim) + '_params.txt'),params[sim][1:])
print('Parameters for this run:\naice = ',aice,'\nasnow = ',asnow,'\nMF = ',MF)
sys.stdout.flush()
# =============================================================================
  
# Set up the debris parameterization:
# =============================================================================
debris_m = debris(debris_parameterization,debris_map,Sfc,cleanice_melt,peak_melt,peak_melt_thickness,transition_thickness,b0,k)
print('Debris parameterization loaded with option:',debris_parameterization)
sys.stdout.flush()
# =============================================================================

# Calculate mean annual precipitation (for refreezing):
# =============================================================================
Pmean = Calculate_Pmean(years,Glacier_ID,Precip_inputs,Temp_inputs,Sfc)
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

    # Get maximum amount of superimposed ice that can form each hydrologic year:
    # =========================================================================
    if year == years[-1]:
        SImax = max_superimposed_ice_finalyear(year, P_array, T_array, timestep, Pmean)
    else:
        SImax = max_superimposed_ice(year, P_array, T_array, timestep, Glacier_ID, Pmean, Precip_inputs, Temp_inputs)

    # Set up output files and trackers for snowpack, potential superimposed ice:
    # =========================================================================
    IceMelt = np.zeros(T_array.shape)
    SnowMelt = np.zeros(T_array.shape)
    RefrozenMelt = np.zeros(T_array.shape)
    MassBal = np.zeros(T_array.shape)
    
    if year == years[0]:
        # Trackers should have +1 extra timestep to account for carry over values for following year:
        Snowpack_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        PotentialSI_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
    else:
        # Save final snowpack and leftover potential superimposed ice from previous year here.
        Snowpack_carryover = np.array(Snowpack_tracker[-1])
        PotentialSI_carryover = np.array(PotentialSI_tracker[-1])
        
        # Reset snowpack for rest of the year to zero
        Snowpack_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        PotentialSI_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        
        # Add carry over values from previous year to beginning of tracker
        Snowpack_tracker[0] = Snowpack_carryover
        PotentialSI_tracker[0] = PotentialSI_carryover

    # Set up timestepping (loop through every timestep in year)
    # =========================================================================
    for timestamp in range(0,len(dates)):
        print(dates[timestamp])
        sys.stdout.flush()
        
        # Get current snowpack, potential SI, and actual SI volumes for this timestep:
        # =====================================================================
        # Renew cold content at beginning of hydrological year:
        if dates[timestamp] == pd.Timestamp(str(year)+'-10-01T00'):
            PotentialSI_tracker[timestamp] += SImax
        else:
            pass
        
        New_snowfall =  np.array(P_array[timestamp,:,:])
        New_snowfall[np.where(T_array[timestamp,:,:] > R2S)] = 0
        Snowpack_tracker[timestamp,:,:] += New_snowfall
        
        Curr_superimposedice = (np.cumsum(RefrozenMelt[:timestamp+1],axis=0)[-1] + np.cumsum(IceMelt[:timestamp+1],axis=0)[-1])
        Curr_superimposedice[np.where(Curr_superimposedice < 0)] = 0
        
        # Calculate Melt:
        # =====================================================================        
        Msnow, Mice, Refreezing, SI_out, SP_out = MassBalance(MF,asnow,aice,T_array[timestamp,:,:],S_array[timestamp,:,:],Snowpack_tracker[timestamp,:,:],PotentialSI_tracker[timestamp,:,:],debris_m,debris_parameterization,Sfc,Curr_superimposedice)
        
        # Calculate refreezing of rain (if any)
        # =====================================================================  
        Rain = np.array(P_array[timestamp,:,:])
        Rain[np.where(T_array[timestamp,:,:] <= R2S)] = 0 # zero the snow fall
    
        RefrozenRain = np.zeros(Rain.shape)
        RefrozenRain[np.where(Rain >= SI_out)] = SI_out[np.where(Rain >= SI_out)] # all potential superimposed ice is used up in refreezing, and then some additional rain happens and runs off
        RefrozenRain[np.where(Rain < SI_out)] = Rain[np.where(Rain < SI_out)] # all rain is refrozen, some potentialSI is leftover
        RefrozenRain[np.where(np.isnan(Sfc))] = np.nan
        SI_out -= RefrozenRain # Subtract the amount of rain that is refrozen from the potentialSI

        # Update output arrays for this timestep:
        # Total Melt = Snow Melt + Ice Melt
        # Net ablation = total melt - refreezing
        # Net balance = Accumulation - net ablation
        # ===================================================================== 
        IceMelt[timestamp,:,:] = Mice
        SnowMelt[timestamp,:,:] = Msnow
        RefrozenMelt[timestamp,:,:] = (Refreezing + RefrozenRain)
        MassBal[timestamp,:,:] = New_snowfall - ((Msnow - RefrozenMelt[timestamp,:,:]) + Mice)  

        # Update snowpack and CC trackers for next timestep:
        # ===================================================================== 
        Snowpack_tracker[timestamp+1,:,:] = SP_out
        PotentialSI_tracker[timestamp+1,:,:] = SI_out   
        
    # Save outputs for the year before starting next year:
    # =========================================================================
    print('Saving model outputs for',year)
    sys.stdout.flush()
    if SaveMBonly == False:
        save_to_netcdf(IceMelt, 'Ice melt', os.path.join(OUTPUT_PATH,'Icemelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 
        save_to_netcdf(SnowMelt, 'Snow melt', os.path.join(OUTPUT_PATH,'Snowmelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 
        save_to_netcdf(RefrozenMelt, 'Refreezing', os.path.join(OUTPUT_PATH,'Refreezing_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 
        save_to_netcdf(Snowpack_tracker[:-1], 'Snowdepth', os.path.join(OUTPUT_PATH,'Snowdepth_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 
        save_to_netcdf(PotentialSI_tracker[:-1], 'Ptau', os.path.join(OUTPUT_PATH,'Ptau_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid)       
    else:
        pass
    
    save_to_netcdf(MassBal, 'Net balance', os.path.join(OUTPUT_PATH,'Netbalance_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid) 

    inT.close()
    inP.close()
    inS.close()

print('RUN COMPLETE!')