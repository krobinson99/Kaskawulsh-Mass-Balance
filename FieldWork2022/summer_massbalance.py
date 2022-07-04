# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:08:22 2022

the purpose of this script is to plot the ablation/accumulation fields for July 15 - Aug 31
(average 2007-2018) to see what we might expect in the field this summer

@author: katierobinson
"""

import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset
#import h5py
import os
import sys
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import KWtributaries

# use baseline case: file:///F:/Mass Balance Model/OUTPUTS/Kaskonly_Baseline

#INITIALIZE SCRIPT#---------------------------------------------------

# inputs
#where should the script get the model outputs from:
FILEPATH = 'F:/Mass Balance Model/OUTPUTS/Kaskonly_Baseline'
File_glacier_in = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/kaskonly.txt'
Case = 'Baseline'
# where should the outputs of this script be stored
OUTPUT_PATH =  'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022'
#NARR_PATH = 'F:\Mass Balance Model\BiasCorrectedInputs\Kaskonly_R2S=1'
#'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs'
domain = 'kaskonly'

summermelt_gen = False #need to re-generate the summersnow.npy / summermelt.npy files? if no: False
summersnow_gen = False

# times
start_year = 2006
end_year = 2018

# params
Params_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/final_params_deb.csv'
Catchment = False

years = []
year = start_year
while year <= end_year:
  years.append(year)
  year = year + 1

##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
#if debris is TRUE then we are working from the kaskonly_deb.txt file,
#if debris is FALSE then we are working from the kaskonly.txt file:
if Catchment == True:
    Ix = glacier[:,4] 
    Iy = glacier[:,5] 
    Ih = glacier[:,6] 
else:
    Ix = glacier[:,3]
    Iy = glacier[:,4] 
    Ih = glacier[:,5]         
        
Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)

# ----------------------------------------------------------------------------
# calculated locations where the value is NaN:
# nanlocs is an array of True/False booleans with shape (218,328)
#File_temp_in = os.path.join(NARR_PATH,'Tempkaskonly_BC_2007.nc')
#File_temp_in = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\Tempkaskonly2007.nc'
File_temp_in = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC\\Tempkaskonly2011.nc'

inT = Dataset(File_temp_in, "r")
T_var = 'Temperature'
T_array = inT.variables[T_var][:]
nanlocs = np.isnan(T_array[0,:,:])


#---------------------------------------------------------------------------
# Get param combination 
params = np.loadtxt(Params_file)
aice = params[0,:]
asnow = params[1,:]
MF = params[2,:]

if summermelt_gen == True:
    maxmelt_list = []
    allruns_Melt = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    #for sim in range(0, len(aice)):
    for year in years[1:]:
        #print('starting sim ' + str(sim))
        print('starting year ' + str(year))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        #for year in years[1:]: # start at index 1 (first year skipped for spin-up)
        for sim in range(0, len(aice)):
            print(sim)
            
            MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Melt'
            MB_array = MB_in.variables[MB_var][:]
            
            #what timesteps correspond to July 15 and Aug 31
            #formula = 8*(DOY-1)
            #july 15 DOY = 196
            # aug 31 DOY = 243
            
            jul15 = 8*(196-1)
            aug31 = 8*(243-1)
        
            N_MB = np.nansum(MB_array[jul15:aug31,:,:], axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            #N_MB[nanlocs] = np.nan
            
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            MB_empty_container[sim,:,:] = N_MB
        
        maxmelt = np.nanmax(MB_empty_container)
        maxmelt_list.append(maxmelt)
        
        #out_file_MB = 'Net_Melt' + str(sim) + '.npy'
        #MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        #np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    #------------------------------------------------------------------------------
    # Concantenate all simulations and take the average mass balance
    
        #run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(MB_empty_container, axis = 0)
        
        run_mean_MB[nanlocs] = np.nan
        
        
        allruns_Melt[np.where(np.asarray(years) == year-1),:] = run_mean_MB
        
    summermelt = np.nanmean(allruns_Melt, axis = 0)
    
    #save output of run means
    out_run_means = 'summermelt_' + Case + '.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, summermelt)
    
    out_run_means = 'allyears_summermelt_' + Case + '.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_Melt)

else:
    pass

summermelt = np.load(os.path.join(OUTPUT_PATH,'summermelt_' + Case + '.npy'))

allyearsmelt = np.load(os.path.join(OUTPUT_PATH,'allyears_summermelt_' + Case + '.npy'))
maxmelt = np.nanmax(allyearsmelt)

for i in range(0,len(allyearsmelt)):
    print(i+2007)
    loc = np.where(allyearsmelt[i] == np.nanmax(allyearsmelt[i]))
    print(np.nanmax(allyearsmelt[i]))
    #print(Zgrid[loc])
    
#max ablation in the elevation zone 1750-1850 m asl
targetelev = np.where((Zgrid > 1750) & (Zgrid < 1850))

for i in range(0,len(allyearsmelt)):
    print(i+2007)
    yearmelt = allyearsmelt[i]
    print(np.nanmax(yearmelt[targetelev]))



#search for minimum ablation in the trunk area only
tribarray = KWtributaries('file:///F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KWTributaries.csv',Zgrid)
nanlocs = np.where(np.isnan(Zgrid))
tribarray[nanlocs] = np.nan

nontrunk = np.where(tribarray < 5)
for i in range(0,len(allruns_Melt)):
    nannontrunk = allruns_Melt[i]
    nannontrunk[nontrunk] = np.nan
    allruns_Melt[i] = nannontrunk
    

minmelt = np.where(allruns_Melt == np.nanmin(allruns_Melt))

if summersnow_gen == True:
    allruns_Snow = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for sim in range(0, len(aice)):
        print('starting sim ' + str(sim))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        for year in years[1:]: # start at index 1 (first year skipped for spin-up)
            print(year)
            
            MB_file_name = 'Accumulation' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Accumulation'
            MB_array = MB_in.variables[MB_var][:]
            
            #what timesteps correspond to July 15 and Aug 31
            #formula = 8*(DOY-1)
            #july 15 DOY = 196
            # aug 31 DOY = 243
            
            jul15 = 8*(196-1)
            aug31 = 8*(243-1)
        
            N_MB = np.nansum(MB_array[jul15:aug31,:,:], axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            #N_MB[nanlocs] = np.nan
            
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            MB_empty_container[np.where(np.asarray(years) == year-1),:] = N_MB
        
        #out_file_MB = 'Net_Melt' + str(sim) + '.npy'
        #MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        #np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    #------------------------------------------------------------------------------
    # Concantenate all simulations and take the average mass balance
    
        #run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(MB_empty_container, axis = 0)
        
        run_mean_MB[nanlocs] = np.nan
        
        allruns_Snow[sim,:,:] = run_mean_MB
        
    summersnow = np.nanmean(allruns_Snow, axis = 0)
        
    #save output of run means
    out_run_means = 'summersnow_' + Case + '.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, summersnow)
else:
    pass

summersnow = np.load(os.path.join(OUTPUT_PATH,'summersnow_' + Case + '.npy'))

###############################################################################
#PLOT OUTPUTS
plt.figure(figsize=(12,7))
plt.contourf(np.flipud(summermelt), cmap = 'YlOrRd', levels = np.linspace(0,5, 21))
legend = plt.colorbar()
legend.ax.set_ylabel('Melt (m w.e.)', rotation=270)
#legend.ax.set_ylabel('Precipitation (m w.e./yr)', rotation=270,x=1.1)
plt.title('Summer (Jul 15 - Aug 31) Ablation (2007-2008)')
#plt.savefig('summmer_melt.png')

plt.figure(figsize=(12,7))
plt.contourf(np.flipud(summersnow), cmap = 'PuRd', levels = np.linspace(0,0.5, 21))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270)
#legend.ax.set_ylabel('Precipitation (m w.e./yr)', rotation=270,x=1.1)
plt.title('Summer (Jul 15 - Aug 31) Accumulation (2007-2008)')
#plt.savefig('summmer_snow.png')


plt.figure(figsize=(12,7))
plt.contourf(np.flipud(allruns_Melt[4]), cmap = 'YlOrRd', levels = np.linspace(0,5, 21))
legend = plt.colorbar()
legend.ax.set_ylabel('Melt (m w.e.)', rotation=270)
#legend.ax.set_ylabel('Precipitation (m w.e./yr)', rotation=270,x=1.1)
plt.title('Summer (Jul 15 - Aug 31) Ablation (2007-2008)')