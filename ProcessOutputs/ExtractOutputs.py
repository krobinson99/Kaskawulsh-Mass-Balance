# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 08:53:53 2021

The purpose of this script is to extract outputs from the MBM runs and plot
relevant metrics

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

#---------------------------------------------------------------------
# WHAT DO WE WANT THE SCRIPT TO DO:
Calculate_Mean_MB = False # calculate mean mass balance across all runs: only run this section if you don't already have an 'allrunMBs_Case' file
Calculate_DailyMB = False

Calculate_NetAblation = False
Calculate_DailyNetAblation = False

Calculate_Accumulation = False
Calculate_DailyAccumulation = False
Calculate_Monthly_Distributed_Snow = False

Calculate_Refreezing = False
Calculate_DailyRefreezing = False

Calculate_Rain = False
Calculate_DailyRain = False
Calculate_Monthly_Distributed_Rain = False

Calculate_DailyIceMelt = True
Calculate_Distributed_IceMelt = True
Calculate_DailySnowMelt = True
Calculate_Distributed_SnowMelt = True

Table_of_MBvalues = True
Print_all_runmeans = False

Calculate_distributed_temp = True
Calculate_daily_temp = True
Calculate_Monthly_Distributed_Temp = True

Calculate_distributed_NARR_precip = True
Calculate_daily_NARR_precip = False
Calculate_daily_BC_Rain = False
Calculate_monthly_distributed_NARR_rain = False
Calculate_monthly_distributed_NARR_precip = False
Calculate_distributed_NARR_NoBCRain = False
Calculate_distributed_NARR_precip_withBC = False

#INITIALIZE SCRIPT#---------------------------------------------------

# inputs
#where should the script get the model outputs from:
FILEPATH = 'D:\Model Runs\Debris Tests\VariableThicknessTest_2007-2018'
File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kaskonly.txt'
# where should the outputs of this script be stored
OUTPUT_PATH =  'D:\Model Runs\Debris Tests\VariableThicknessTest_2007-2018\ProcessedOutputs'
NARR_PATH = 'F:\Mass Balance Model\BiasCorrectedInputs\Kaskonly_R2S=1'
#'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs'
domain = 'kaskonly'

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
    Ih = glacier[:,2]        
        
Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)

# ----------------------------------------------------------------------------
# calculated locations where the value is NaN:
# nanlocs is an array of True/False booleans with shape (218,328)
#File_temp_in = os.path.join(NARR_PATH,'Tempkaskonly_BC_2007.nc')
#File_temp_in = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\Tempkaskonly2007.nc'
File_temp_in = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC\\Tempkaskonly2011.nc'
#File_temp_in = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Catchment_R2S=1\\Tempkaskonly_BC_2007.nc'
inT = Dataset(File_temp_in, "r")
T_var = 'Temperature'
T_array = inT.variables[T_var][:]
nanlocs = np.isnan(T_array[0,:,:])
#T_array_shifted = T_array - 0.9
#T_array_shifted_distributed = np.mean(T_array_shifted, axis=0)
#plt.contourf(np.flipud(np.mean(T_array_shifted, axis = 0)), cmap = 'RdYlBu', levels = np.linspace(-12,7, 20))
#T_array_file = 'Temperature_shifted-0.9.npy'
#Tarray = os.path.join(OUTPUT_PATH,T_array_file)
#np.save(Tarray, T_array_shifted_distributed)

#---------------------------------------------------------------------------
# Get param combination 
params = np.loadtxt(Params_file)
aice = params[0,:]
asnow = params[1,:]
MF = params[2,:]

#just adding this to make the param list shorter -- delete when done. 
#aice = aice[:23]
#asnow = asnow[:23]
#MF = MF[:23]


# Creat Zgrid as a .npy array
out_file_Zgrid = 'Elevation_grid.npy'
Zgrid_out_path = os.path.join(OUTPUT_PATH,out_file_Zgrid)
np.save(Zgrid_out_path, Zgrid) #results in a file caled Net_MB1.npy for each simulation

out_file_Xgrid = 'X_grid.npy'
Xgrid_out_path = os.path.join(OUTPUT_PATH,out_file_Xgrid)
np.save(Xgrid_out_path, Xgrid) #results in a file caled Net_MB1.npy for each simulation

out_file_Ygrid = 'y_grid.npy'
Ygrid_out_path = os.path.join(OUTPUT_PATH,out_file_Ygrid)
np.save(Ygrid_out_path, Ygrid) #results in a file caled Net_MB1.npy for each simulation

# ---------------------------------------------------------------------------
# Concatenate all years from each model run  

if Calculate_Mean_MB == True:
    print('Extracting Mass Balance:')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
    for sim in range(0, len(aice)):
        print('starting sim ' + str(sim))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        for year in years[1:]: # start at index 1 (first year skipped for spin-up)
            print(year)
            
            MB_file_name = 'MB' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'MB'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            
            # insert the net mass balance for the year into the empty MB container for that year
            MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
        
        out_file_MB = 'Net_MB' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    # Concantenate all simulations and take the average mass balance
    
        run_MB = np.load(MB_out_path)
        
        run_mean_MB = np.nanmean(run_MB, axis = 0)
        run_mean_MB[nanlocs] = np.nan
        
        allruns_MB[sim,:,:] = run_mean_MB
        
    #save output of run means
    out_run_means = 'allrunMBs.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_MB)
    print('Mean MB calculated for all runs: Outputs are saved in: '+ str(out_allruns))
else:
    print('Calculate mean mass balances SKIPPED')

if Print_all_runmeans == True:
    
    allyearmeans = []
    for sim in range(0,len(aice)):
        out_file_MB = 'Net_MB' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        runMB = np.load(MB_out_path)
        #print(runMB.shape)
        meanMBx = np.nanmean(runMB,axis=1)
        #print(meanMBx.shape)
        meanMBy = np.nanmean(meanMBx,axis = 1)
        #print(meanMBy)
        allyearmeans.append(np.nanmean(meanMBy))
        
    

#-------------------------------------------------------------------------------

if Calculate_NetAblation == True:
    print('Extracting Net Melt:')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for sim in range(0, len(aice)):
        print('starting sim ' + str(sim))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        for year in years[1:]: # start at index 1 (first year skipped for spin-up)
            print(year)
            
            MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Melt'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            #Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            #Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            #Daysum_array = np.nanmean(Daysum_array, axis = 1)
            #Daysum_array = np.nanmean(Daysum_array, axis = 1) #mean daily net melt across whole domain
            
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            MB_empty_container[np.where(np.asarray(years) == year-1),:] = N_MB
        
        out_file_MB = 'Net_Melt' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    #------------------------------------------------------------------------------
    # Concantenate all simulations and take the average mass balance
    
        run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(run_MB, axis = 0)
        
        run_mean_MB[nanlocs] = np.nan
        
        allruns_MB[sim,:,:] = run_mean_MB
        
    #save output of run means
    out_run_means = 'allrunNetMelts.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_MB)
    print('Net Melt calculated for all runs: Outputs are saved in: '+ str(out_allruns))
else:
    print('Calculate Net Melt SKIPPED')
    
#----------------------------------------------------------------------------
    
if Calculate_Accumulation == True:
    print('Extracting Accumulation:')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
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
            #print('variable in')
            
            N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            
            # insert the net mass balance for the year into the empty MB container for that year
            MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
        
        out_file_MB = 'Net_Accumulation' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    #------------------------------------------------------------------------------
    # Concantenate all simulations and take the average mass balance
    
        run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(run_MB, axis = 0) #calculates mean across one simulations for all time
        run_mean_MB[nanlocs] = np.nan
        
        allruns_MB[sim,:,:] = run_mean_MB
        
    #save output of run means
    out_run_means = 'allrunAccumulations.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_MB)
    print('Accumulation calculated for all runs: Outputs are saved in: '+ str(out_allruns))
else:
    print('Calculate Accumulation SKIPPED')
    
#-------------------------------------------------------------------------------
   
if Calculate_Rain == True:
    print('Extracting Rain:')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for sim in range(0, len(aice)):
        print('starting sim ' + str(sim))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        for year in years[1:]: # start at index 1 (first year skipped for spin-up)
            print(year)
            
            MB_file_name = 'Rain' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Rain'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            
            # insert the net mass balance for the year into the empty MB container for that year
            MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
        
        out_file_MB = 'Net_Rain' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    

    # Concantenate all simulations and take the average mass balance
    
        run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(run_MB, axis = 0) #calculates mean across one simulations for all time
        run_mean_MB[nanlocs] = np.nan
        
        allruns_MB[sim,:,:] = run_mean_MB
        
    #save output of run means
    out_run_means = 'allrunRains.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_MB)
    print('Rain calculated for all runs: Outputs are saved in: '+ str(out_allruns))
else:
    print('Calculate Rain SKIPPED')


#--------------------------------------------------------------------------------
   
if Calculate_Refreezing == True:
    print('Extracting Distributed Refreezing:')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for sim in range(0, len(aice)):
        print('starting sim ' + str(sim))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        for year in years[1:]: # start at index 1 (first year skipped for spin-up)
            print(year)
            
            MB_file_name = 'Refreezing' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Refreezing'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            
            # insert the net mass balance for the year into the empty MB container for that year
            MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
        
        out_file_MB = 'Refreezing' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    # Concantenate all simulations and take the average mass balance
    
        run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(run_MB, axis = 0) #calculates mean across one simulations for all time
        run_mean_MB[nanlocs] = np.nan
        
        allruns_MB[sim,:,:] = run_mean_MB
        
    #save output of run means
    out_run_means = 'allrunRefreezings.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_MB)
    print('Refreezing calculated for all runs: Outputs are saved in: '+ str(out_allruns))
else:
    print('Calculate Refreezing SKIPPED')
    
    
#---------------------------------------------------------------------------------
    
if Calculate_DailyNetAblation == True:
    print('Extracting Daily Net Melt by YEAR:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
            
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Melt'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            Daysum_meanx = np.nanmean(Daysum_array,axis=1)
            Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            Yearly_NetMelt_container[sim,:] = Daysum_meany
        
        out_file_MB = 'Daily_Net_Melt' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Yearly_NetMelt_container) #results in a file caled Net_MB1.npy for each simulation
    
    print('Daily Net Melt calculated for all years')
else:
    print('Calculate Daily Net Melt SKIPPED')
    
#------------------------------------------------------------------------------    

if Calculate_DailyMB == True:
    print('Extracting Daily Mass Balance by YEAR:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'MB' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'MB'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            Year_Sum = np.nansum(MB_array,axis=0)
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            Daysum_meanx = np.nanmean(Daysum_array,axis=1)
            Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
            #Daysum_array = np.nanmean(Daysum_array, axis = 1)
            #Daysum_array = np.nanmean(Daysum_array, axis = 1) #mean daily net melt across whole domain
            
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            Yearly_NetMelt_container[sim,:] = Daysum_meany
        
        out_file_MB = 'Daily_MB' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Yearly_NetMelt_container) #results in a file caled Net_MB1.npy for each simulation
    
    print('Daily Mass Balance calculated for all years')
else:
    print('Calculate Daily Mass Balance SKIPPED')

#------------------------------------------------------------------------------
if Calculate_DailyAccumulation == True:
    print('Extracting Daily Accumulation by YEAR:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'Accumulation' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Accumulation'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
           
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            Daysum_meanx = np.nanmean(Daysum_array,axis=1)
            Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            Yearly_NetMelt_container[sim,:] = Daysum_meany
        
        out_file_MB = 'Daily_Accumulation' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Yearly_NetMelt_container) #results in a file caled Net_MB1.npy for each simulation
    
    print('Daily Accumulation calculated for all years')
else:
    print('Calculate Daily Accumulation SKIPPED')
    
#------------------------------------------------------------------------------
if Calculate_DailyRain == True:
    print('Extracting Daily Rain by YEAR:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'Rain' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Rain'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            Daysum_meanx = np.nanmean(Daysum_array,axis=1)
            Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            Yearly_NetMelt_container[sim,:] = Daysum_meany
        
        out_file_MB = 'Daily_Rain' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Yearly_NetMelt_container) #results in a file caled Net_MB1.npy for each simulation
    
    print('Daily Rain calculated for all years')
else:
    print('Calculate Daily Rain SKIPPED')


#-------------------------------------------------------------------------------
if Calculate_DailyRefreezing == True:
    print('Extracting Daily Refreezing by YEAR:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            MB_file_name = 'Refreezing' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Refreezing'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
           
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            Daysum_meanx = np.nanmean(Daysum_array,axis=1)
            Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            Yearly_NetMelt_container[sim,:] = Daysum_meany
        
        out_file_MB = 'Refreezing_' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Yearly_NetMelt_container) #results in a file caled Net_MB1.npy for each simulation
    
    print('Daily Refreezing calculated for all years')
else:
    print('Calculate Refreezing SKIPPED')
    

if Table_of_MBvalues == True:
    print('Calculating table of MB values')
    
    # Set up empty array to store all mean MB values in:
    MB_values = np.zeros((len(aice) + 1,len(aice) + 1))
    year = 2007
    sim = int(1.0)
    for i in range(1,len(MB_values)):
        MB_values[i][0] = year
        MB_values[0][i] = sim
        sim += 1
        year += 1
    
    def MB_byYearandSim(year,sim):
        # Import the model output file for a specific year and simulation
        MB_file_name = 'MB' + str(year) + str(sim) + '.nc'
        MB_File_in = os.path.join(FILEPATH,MB_file_name)
        MB_in = Dataset(MB_File_in, "r")
        MB_var = 'MB'
        MB_array = MB_in.variables[MB_var][:]
        
        # calculate the MEAN MASS BALANCE over all time and space in this file:
        Total_MB = np.nansum(MB_array, axis =0) #sum over all time
        Mean_mb = np.nanmean(Total_MB, axis = 0) #average over all x
        Mean_MB = np.nanmean(Mean_mb , axis = 0) # average over all y
        
        
        return Mean_MB
    
    i = 1
    for year in years[1:]:
        for sim in range(0, len(aice)):
            print('year: ' + str(year) + '...' + ' sim: ' + str(sim))
            #print(MB_byYearandSim(year,sim))
            MB_values[year - 2006][sim+1] = MB_byYearandSim(year,sim)
    
    np.set_printoptions(suppress=True)
    MB_Table = np.round(MB_values,2)
    
    Table_filename = 'MB_Table.npy' 
    Table_filepath = os.path.join(OUTPUT_PATH,Table_filename)
    np.save(Table_filepath,MB_Table)

    # calculate mean by by row (ie. by year)
    year_means = []
    years = []
    sim_means = []
    for i in range(1,len(aice)+1):
        yr = MB_Table[i][0]
        years.append(yr)
        ymean = round(np.mean(MB_Table[i][1:]),2)
        year_means.append(ymean)
        
        smean = round(np.mean(MB_Table[1:,i]),2)
        sim_means.append(smean)
else:
    pass

#---------------------------------------------------------------------------------
    
if Calculate_DailyIceMelt == True:
    print('Extracting Daily Ice Melt by YEAR:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'IceMelt' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'IceMelt'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            Daysum_meanx = np.nanmean(Daysum_array,axis=1)
            Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            Yearly_NetMelt_container[sim,:] = Daysum_meany
        
        out_file_MB = 'Daily_IceMelt' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Yearly_NetMelt_container) 
    
    print('Daily Ice Melt calculated for all years')
else:
    print('Calculate Daily Ice Melt SKIPPED')
    
#------------------------------------------------------------------------------    
    
if Calculate_DailySnowMelt == True:
    print('Extracting Daily Snow Melt by YEAR:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'SnowMelt' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'SnowMelt'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            Daysum_meanx = np.nanmean(Daysum_array,axis=1)
            Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            Yearly_NetMelt_container[sim,:] = Daysum_meany
        
        out_file_MB = 'Daily_SnowMelt' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Yearly_NetMelt_container) #results in a file caled Net_MB1.npy for each simulation
    
    print('Daily Snow Melt calculated for all years')
else:
    print('Calculate Daily Snow Melt SKIPPED')    

#---------------------------------------------------------------------------------
if Calculate_Monthly_Distributed_Rain == True:
    print('Extracting Monthly distributed rain:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation

    allyears_jan_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_feb_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_mar_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_april_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_may_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_jun_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_july_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_aug_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_sept_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_oct_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_nov_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_dec_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    for year in years[1:]:
        #Monthly_container = np.empty(len(aice),12, len(T_array[1]),len(T_array[1,1,:])) #num of sims, num of months, x, y
        print(year)
        
        sims_jan_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_feb_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_mar_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_april_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_may_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_june_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_july_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_aug_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_sept_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_oct_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_nov_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_dec_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'Rain' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Rain'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            
            jan = Daysum_array[0:31,:,:]
            feb = Daysum_array[31:59,:,:]
            march = Daysum_array[59:90,:,:]
            april = Daysum_array[90:120,:,:]
            may = Daysum_array[120:151,:,:]
            june = Daysum_array[151:181,:,:]
            july = Daysum_array[181:212,:,:]
            aug = Daysum_array[212:243,:,:]
            sept = Daysum_array[243:273,:,:]
            octo = Daysum_array[273:304,:,:]
            nov = Daysum_array[304:334,:,:]
            dec = Daysum_array[334:365,:,:]   
            
            jan_xy = np.nanmean(jan,axis=0)
            feb_xy = np.nanmean(feb,axis=0)
            march_xy = np.nanmean(march,axis=0)
            april_xy = np.nanmean(april,axis=0)
            may_xy = np.nanmean(may,axis=0)
            june_xy = np.nanmean(june,axis=0)
            july_xy = np.nanmean(july,axis=0)
            aug_xy = np.nanmean(aug,axis=0)
            sept_xy = np.nanmean(sept,axis=0)
            oct_xy = np.nanmean(octo,axis=0)
            nov_xy = np.nanmean(nov,axis=0)
            dec_xy = np.nanmean(dec,axis=0)
            
            sims_jan_container[sim,:,:] = jan_xy
            sims_feb_container[sim,:,:] = feb_xy
            sims_mar_container[sim,:,:] = march_xy
            sims_april_container[sim,:,:] = april_xy
            sims_may_container[sim,:,:] = may_xy
            sims_june_container[sim,:,:] = june_xy
            sims_july_container[sim,:,:] = july_xy
            sims_aug_container[sim,:,:] = aug_xy
            sims_sept_container[sim,:,:] = sept_xy
            sims_oct_container[sim,:,:] = oct_xy
            sims_nov_container[sim,:,:] = nov_xy
            sims_dec_container[sim,:,:] = dec_xy
            
        allsims_avg_jan = np.nanmean(sims_jan_container,axis=0)
        allsims_avg_feb = np.nanmean(sims_feb_container,axis=0)
        allsims_avg_mar = np.nanmean(sims_mar_container,axis=0)
        allsims_avg_april = np.nanmean(sims_april_container,axis=0)
        allsims_avg_may = np.nanmean(sims_may_container,axis=0)
        allsims_avg_jun = np.nanmean(sims_june_container,axis=0)
        allsims_avg_july = np.nanmean(sims_july_container,axis=0)
        allsims_avg_aug = np.nanmean(sims_aug_container,axis=0)
        allsims_avg_sept = np.nanmean(sims_sept_container,axis=0)
        allsims_avg_oct = np.nanmean(sims_oct_container,axis=0)
        allsims_avg_nov = np.nanmean(sims_nov_container,axis=0)
        allsims_avg_dec = np.nanmean(sims_dec_container,axis=0)
        
        allyears_jan_container[year-2007,:,:] = allsims_avg_jan
        allyears_feb_container[year-2007,:,:] = allsims_avg_feb
        allyears_mar_container[year-2007,:,:] = allsims_avg_mar
        allyears_april_container[year-2007,:,:] = allsims_avg_april
        allyears_may_container[year-2007,:,:] = allsims_avg_may
        allyears_jun_container[year-2007,:,:] = allsims_avg_jun
        allyears_july_container[year-2007,:,:] = allsims_avg_july
        allyears_aug_container[year-2007,:,:] = allsims_avg_aug
        allyears_sept_container[year-2007,:,:] = allsims_avg_sept
        allyears_oct_container[year-2007,:,:] = allsims_avg_oct
        allyears_nov_container[year-2007,:,:] = allsims_avg_nov
        allyears_dec_container[year-2007,:,:] = allsims_avg_dec
    
    avg_jan = np.nanmean(allyears_jan_container,axis=0)
    avg_feb = np.nanmean(allyears_feb_container,axis=0)
    avg_mar = np.nanmean(allyears_mar_container,axis=0)
    avg_april = np.nanmean(allyears_april_container,axis=0)
    avg_may = np.nanmean(allyears_may_container,axis=0)
    avg_jun = np.nanmean(allyears_jun_container,axis=0)
    avg_july = np.nanmean(allyears_july_container,axis=0)
    avg_aug = np.nanmean(allyears_aug_container,axis=0)
    avg_sept = np.nanmean(allyears_sept_container,axis=0)
    avg_oct = np.nanmean(allyears_oct_container,axis=0)
    avg_nov = np.nanmean(allyears_nov_container,axis=0)
    avg_dec = np.nanmean(allyears_dec_container,axis=0)
    
    jan_outpath = os.path.join(OUTPUT_PATH,'jan_distributedRain.npy')
    feb_outpath = os.path.join(OUTPUT_PATH,'feb_distributedRain.npy')
    march_outpath = os.path.join(OUTPUT_PATH,'march_distributedRain.npy')
    april_outpath = os.path.join(OUTPUT_PATH,'april_distributedRain.npy')
    may_outpath = os.path.join(OUTPUT_PATH,'may_distributedRain.npy')
    june_outpath = os.path.join(OUTPUT_PATH,'june_distributedRain.npy')
    july_outpath = os.path.join(OUTPUT_PATH,'july_distributedRain.npy')
    aug_outpath = os.path.join(OUTPUT_PATH,'aug_distributedRain.npy')
    sept_outpath = os.path.join(OUTPUT_PATH,'sept_distributedRain.npy')
    oct_outpath = os.path.join(OUTPUT_PATH,'oct_distributedRain.npy')
    nov_outpath = os.path.join(OUTPUT_PATH,'nov_distributedRain.npy')
    dec_outpath = os.path.join(OUTPUT_PATH,'dec_distributedRain.npy')
    
    np.save(jan_outpath,avg_jan)
    np.save(feb_outpath,avg_feb)
    np.save(march_outpath,avg_mar)
    np.save(april_outpath,avg_april)
    np.save(may_outpath,avg_may)
    np.save(june_outpath,avg_jun)
    np.save(july_outpath,avg_july)
    np.save(aug_outpath,avg_aug)
    np.save(sept_outpath,avg_sept)
    np.save(oct_outpath,avg_oct)
    np.save(nov_outpath,avg_nov)
    np.save(dec_outpath,avg_dec)
    
    print('distributed rain calculated for all months')
else:
    print('Calculate Distributed Monthly Rain SKIPPED')
    
#---------------------------------------------------------------------------------------------------------------
if Calculate_distributed_temp == True:
    print('extracting distributed temp')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
    MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
    for year in years[1:]: # start at index 1 (first year skipped for spin-up)
        print(year)
        
        MB_file_name = 'Temp' + str(domain) + str(year) + '.nc'
        MB_File_in = os.path.join(NARR_PATH,MB_file_name)

        MB_in = Dataset(MB_File_in, "r")
        #print('file in')
        MB_var = 'Temperature'
        T_array = MB_in.variables[MB_var][:]
        
        # apply temp shift
        MB_array = T_array
        
        N_MB = np.nanmean(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        
        # insert the net mass balance for the year into the empty MB container for that year
        MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.

    Temp_allyears_mean = np.nanmean(MB_empty_container, axis =0)
    out_file_MB = 'Mean_Distributed_Temp.npy'
    MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
    np.save(MB_out_path, Temp_allyears_mean) #results in a file caled Net_MB1.npy for each simulation

else:
    pass
#------------------------------------------------------------------------------
if Calculate_daily_temp == True:
    print('calculating daily temperatures')
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        Temp_file_name = 'Temp' + str(domain) + '_BC_' + str(year) + '.nc'
        Temp_file_in = os.path.join(NARR_PATH,Temp_file_name)
        
        MB_in = Dataset(Temp_file_in, "r")
        #print('file in')
        MB_var = 'Temperature'
        T_array = MB_in.variables[MB_var][:]
        
        # apply tempshift 
        MB_array = T_array
        
        Day_Sum = [np.nanmean(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
        Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
        for day in range(0,len(Daysum_array)):
            Daysum_array[day,:,:][nanlocs] = np.nan
                
        Daysum_meanx = np.nanmean(Daysum_array,axis=1)
        Daysum_meany = np.nanmean(Daysum_meanx,axis=1)
        
        #save daily distributed temp
        out_file_dailydistributedT = 'DailyDistributedTemp' + str(year) + '.npy'
        File_out_path = os.path.join(OUTPUT_PATH,out_file_dailydistributedT) 
        np.save(File_out_path, Daysum_array)
        
        #save daily temp timeseries
        out_file_dailyTtimeseries = 'DailyTimeseriesTemp' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_dailyTtimeseries)  
        np.save(MB_out_path, Daysum_meany) 
    print('distributed temp saved')
else:
    pass

#-----------------------------------------------------------------------------
if Calculate_Distributed_IceMelt == True:
    print('Extracting Ice Melt:')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for sim in range(0, len(aice)):
        print('starting sim ' + str(sim))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        for year in years[1:]: # start at index 1 (first year skipped for spin-up)
            print(year)
            
            MB_file_name = 'IceMelt' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'IceMelt'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            #Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            #Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            #Daysum_array = np.nanmean(Daysum_array, axis = 1)
            #Daysum_array = np.nanmean(Daysum_array, axis = 1) #mean daily net melt across whole domain
            
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            MB_empty_container[np.where(np.asarray(years) == year-1),:] = N_MB
        
        out_file_MB = 'IceMelt' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    #------------------------------------------------------------------------------
    # Concantenate all simulations and take the average mass balance
    
        run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(run_MB, axis = 0)
        
        run_mean_MB[nanlocs] = np.nan
        
        allruns_MB[sim,:,:] = run_mean_MB
        
    #save output of run means
    out_run_means = 'allrunIceMelts.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_MB)
    print('Ice Melt calculated for all runs: Outputs are saved in: '+ str(out_allruns))
else:
    print('Calculate Ice Melt SKIPPED')

#-------------------------------------------------------------------------------

if Calculate_Distributed_SnowMelt == True:
    print('Extracting Snow Melt:')
    allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation
    for sim in range(0, len(aice)):
        print('starting sim ' + str(sim))
    
        #create empty array with shape [t,x,y] = number of years, MB_array dims (ie. (12,218,328) for kaskonly)
        MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
        
        for year in years[1:]: # start at index 1 (first year skipped for spin-up)
            print(year)
            
            MB_file_name = 'SnowMelt' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'SnowMelt'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            #Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            #Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            #Daysum_array = np.nanmean(Daysum_array, axis = 1)
            #Daysum_array = np.nanmean(Daysum_array, axis = 1) #mean daily net melt across whole domain
            
            # insert the net mass balance for the year into the empty MB container for that year
            #MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.
            MB_empty_container[np.where(np.asarray(years) == year-1),:] = N_MB
        
        out_file_MB = 'SnowMelt' + str(sim) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, MB_empty_container) #results in a file caled Net_MB1.npy for each simulation
    
    #------------------------------------------------------------------------------
    # Concantenate all simulations and take the average mass balance
    
        run_MB = np.load(MB_out_path) # loads the file w each full run in it
        
        run_mean_MB = np.nanmean(run_MB, axis = 0)
        
        run_mean_MB[nanlocs] = np.nan
        
        allruns_MB[sim,:,:] = run_mean_MB
        
    #save output of run means
    out_run_means = 'allrunSnowMelts.npy'
    out_allruns = os.path.join(OUTPUT_PATH,out_run_means)
    np.save(out_allruns, allruns_MB)
    print('Snow Melt calculated for all runs: Outputs are saved in: '+ str(out_allruns))
else:
    print('Calculate Snow Melt SKIPPED')
    
#------------------------------------------------------------------------------
if Calculate_distributed_NARR_precip == True:
    print('extracting distributed precip')

    MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
    for year in years[1:]: # start at index 1 (first year skipped for spin-up)
        print(year)
        
        MB_file_name = 'netSnow' + str(domain) + str(year) + '.nc'
        MB_File_in = os.path.join(NARR_PATH,MB_file_name)

        MB_in = Dataset(MB_File_in, "r")
        #print('file in')
        MB_var = 'Temperature'
        MB_array = MB_in.variables[MB_var][:]
        
        
        N_MB = np.sum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        
        # insert the net mass balance for the year into the empty MB container for that year
        MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.

    Temp_allyears_mean = np.nanmean(MB_empty_container, axis =0)
    out_file_MB = 'allrunPrecipitation_' + 'DownscaledNoBC' + '.npy'
    MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
    np.save(MB_out_path, Temp_allyears_mean) #results in a file caled Net_MB1.npy for each simulation

else:
    pass
#------------------------------------------------------------------------------
if Calculate_daily_NARR_precip == True:
    print('calculating daily temperatures')
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        Temp_file_name = 'netSnow' + str(domain) + str(year) + '.nc'
        Temp_file_in = os.path.join(NARR_PATH,Temp_file_name)
        
        MB_in = Dataset(Temp_file_in, "r")
        #print('file in')
        MB_var = 'Net snow'
        T_array = MB_in.variables[MB_var][:]
        
        # apply tempshift 
        MB_array = T_array
        
        Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
        Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
        for day in range(0,len(Daysum_array)):
            Daysum_array[day,:,:][nanlocs] = np.nan
                
        Daysum_meanx = np.nanmean(Daysum_array,axis=1)
        Daysum_meany = np.nanmean(Daysum_meanx,axis=1)

        out_file_MB = 'DailyNARR_Prcp' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Daysum_meany) 
    print('daily prcp saved')
else:
    pass

if Calculate_Monthly_Distributed_Temp == True:
    print('Extracting Monthly distributed temp:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation

    allyears_jan_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_feb_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_mar_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_april_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_may_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_jun_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_july_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_aug_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_sept_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_oct_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_nov_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_dec_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    
    for year in years[1:]:
        #Monthly_container = np.empty(len(aice),12, len(T_array[1]),len(T_array[1,1,:])) #num of sims, num of months, x, y
        print(year)
        
        Temp_file_name = 'Temp' + str(domain) + str(year) + '.nc'
        Temp_file_in = os.path.join(NARR_PATH,Temp_file_name)
        
        MB_in = Dataset(Temp_file_in, "r")
        MB_var = 'Temperature'
        T_array = MB_in.variables[MB_var][:]
        
        #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        Day_mean = [np.nanmean(T_array[i:i+8], axis = 0) for i in range(0, len(T_array), 8)]
        Daysum_array = np.array(Day_mean) # convert daily net melt to an array
        
        for day in range(0,len(Daysum_array)):
            Daysum_array[day,:,:][nanlocs] = np.nan
        
        jan = Daysum_array[0:31,:,:]
        feb = Daysum_array[31:59,:,:]
        march = Daysum_array[59:90,:,:]
        april = Daysum_array[90:120,:,:]
        may = Daysum_array[120:151,:,:]
        june = Daysum_array[151:181,:,:]
        july = Daysum_array[181:212,:,:]
        aug = Daysum_array[212:243,:,:]
        sept = Daysum_array[243:273,:,:]
        octo = Daysum_array[273:304,:,:]
        nov = Daysum_array[304:334,:,:]
        dec = Daysum_array[334:365,:,:]   
        
        jan_xy = np.nanmean(jan,axis=0)
        feb_xy = np.nanmean(feb,axis=0)
        march_xy = np.nanmean(march,axis=0)
        april_xy = np.nanmean(april,axis=0)
        may_xy = np.nanmean(may,axis=0)
        june_xy = np.nanmean(june,axis=0)
        july_xy = np.nanmean(july,axis=0)
        aug_xy = np.nanmean(aug,axis=0)
        sept_xy = np.nanmean(sept,axis=0)
        oct_xy = np.nanmean(octo,axis=0)
        nov_xy = np.nanmean(nov,axis=0)
        dec_xy = np.nanmean(dec,axis=0)
        
        allyears_jan_container[year-2007,:,:] = jan_xy
        allyears_feb_container[year-2007,:,:] = feb_xy
        allyears_mar_container[year-2007,:,:] = march_xy
        allyears_april_container[year-2007,:,:] = april_xy
        allyears_may_container[year-2007,:,:] = may_xy
        allyears_jun_container[year-2007,:,:] = june_xy
        allyears_july_container[year-2007,:,:] = july_xy
        allyears_aug_container[year-2007,:,:] = aug_xy
        allyears_sept_container[year-2007,:,:] = sept_xy
        allyears_oct_container[year-2007,:,:] = oct_xy
        allyears_nov_container[year-2007,:,:] = nov_xy
        allyears_dec_container[year-2007,:,:] = dec_xy

    avg_jan = np.nanmean(allyears_jan_container,axis=0)
    avg_feb = np.nanmean(allyears_feb_container,axis=0)
    avg_mar = np.nanmean(allyears_mar_container,axis=0)
    avg_april = np.nanmean(allyears_april_container,axis=0)
    avg_may = np.nanmean(allyears_may_container,axis=0)
    avg_jun = np.nanmean(allyears_jun_container,axis=0)
    avg_july = np.nanmean(allyears_july_container,axis=0)
    avg_aug = np.nanmean(allyears_aug_container,axis=0)
    avg_sept = np.nanmean(allyears_sept_container,axis=0)
    avg_oct = np.nanmean(allyears_oct_container,axis=0)
    avg_nov = np.nanmean(allyears_nov_container,axis=0)
    avg_dec = np.nanmean(allyears_dec_container,axis=0)
    
    jan_outpath = os.path.join(OUTPUT_PATH,'jan_distributedTemp.npy')
    feb_outpath = os.path.join(OUTPUT_PATH,'feb_distributedTemp.npy')
    march_outpath = os.path.join(OUTPUT_PATH,'march_distributedTemp.npy')
    april_outpath = os.path.join(OUTPUT_PATH,'april_distributedTemp.npy')
    may_outpath = os.path.join(OUTPUT_PATH,'may_distributedTemp.npy')
    june_outpath = os.path.join(OUTPUT_PATH,'june_distributedTemp.npy')
    july_outpath = os.path.join(OUTPUT_PATH,'july_distributedTemp.npy')
    aug_outpath = os.path.join(OUTPUT_PATH,'aug_distributedTemp.npy')
    sept_outpath = os.path.join(OUTPUT_PATH,'sept_distributedTemp.npy')
    oct_outpath = os.path.join(OUTPUT_PATH,'oct_distributedTemp.npy')
    nov_outpath = os.path.join(OUTPUT_PATH,'nov_distributedTemp.npy')
    dec_outpath = os.path.join(OUTPUT_PATH,'dec_distributedTemp.npy')
    
    np.save(jan_outpath,avg_jan)
    np.save(feb_outpath,avg_feb)
    np.save(march_outpath,avg_mar)
    np.save(april_outpath,avg_april)
    np.save(may_outpath,avg_may)
    np.save(june_outpath,avg_jun)
    np.save(july_outpath,avg_july)
    np.save(aug_outpath,avg_aug)
    np.save(sept_outpath,avg_sept)
    np.save(oct_outpath,avg_oct)
    np.save(nov_outpath,avg_nov)
    np.save(dec_outpath,avg_dec)
    
    print('distributed temp calculated for all months')
else:
    print('Calculate Distributed Monthly Temp SKIPPED')
    
#-------------------------------------------------------------------------------
    
if Calculate_monthly_distributed_NARR_rain == True:
    print('Extracting Monthly distributed rain:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation

    allyears_jan_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_feb_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_mar_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_april_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_may_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_jun_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_july_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_aug_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_sept_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_oct_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_nov_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_dec_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    
    for year in years[1:]:
        #Monthly_container = np.empty(len(aice),12, len(T_array[1]),len(T_array[1,1,:])) #num of sims, num of months, x, y
        print(year)
        
        Prcp_file_name = 'netSnow' + str(domain) + str(year) + '.nc'
        Prcp_file_in = os.path.join(NARR_PATH,Prcp_file_name)
        
        MB_in = Dataset(Prcp_file_in, "r")
        #MB_var = 'Net snow'
        MB_var = 'Net snow'
        P_array = MB_in.variables[MB_var][:]
        
        Temp_file_name = 'Temp' + str(domain) + str(year) + '.nc'
        Temp_file_in = os.path.join(NARR_PATH,Temp_file_name)
        
        MB_in = Dataset(Temp_file_in, "r")
        MB_var = 'Temperature'
        T_array = MB_in.variables[MB_var][:]
        
        R2S = 1.0
        snowlocs = np.where(T_array <= R2S)
        P_array[snowlocs] = 0
        
        #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        Day_sum = [np.nansum(P_array[i:i+8], axis = 0) for i in range(0, len(P_array), 8)]
        Daysum_array = np.array(Day_sum) # convert daily net melt to an array
        
        for day in range(0,len(Daysum_array)):
            Daysum_array[day,:,:][nanlocs] = np.nan
        
        jan = Daysum_array[0:31,:,:]
        feb = Daysum_array[31:59,:,:]
        march = Daysum_array[59:90,:,:]
        april = Daysum_array[90:120,:,:]
        may = Daysum_array[120:151,:,:]
        june = Daysum_array[151:181,:,:]
        july = Daysum_array[181:212,:,:]
        aug = Daysum_array[212:243,:,:]
        sept = Daysum_array[243:273,:,:]
        octo = Daysum_array[273:304,:,:]
        nov = Daysum_array[304:334,:,:]
        dec = Daysum_array[334:365,:,:]   
        
        jan_xy = np.nanmean(jan,axis=0)
        feb_xy = np.nanmean(feb,axis=0)
        march_xy = np.nanmean(march,axis=0)
        april_xy = np.nanmean(april,axis=0)
        may_xy = np.nanmean(may,axis=0)
        june_xy = np.nanmean(june,axis=0)
        july_xy = np.nanmean(july,axis=0)
        aug_xy = np.nanmean(aug,axis=0)
        sept_xy = np.nanmean(sept,axis=0)
        oct_xy = np.nanmean(octo,axis=0)
        nov_xy = np.nanmean(nov,axis=0)
        dec_xy = np.nanmean(dec,axis=0)
        
        allyears_jan_container[year-2007,:,:] = jan_xy
        allyears_feb_container[year-2007,:,:] = feb_xy
        allyears_mar_container[year-2007,:,:] = march_xy
        allyears_april_container[year-2007,:,:] = april_xy
        allyears_may_container[year-2007,:,:] = may_xy
        allyears_jun_container[year-2007,:,:] = june_xy
        allyears_july_container[year-2007,:,:] = july_xy
        allyears_aug_container[year-2007,:,:] = aug_xy
        allyears_sept_container[year-2007,:,:] = sept_xy
        allyears_oct_container[year-2007,:,:] = oct_xy
        allyears_nov_container[year-2007,:,:] = nov_xy
        allyears_dec_container[year-2007,:,:] = dec_xy

    avg_jan = np.nanmean(allyears_jan_container,axis=0)
    avg_feb = np.nanmean(allyears_feb_container,axis=0)
    avg_mar = np.nanmean(allyears_mar_container,axis=0)
    avg_april = np.nanmean(allyears_april_container,axis=0)
    avg_may = np.nanmean(allyears_may_container,axis=0)
    avg_jun = np.nanmean(allyears_jun_container,axis=0)
    avg_july = np.nanmean(allyears_july_container,axis=0)
    avg_aug = np.nanmean(allyears_aug_container,axis=0)
    avg_sept = np.nanmean(allyears_sept_container,axis=0)
    avg_oct = np.nanmean(allyears_oct_container,axis=0)
    avg_nov = np.nanmean(allyears_nov_container,axis=0)
    avg_dec = np.nanmean(allyears_dec_container,axis=0)
    
    jan_outpath = os.path.join(OUTPUT_PATH,'jan_distributedRain.npy')
    feb_outpath = os.path.join(OUTPUT_PATH,'feb_distributedRain.npy')
    march_outpath = os.path.join(OUTPUT_PATH,'march_distributedRain.npy')
    april_outpath = os.path.join(OUTPUT_PATH,'april_distributedRain.npy')
    may_outpath = os.path.join(OUTPUT_PATH,'may_distributedRain.npy')
    june_outpath = os.path.join(OUTPUT_PATH,'june_distributedRain.npy')
    july_outpath = os.path.join(OUTPUT_PATH,'july_distributedRain.npy')
    aug_outpath = os.path.join(OUTPUT_PATH,'aug_distributedRain.npy')
    sept_outpath = os.path.join(OUTPUT_PATH,'sept_distributedRain.npy')
    oct_outpath = os.path.join(OUTPUT_PATH,'oct_distributedRain.npy')
    nov_outpath = os.path.join(OUTPUT_PATH,'nov_distributedRain.npy')
    dec_outpath = os.path.join(OUTPUT_PATH,'dec_distributedRain.npy')
    
    np.save(jan_outpath,avg_jan)
    np.save(feb_outpath,avg_feb)
    np.save(march_outpath,avg_mar)
    np.save(april_outpath,avg_april)
    np.save(may_outpath,avg_may)
    np.save(june_outpath,avg_jun)
    np.save(july_outpath,avg_july)
    np.save(aug_outpath,avg_aug)
    np.save(sept_outpath,avg_sept)
    np.save(oct_outpath,avg_oct)
    np.save(nov_outpath,avg_nov)
    np.save(dec_outpath,avg_dec)
    
    print('distributed NARR rain calculated for all months')
else:
    print('Calculate Distributed Monthly NARR rain SKIPPED')
    
#-------------------------------------------------------------------------------
    
if Calculate_monthly_distributed_NARR_precip == True:
    print('Extracting Monthly distributed precip:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation

    allyears_jan_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_feb_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_mar_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_april_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_may_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_jun_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_july_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_aug_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_sept_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_oct_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_nov_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_dec_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    
    for year in years[1:]:
        #Monthly_container = np.empty(len(aice),12, len(T_array[1]),len(T_array[1,1,:])) #num of sims, num of months, x, y
        print(year)
        
        Prcp_file_name = 'netSnow' + str(domain) + str(year) + '.nc'
        Prcp_file_in = os.path.join(NARR_PATH,Prcp_file_name)
        
        MB_in = Dataset(Prcp_file_in, "r")
        MB_var = 'Net snow'
        P_array = MB_in.variables[MB_var][:]
        
        #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        Day_sum = [np.nansum(P_array[i:i+8], axis = 0) for i in range(0, len(P_array), 8)]
        Daysum_array = np.array(Day_sum) # convert daily net melt to an array
        
        for day in range(0,len(Daysum_array)):
            Daysum_array[day,:,:][nanlocs] = np.nan
        
        jan = Daysum_array[0:31,:,:]
        feb = Daysum_array[31:59,:,:]
        march = Daysum_array[59:90,:,:]
        april = Daysum_array[90:120,:,:]
        may = Daysum_array[120:151,:,:]
        june = Daysum_array[151:181,:,:]
        july = Daysum_array[181:212,:,:]
        aug = Daysum_array[212:243,:,:]
        sept = Daysum_array[243:273,:,:]
        octo = Daysum_array[273:304,:,:]
        nov = Daysum_array[304:334,:,:]
        dec = Daysum_array[334:365,:,:]   
        
        jan_xy = np.nanmean(jan,axis=0)
        feb_xy = np.nanmean(feb,axis=0)
        march_xy = np.nanmean(march,axis=0)
        april_xy = np.nanmean(april,axis=0)
        may_xy = np.nanmean(may,axis=0)
        june_xy = np.nanmean(june,axis=0)
        july_xy = np.nanmean(july,axis=0)
        aug_xy = np.nanmean(aug,axis=0)
        sept_xy = np.nanmean(sept,axis=0)
        oct_xy = np.nanmean(octo,axis=0)
        nov_xy = np.nanmean(nov,axis=0)
        dec_xy = np.nanmean(dec,axis=0)
        
        allyears_jan_container[year-2007,:,:] = jan_xy
        allyears_feb_container[year-2007,:,:] = feb_xy
        allyears_mar_container[year-2007,:,:] = march_xy
        allyears_april_container[year-2007,:,:] = april_xy
        allyears_may_container[year-2007,:,:] = may_xy
        allyears_jun_container[year-2007,:,:] = june_xy
        allyears_july_container[year-2007,:,:] = july_xy
        allyears_aug_container[year-2007,:,:] = aug_xy
        allyears_sept_container[year-2007,:,:] = sept_xy
        allyears_oct_container[year-2007,:,:] = oct_xy
        allyears_nov_container[year-2007,:,:] = nov_xy
        allyears_dec_container[year-2007,:,:] = dec_xy

    avg_jan = np.nanmean(allyears_jan_container,axis=0)
    avg_feb = np.nanmean(allyears_feb_container,axis=0)
    avg_mar = np.nanmean(allyears_mar_container,axis=0)
    avg_april = np.nanmean(allyears_april_container,axis=0)
    avg_may = np.nanmean(allyears_may_container,axis=0)
    avg_jun = np.nanmean(allyears_jun_container,axis=0)
    avg_july = np.nanmean(allyears_july_container,axis=0)
    avg_aug = np.nanmean(allyears_aug_container,axis=0)
    avg_sept = np.nanmean(allyears_sept_container,axis=0)
    avg_oct = np.nanmean(allyears_oct_container,axis=0)
    avg_nov = np.nanmean(allyears_nov_container,axis=0)
    avg_dec = np.nanmean(allyears_dec_container,axis=0)
    
    jan_outpath = os.path.join(OUTPUT_PATH,'jan_distributedPrecip.npy')
    feb_outpath = os.path.join(OUTPUT_PATH,'feb_distributedPrecip.npy')
    march_outpath = os.path.join(OUTPUT_PATH,'march_distributedPrecip.npy')
    april_outpath = os.path.join(OUTPUT_PATH,'april_distributedPrecip.npy')
    may_outpath = os.path.join(OUTPUT_PATH,'may_distributedPrecip.npy')
    june_outpath = os.path.join(OUTPUT_PATH,'june_distributedPrecip.npy')
    july_outpath = os.path.join(OUTPUT_PATH,'july_distributedPrecip.npy')
    aug_outpath = os.path.join(OUTPUT_PATH,'aug_distributedPrecip.npy')
    sept_outpath = os.path.join(OUTPUT_PATH,'sept_distributedPrecip.npy')
    oct_outpath = os.path.join(OUTPUT_PATH,'oct_distributedPrecip.npy')
    nov_outpath = os.path.join(OUTPUT_PATH,'nov_distributedPrecip.npy')
    dec_outpath = os.path.join(OUTPUT_PATH,'dec_distributedPrecip.npy')
    
    np.save(jan_outpath,avg_jan)
    np.save(feb_outpath,avg_feb)
    np.save(march_outpath,avg_mar)
    np.save(april_outpath,avg_april)
    np.save(may_outpath,avg_may)
    np.save(june_outpath,avg_jun)
    np.save(july_outpath,avg_july)
    np.save(aug_outpath,avg_aug)
    np.save(sept_outpath,avg_sept)
    np.save(oct_outpath,avg_oct)
    np.save(nov_outpath,avg_nov)
    np.save(dec_outpath,avg_dec)
    
    print('distributed NARR precip calculated for all months')
else:
    print('Calculate Distributed Monthly NARR precip SKIPPED')
    
#---------------------------------------------------------------------------------
if Calculate_Monthly_Distributed_Snow == True:
    print('Extracting Monthly distributed snow:')
    #allruns_MB = np.empty((len(aice),len(Zgrid),len(Zgrid[0]))) #an empty array with a 'slot'for each simulation

    allyears_jan_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_feb_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_mar_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_april_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_may_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_jun_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_july_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_aug_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_sept_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_oct_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_nov_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    allyears_dec_container = np.empty((len(years[1:]),len(Zgrid),len(Zgrid[0])))
    for year in years[1:]:
        #Monthly_container = np.empty(len(aice),12, len(T_array[1]),len(T_array[1,1,:])) #num of sims, num of months, x, y
        print(year)
        
        sims_jan_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_feb_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_mar_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_april_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_may_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_june_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_july_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_aug_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_sept_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_oct_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_nov_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        sims_dec_container = np.empty((len(aice),len(Zgrid),len(Zgrid[0])))
        for sim in range(0, len(aice)):
            print('starting sim ' + str(sim))
        
            #MB_file_name = 'NetMelt' + str(year) + str(sim) + '.nc'
            MB_file_name = 'Accumulation' + str(year) + str(sim) + '.nc'
            MB_File_in = os.path.join(FILEPATH,MB_file_name)
    
            MB_in = Dataset(MB_File_in, "r")
            #print('file in')
            MB_var = 'Accumulation'
            MB_array = MB_in.variables[MB_var][:]
            #print('variable in')
            
            #N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
            Day_Sum = [np.nansum(MB_array[i:i+8], axis = 0) for i in range(0, len(MB_array), 8)]
            Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
            for day in range(0,len(Daysum_array)):
                Daysum_array[day,:,:][nanlocs] = np.nan
                
            
            jan = Daysum_array[0:31,:,:]
            feb = Daysum_array[31:59,:,:]
            march = Daysum_array[59:90,:,:]
            april = Daysum_array[90:120,:,:]
            may = Daysum_array[120:151,:,:]
            june = Daysum_array[151:181,:,:]
            july = Daysum_array[181:212,:,:]
            aug = Daysum_array[212:243,:,:]
            sept = Daysum_array[243:273,:,:]
            octo = Daysum_array[273:304,:,:]
            nov = Daysum_array[304:334,:,:]
            dec = Daysum_array[334:365,:,:]   
            
            jan_xy = np.nanmean(jan,axis=0)
            feb_xy = np.nanmean(feb,axis=0)
            march_xy = np.nanmean(march,axis=0)
            april_xy = np.nanmean(april,axis=0)
            may_xy = np.nanmean(may,axis=0)
            june_xy = np.nanmean(june,axis=0)
            july_xy = np.nanmean(july,axis=0)
            aug_xy = np.nanmean(aug,axis=0)
            sept_xy = np.nanmean(sept,axis=0)
            oct_xy = np.nanmean(octo,axis=0)
            nov_xy = np.nanmean(nov,axis=0)
            dec_xy = np.nanmean(dec,axis=0)
            
            sims_jan_container[sim,:,:] = jan_xy
            sims_feb_container[sim,:,:] = feb_xy
            sims_mar_container[sim,:,:] = march_xy
            sims_april_container[sim,:,:] = april_xy
            sims_may_container[sim,:,:] = may_xy
            sims_june_container[sim,:,:] = june_xy
            sims_july_container[sim,:,:] = july_xy
            sims_aug_container[sim,:,:] = aug_xy
            sims_sept_container[sim,:,:] = sept_xy
            sims_oct_container[sim,:,:] = oct_xy
            sims_nov_container[sim,:,:] = nov_xy
            sims_dec_container[sim,:,:] = dec_xy
            
        allsims_avg_jan = np.nanmean(sims_jan_container,axis=0)
        allsims_avg_feb = np.nanmean(sims_feb_container,axis=0)
        allsims_avg_mar = np.nanmean(sims_mar_container,axis=0)
        allsims_avg_april = np.nanmean(sims_april_container,axis=0)
        allsims_avg_may = np.nanmean(sims_may_container,axis=0)
        allsims_avg_jun = np.nanmean(sims_june_container,axis=0)
        allsims_avg_july = np.nanmean(sims_july_container,axis=0)
        allsims_avg_aug = np.nanmean(sims_aug_container,axis=0)
        allsims_avg_sept = np.nanmean(sims_sept_container,axis=0)
        allsims_avg_oct = np.nanmean(sims_oct_container,axis=0)
        allsims_avg_nov = np.nanmean(sims_nov_container,axis=0)
        allsims_avg_dec = np.nanmean(sims_dec_container,axis=0)
        
        allyears_jan_container[year-2007,:,:] = allsims_avg_jan
        allyears_feb_container[year-2007,:,:] = allsims_avg_feb
        allyears_mar_container[year-2007,:,:] = allsims_avg_mar
        allyears_april_container[year-2007,:,:] = allsims_avg_april
        allyears_may_container[year-2007,:,:] = allsims_avg_may
        allyears_jun_container[year-2007,:,:] = allsims_avg_jun
        allyears_july_container[year-2007,:,:] = allsims_avg_july
        allyears_aug_container[year-2007,:,:] = allsims_avg_aug
        allyears_sept_container[year-2007,:,:] = allsims_avg_sept
        allyears_oct_container[year-2007,:,:] = allsims_avg_oct
        allyears_nov_container[year-2007,:,:] = allsims_avg_nov
        allyears_dec_container[year-2007,:,:] = allsims_avg_dec
    
    avg_jan = np.nanmean(allyears_jan_container,axis=0)
    avg_feb = np.nanmean(allyears_feb_container,axis=0)
    avg_mar = np.nanmean(allyears_mar_container,axis=0)
    avg_april = np.nanmean(allyears_april_container,axis=0)
    avg_may = np.nanmean(allyears_may_container,axis=0)
    avg_jun = np.nanmean(allyears_jun_container,axis=0)
    avg_july = np.nanmean(allyears_july_container,axis=0)
    avg_aug = np.nanmean(allyears_aug_container,axis=0)
    avg_sept = np.nanmean(allyears_sept_container,axis=0)
    avg_oct = np.nanmean(allyears_oct_container,axis=0)
    avg_nov = np.nanmean(allyears_nov_container,axis=0)
    avg_dec = np.nanmean(allyears_dec_container,axis=0)
    
    jan_outpath = os.path.join(OUTPUT_PATH,'jan_distributedSnow.npy')
    feb_outpath = os.path.join(OUTPUT_PATH,'feb_distributedSnow.npy')
    march_outpath = os.path.join(OUTPUT_PATH,'march_distributedSnow.npy')
    april_outpath = os.path.join(OUTPUT_PATH,'april_distributedSnow.npy')
    may_outpath = os.path.join(OUTPUT_PATH,'may_distributedSnow.npy')
    june_outpath = os.path.join(OUTPUT_PATH,'june_distributedSnow.npy')
    july_outpath = os.path.join(OUTPUT_PATH,'july_distributedSnow.npy')
    aug_outpath = os.path.join(OUTPUT_PATH,'aug_distributedSnow.npy')
    sept_outpath = os.path.join(OUTPUT_PATH,'sept_distributedSnow.npy')
    oct_outpath = os.path.join(OUTPUT_PATH,'oct_distributedSnow.npy')
    nov_outpath = os.path.join(OUTPUT_PATH,'nov_distributedSnow.npy')
    dec_outpath = os.path.join(OUTPUT_PATH,'dec_distributedSnow.npy')
    
    np.save(jan_outpath,avg_jan)
    np.save(feb_outpath,avg_feb)
    np.save(march_outpath,avg_mar)
    np.save(april_outpath,avg_april)
    np.save(may_outpath,avg_may)
    np.save(june_outpath,avg_jun)
    np.save(july_outpath,avg_july)
    np.save(aug_outpath,avg_aug)
    np.save(sept_outpath,avg_sept)
    np.save(oct_outpath,avg_oct)
    np.save(nov_outpath,avg_nov)
    np.save(dec_outpath,avg_dec)
    
    print('distributed snow calculated for all months')
else:
    print('Calculate Distributed Monthly Snow SKIPPED')
    
#------------------------------------------------------------------------------
if Calculate_daily_BC_Rain == True:
    print('calculating daily bias corrected RAIN')
    for year in years[1:]:
        if year == 2008 or year == 2012 or year == 2016:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)+1))
        else:
            Yearly_NetMelt_container = np.empty((len(aice),int(len(T_array)/8)))
        print(year)
        
        Temp_file_name = 'Rain' + str(domain) + '_BC_' + str(year) +'.nc'
        Temp_file_in = os.path.join(NARR_PATH,Temp_file_name)
        
        MB_in = Dataset(Temp_file_in, "r")
        #print('file in')
        MB_var = 'Rain'
        R_array = MB_in.variables[MB_var][:]
        
        
        Day_Sum = [np.nansum(R_array[i:i+8], axis = 0) for i in range(0, len(R_array), 8)]
        Daysum_array = np.array(Day_Sum) # convert daily net melt to an array
            
        for day in range(0,len(Daysum_array)):
            Daysum_array[day,:,:][nanlocs] = np.nan
                
        Daysum_meanx = np.nanmean(Daysum_array,axis=1)
        Daysum_meany = np.nanmean(Daysum_meanx,axis=1)

        out_file_MB = 'Daily_Rain' + str(year) + '.npy'
        MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
        np.save(MB_out_path, Daysum_meany) 
    print('daily prcp saved')
else:
    pass

if Calculate_distributed_NARR_NoBCRain == True:
    print('extracting distributed precip')

    MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
    for year in years[1:]: # start at index 1 (first year skipped for spin-up)
        print(year)
        
        MB_file_name = 'Rain' + str(domain) + str(year) + '.nc'
        MB_File_in = os.path.join(NARR_PATH,MB_file_name)

        MB_in = Dataset(MB_File_in, "r")
        #print('file in')
        MB_var = 'Rain'
        MB_array = MB_in.variables[MB_var][:]
        
        
        N_MB = np.sum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        
        # insert the net mass balance for the year into the empty MB container for that year
        MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.

    Temp_allyears_mean = np.nanmean(MB_empty_container, axis =0)
    out_file_MB = 'allrunRain_' + 'DownscaledNoBC' + '.npy'
    MB_out_path = os.path.join(NARR_PATH,out_file_MB)  
    np.save(MB_out_path, Temp_allyears_mean) #results in a file caled Net_MB1.npy for each simulation

else:
    pass
    
#------------------------------------------------------------------------------
if Calculate_distributed_NARR_precip_withBC == True:
    print('extracting distributed precip')

    MB_empty_container = np.empty(((len(years)-1),len(Zgrid),len(Zgrid[0])))
    for year in years[1:]: # start at index 1 (first year skipped for spin-up)
        print(year)
        
        MB_file_name = 'Prcp' + str(domain) + '_BC_' + str(year) + '.nc'
        MB_File_in = os.path.join(NARR_PATH,MB_file_name)

        MB_in = Dataset(MB_File_in, "r")
        #print('file in')
        MB_var = 'Net snow'
        MB_array = MB_in.variables[MB_var][:]
        
        
        N_MB = np.sum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        
        # insert the net mass balance for the year into the empty MB container for that year
        MB_empty_container[np.where(np.asarray(years) == year-1),:,:] = N_MB # np.asarray Converts the input to an array.

    Temp_allyears_mean = np.nanmean(MB_empty_container, axis =0)
    out_file_MB = 'allrunPrecipitation_' + 'DownscaledBC' + '.npy'
    MB_out_path = os.path.join(OUTPUT_PATH,out_file_MB)  
    np.save(MB_out_path, Temp_allyears_mean) #results in a file caled Net_MB1.npy for each simulation

else:
    pass
#------------------------------------------------------------------------------