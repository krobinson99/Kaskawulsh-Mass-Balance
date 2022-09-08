# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 09:22:23 2021

@author: katierobinson
"""
#NOTE: input files (downscaled temp and precip) should already be in NETCDF
# form before the bias correction can be applied

import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import sys
import os
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
import Model_functions_ver4
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import netcdf_container_gen


#initialize model --------------------------------------------------------
BiasCorrect_Temp = True 
# NOTE: TEMP BIAS CORRECTION MUST HAPPEN FIRST : since R2S threshold is based on bias corrected T_array
BiasCorrect_Prcp = True
BiasCorrect_Rain = False #take daily rain arrays and bias correct them w\ the accumulation elevation dependent BC

R2S = 1.0
start_year = 2021
end_year = 2022

Glacier_id = 'kaskonly'
File_glacier_in = 'kaskonly.txt'
considering_catchment = False
#INPUT_PATH = 'F:\\Mass Balance Model\\DownscaledNARR' #where are the downscaled NON-BIAS CORRECTED files located
INPUT_PATH = 'F:\\Mass Balance Model\\DownscaledNARR_1979-2022\\Constant_Z'
OUTPUT_PATH = 'F:\\Mass Balance Model\\BiasCorrectedNARR_1979-2022\\Constant_Z'
#OUTPUT_PATH = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'
#RAIN_PATH = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
File_suffix = '.nc'

years = []
year = start_year
while year <= end_year:
  years.append(year)
  year = year + 1

##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')

if considering_catchment == True:
    Ix = glacier[:,4]
    Iy = glacier[:,5] 
    Ih = glacier[:,6]    
else:
    Ix = glacier[:,3]
    Iy = glacier[:,4] 
    Ih = glacier[:,5]         
        
Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#Set up values for bias correction; given array are DT values developped indepedently of
#model for time period, interpolated linearly between months at a daily 
#resolution 
#DP values are elevation dependent bias correction factors for accumulation 

DTval = interp1d(np.asarray([1.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,366.]), \
                 np.asarray([-4.86,-3.98,-2.1,0.03,0.89,1.42,0.94,0.32,-0.32,-1.72,-4.37,-5.22, -4.86]), \
                 kind = 'linear')   
    
DPval = interp1d(np.asarray([500, 900, 1100, 1300, 1500, 1700, 2100, 2300, 2500, \
        2700, 3100, 3700, 5000]), np.asarray([ 1.28530635,  1.26180458,  1.23360247,  1.25240269,  1.28293036,
        1.26320929,  1.24449735,  1.30082501,  2.46677552,  4.24804321,
        5.45333788,  5.45333788, 5.45333788]), kind = 'linear')

#--------------------------------------------------------------------------
#####Apply TEMPERATURE bias correction first:########
if BiasCorrect_Temp == True:
    print('Beginning Temperature Bias Correction')
    for year in years:
        print('Starting year: ' + str(year))
        
        # call in the downscaled temp file (must be already converted to netcdf)
        Temp_downscaled_filename = 'Temp' + Glacier_id + str(year) + '.nc' #can't use this bc it has already been bias corrected
        Temp_input = os.path.join(INPUT_PATH,Temp_downscaled_filename)
        inT = Dataset(Temp_input,'r')
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:] #T_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        # prepare the output file (to contain downscaled and Bias Corrected temp in netcdf form)
        Temp_output_name = 'Temp' + Glacier_id + '_BC_' + str(year)
        Temp_output_path = os.path.join(OUTPUT_PATH,Temp_output_name)
        
        # calculate bias correction factor for each timestep(every 3 hours) and apply it to the T_array
        index = 0
        timesteps = []
        while index < (len(T_array)):
            timesteps = np.linspace(1,(len(T_array)/8),len(T_array)) #applies bias correction at each 3hourly step
            BC_factor = DTval(timesteps[index])
            T_array[index][:][:] = T_array[index][:][:] + BC_factor #apply the bias correction here
            index+=1

    #save the new T_array to a netcdf file#------------------------------------------
        netcdf_container_gen(T_array, 'Temperature', Temp_output_path, File_suffix, ybounds, xbounds, year)
        print('Temperature bias corrected for ' + str(year))
        
    print('ALL YEARS BIAS CORRECTED FOR TEMPERATURE')
else:
    print('Temperature NOT bias corrected')
#------------------------------------------------------------------------------
###Begin Accumulation bias correction
#--------------------------------------------------------------------------
if BiasCorrect_Prcp == True:
    print('Beginning Accumulation Bias Correction')
    for year in years:
        print('Starting year: ' + str(year))
        
        # call in the downscaled temp file (must be already converted to netcdf)
        Just_Downscaled_Temp = 'Temp' + Glacier_id + '_BC_' + str(year) + '.nc'
        Downscaled_Temp = os.path.join(OUTPUT_PATH,Just_Downscaled_Temp)
        inT = Dataset(Downscaled_Temp,'r')
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:] #T_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        rainlocs = np.where(T_array > R2S)
        
        # call in the downscaled accumulation file (must be already converted to netcdf)
        Prcp_downscaled_filename = 'netSnow' + Glacier_id + str(year) + '.nc' #can't use this bc it has already been bias corrected
        Prcp_input = os.path.join(INPUT_PATH,Prcp_downscaled_filename)
        inP = Dataset(Prcp_input,'r')
        P_var = 'Temperature'
        P_array = inP.variables[P_var][:] #P_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        # get rain locations based on T_array
        
        Precip_uncorrected_array = np.zeros(P_array.shape)
        Precip_uncorrected_array[:,:,:] = P_array[:,:,:]
        Rain_array = np.zeros(P_array.shape)
        
        #print(T_array[1440,100,100])
        #print(P_array[1440,100,100]) 
        #print(Precip_uncorrected_array[1440,100,100])
        
        # prepare the output file (to contain downscaled and Bias Corrected temp in netcdf form)
        Prcp_output_name = 'Prcp' + Glacier_id + '_BC_' + str(year) 
        Prcp_output_path = os.path.join(OUTPUT_PATH,Prcp_output_name)
        
        Rain_output_name = 'Rain' + Glacier_id + str(year)
        Rain_output_path = os.path.join(OUTPUT_PATH,Rain_output_name)
        
        # Calculate bias correction factor for each elevation and apply it to the P_array
        for x in range(len(Zgrid)):
            for y in range(len(Zgrid[0])):
                z = Zgrid[x][y] #find out the elevation at each x,y coordinate
                DeltaC = DPval(z) # find the Bias Correction factor (DeltaC)
                P_array[:,x,y] = P_array[:,x,y]*(DeltaC)
        
        #print(P_array[1440,100,100])
        # replace P_array with un BC'd rainlocs
        P_array[rainlocs] = Precip_uncorrected_array[rainlocs]
        Rain_array[rainlocs] = Precip_uncorrected_array[rainlocs] 
        #print(P_array[1440,100,100]) 

                
    #save the new P_array to a netcdf file#------------------------------------------
        netcdf_container_gen(P_array, 'Net snow', Prcp_output_path, File_suffix, ybounds, xbounds, year)
        netcdf_container_gen(Rain_array, 'Rain', Rain_output_path, File_suffix, ybounds, xbounds, year)
        print('Accumulation bias corrected for ' + str(year))
        
    print('ALL YEARS BIAS CORRECTED FOR Precipitation')
else:
    print('Precipitation NOT bias corrected')

#------------------------------------------------------------------------------
###Begin RAIN bias correction
#--------------------------------------------------------------------------
if BiasCorrect_Rain == True:
    print('Beginning Rain Bias Correction')
    for year in years:
        print('Starting year: ' + str(year))
        
        # call in the downscaled temp file (must be already converted to netcdf)
        Just_Downscaled_Temp = 'Temp' + Glacier_id + '_BC_' + str(year) + '.nc'
        Downscaled_Temp = os.path.join(OUTPUT_PATH,Just_Downscaled_Temp)
        inT = Dataset(Downscaled_Temp,'r')
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:] #T_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        rainlocs = np.where(T_array > R2S)
        
        # call in the downscaled accumulation file (must be already converted to netcdf)
        Prcp_downscaled_filename = 'netSnow' + Glacier_id + str(year) + '.nc' #can't use this bc it has already been bias corrected
        Prcp_input = os.path.join(INPUT_PATH,Prcp_downscaled_filename)
        inP = Dataset(Prcp_input,'r')
        P_var = 'Temperature'
        P_array = inP.variables[P_var][:] #P_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        # get rain locations based on T_array
        
        Rain_array = np.zeros(P_array.shape)
        
        # prepare the output file (to contain downscaled and Bias Corrected temp in netcdf form)
        Prcp_output_name = 'Rain' + Glacier_id + '_BC_' + str(year) 
        Prcp_output_path = os.path.join(OUTPUT_PATH,Prcp_output_name)
        
        # Calculate bias correction factor for each elevation and apply it to the P_array
        for x in range(len(Zgrid)):
            for y in range(len(Zgrid[0])):
                z = Zgrid[x][y] #find out the elevation at each x,y coordinate
                DeltaC = DPval(z) # find the Bias Correction factor (DeltaC)
                P_array[:,x,y] = P_array[:,x,y]*(DeltaC)
        
        #print(P_array[1440,100,100])
        # replace P_array with un BC'd rainlocs
        Rain_array[rainlocs] = P_array[rainlocs]
        #print(P_array[1440,100,100]) 

                
    #save the new P_array to a netcdf file#------------------------------------------
        netcdf_container_gen(Rain_array, 'Rain', Prcp_output_path, File_suffix, ybounds, xbounds, year)
        print('Rain bias corrected for ' + str(year))
        
    print('ALL YEARS BIAS CORRECTED FOR Rain')
else:
    print('Rain NOT bias corrected')