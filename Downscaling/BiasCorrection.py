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
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import netcdf_container_gen
from Model_functions_ver4 import model_domain


#initialize model --------------------------------------------------------
BiasCorrect_Temp = True 
# NOTE: TEMP BIAS CORRECTION MUST HAPPEN FIRST : since R2S threshold is based on bias corrected T_array
BiasCorrect_Prcp = True
BiasCorrect_Rain = False #take daily rain arrays and bias correct them w\ the accumulation elevation dependent BC

R2S = 1.0
start_year = 2007
end_year = 2018

Glacier_name = 'kaskawulsh'
Glacier_id = 'kaskonly'
catchment = True
INPUT_PATH = 'D:/Downscaled_files/Kaskonly_staticZ/missing_trib/DownscaledNARR_2006-2018' #where are the downscaled NON-BIAS CORRECTED files located
OUTPUT_PATH = 'D:/BiasCorrected_files/Kaskonly/MissingCAtrib/NewBiasCorrection_2007-2020'

File_suffix = '.nc'
reffiles = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Ref_files'

years = []
year = start_year
while year <= end_year:
  years.append(year)
  year = year + 1

##-------Turn glacier grid vectors into 3D grids--------------------##
Zgrid, Xgrid, Ygrid, xbounds, ybounds, Sfc = model_domain(catchment)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#Set up values for bias correction; given array are DT values developped indepedently of
#model for time period, interpolated linearly between months at a daily 
#resolution 
#DP values are elevation dependent bias correction factors for accumulation 

DTval = interp1d(np.asarray([1.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,366.]), \
                 np.asarray([-4.86,-3.98,-2.1,0.03,0.89,1.42,0.94,0.32,-0.32,-1.72,-4.37,-5.22, -4.86]), \
                 kind = 'linear')   
    
#DPval = interp1d(np.asarray([500, 900, 1100, 1300, 1500, 1700, 2100, 2300, 2500, \
#        2700, 3100, 3700, 5000]), np.asarray([ 1.28530635,  1.26180458,  1.23360247,  1.25240269,  1.28293036,
#        1.26320929,  1.24449735,  1.30082501,  2.46677552,  4.24804321,
#        5.45333788,  5.45333788, 5.45333788]), kind = 'linear')
  
# This is the new bias correction (DP_bin450 from Recreate_BiasCorrection.py)  
DPval = interp1d(np.asarray([500, 1225.0, 1675.0, 2125.0, 2575.0, 5000]), np.asarray([1.50821731,1.50821731,1.27023925,1.29709692,2.42742537,2.42742537]), kind = 'linear')    
    
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
        #T_var = 'Precipitation'
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:] #T_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        # prepare the output file (to contain downscaled and Bias Corrected temp in netcdf form)
        Temp_output_name = 'Temp_' + Glacier_name + '_BC_' + str(year)
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
        netcdf_container_gen(T_array, 'Temperature', Temp_output_path, File_suffix, ybounds, xbounds, year, reffiles)
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
        Just_Downscaled_Temp = 'Temp_' + Glacier_name + '_BC_' + str(year) + '.nc'
        Downscaled_Temp = os.path.join(OUTPUT_PATH,Just_Downscaled_Temp)
        inT = Dataset(Downscaled_Temp,'r')
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:] #T_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        rainlocs = np.where(T_array > R2S)
        
        # call in the downscaled accumulation file (must be already converted to netcdf)
        Prcp_downscaled_filename = 'netSnow' + 'kaskonly' + str(year) + '.nc' #can't use this bc it has already been bias corrected
        Prcp_input = os.path.join(INPUT_PATH,Prcp_downscaled_filename)
        inP = Dataset(Prcp_input,'r')
        P_var = 'Precipitation'
        P_array = inP.variables[P_var][:] #P_array has shape (2920,218,328) = (t,x,y)
        sys.stdout.flush()
        
        # fix one anomalous timestep in 2014
        if year == 2014:
            P_array[1744][np.where(P_array[1744]>0)] = 0
        
        # get rain locations based on T_array
        
        Precip_uncorrected_array = np.zeros(P_array.shape)
        Precip_uncorrected_array[:,:,:] = P_array[:,:,:]
        Rain_array = np.zeros(P_array.shape)
        
        #print(T_array[1440,100,100])
        #print(P_array[1440,100,100]) 
        #print(Precip_uncorrected_array[1440,100,100])
        
        # prepare the output file (to contain downscaled and Bias Corrected temp in netcdf form)
        Prcp_output_name = 'Snow_' + Glacier_name + '_BC_' + str(year) 
        Prcp_output_path = os.path.join(OUTPUT_PATH,Prcp_output_name)
        
        Rain_output_name = 'Rain_' + Glacier_name + '_' + str(year)
        Rain_output_path = os.path.join(OUTPUT_PATH,Rain_output_name)
        
        # Calculate bias correction factor for each elevation and apply it to the P_array
        for x in range(len(Zgrid)):
            for y in range(len(Zgrid[0])):
                z = Zgrid[x][y] #find out the elevation at each x,y coordinate
                DeltaC = DPval(z) # find the Bias Correction factor (DeltaC)
                P_array[:,x,y] = P_array[:,x,y]*(DeltaC)
        
        #print(P_array[1440,100,100])
        # replace P_array with un-BC'd rainlocs
        P_array[rainlocs] = Precip_uncorrected_array[rainlocs]
        Rain_array[rainlocs] = Precip_uncorrected_array[rainlocs] 
        #print(P_array[1440,100,100]) 

                
    #save the new P_array to a netcdf file#------------------------------------------
        netcdf_container_gen(P_array, 'Snow', Prcp_output_path, File_suffix, ybounds, xbounds, year, reffiles)
        #netcdf_container_gen(Rain_array, 'Rain', Rain_output_path, File_suffix, ybounds, xbounds, year, reffiles)
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
        Just_Downscaled_Temp = 'Temp_' + Glacier_name + '_BC_' + str(year) + '.nc'
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
        netcdf_container_gen(Rain_array, 'Rain', Prcp_output_path, File_suffix, ybounds, xbounds, year, reffiles)
        print('Rain bias corrected for ' + str(year))
        
    print('ALL YEARS BIAS CORRECTED FOR Rain')
else:
    print('Rain NOT bias corrected')