# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:02:43 2021

@author: agribbon
"""
######IMPORTANT########
#                     #
# RUN MBMnamelist.py  # 
# BEFORE running this #
#      script!        #
#                     #
#######################

#This is a copy of the script MB_MODEL_FINALRUNS_balancefluxes.py (originally by Erik Young) with revisions by Katie Robinson 

#import libraries
import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime as dt
from datetime import datetime
import sys
import os
# import model functions
from Model_functions_ver4 import MB_vectorized_discreteSnI
from Model_functions_ver4 import cold_content
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import get_meanSP
from Model_functions_ver4 import netcdf_container_gen
# import debris functions
sys.path.insert(1,'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\DebrisThickness')
from DebrisFunctions import MeltFactors
# import from the config file
from MBMnamelist import glacier_id
from MBMnamelist import params_filename
from MBMnamelist import start_year
from MBMnamelist import start_day
from MBMnamelist import end_year
from MBMnamelist import debris
from MBMnamelist import time_step
from MBMnamelist import Output_path
from MBMnamelist import Rain_to_snow as R2S
from MBMnamelist import Refreezing
from MBMnamelist import T_inputs
from MBMnamelist import P_inputs
from MBMnamelist import SR_inputs
from MBMnamelist import Temp_shift
from MBMnamelist import temp_shift_factor
from MBMnamelist import Bias_CorrectionT as BC_T
from MBMnamelist import Bias_CorrectionP as BC_P
from MBMnamelist import Considering_Catchment
from MBMnamelist import transition_thickness as tt
from MBMnamelist import debris_treatment
from MBMnamelist import debris_thickness_map

#from MBMnamelist import *

#Initialize model (get parameterizations from the namelist)
sim = -1
glacier_ID = glacier_id #namelist!!
OUTPUT_PATH = Output_path #namelist
#R2S = Rain_to_snow #rain to snow melt threshold
REFREEZING = Refreezing
#BC_T = Bias_CorrectionT
#BC_P = Bias_CorrectionP
considering_catchment = Considering_Catchment
#where are the downscaled/bias corrected inputs stored
Temp_input_path = T_inputs
Precip_input_path = P_inputs
Srad_input_path = SR_inputs

#Load parameters 
params = np.loadtxt(params_filename) #namelist!!
aice_p = params[0,:]
asnow_p = params[1,:]
MF_p = params[2,:]

## set up time range ###
years = []
year = start_year
while year <= end_year:
  years.append(year)
  year = year + 1
  
#-------------------------#
#calculate mean snowpack for cold content before looping
DHval = get_meanSP(years,glacier_ID,R2S,BC_T,BC_P,Temp_input_path,Precip_input_path) #change this func to match shape of domain
#-------------------------#

#-------Set up glacier vectors-----------------------

if considering_catchment == True:
    File_glacier_in = 'kask_catchment.txt'
else:
    File_glacier_in = glacier_ID + '_deb.txt'

glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')

if considering_catchment == True:
    if debris == True:
        Ix = glacier[:,4] 
        Iy = glacier[:,5] 
        Ih = glacier[:,6]
        sfc_type = glacier[:,8]
        debris_array = glacier[:,9]
    else:      
        Ix = glacier[:,4] 
        Iy = glacier[:,5] 
        Ih = glacier[:,6]
        sfc_type = glacier[:,8]
else:
    if debris == True:
        Ix = glacier[:,3] 
        Iy = glacier[:,4] 
        Ih = glacier[:,2] 
        debris_array = glacier[:,6]
    else:  
        Ix = glacier[:,3] 
        Iy = glacier[:,4] 
        Ih = glacier[:,2]

#-------Turn vectors into 3D gridded inputs--------------------

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs = np.where(np.isnan(Zgrid))

#Setup debris mask for use in radiation parameters
if debris == True:
    if debris_treatment == 'Boolean':
        debris_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, debris_array)
        debris_m = np.zeros(debris_grid.shape)
        debris_m[np.where(debris_grid > 100)] = 0.
        debris_m[np.where(debris_grid <= 100)] = 1.
    elif debris_treatment == 'Variable Thickness':
        #EDIT!! get melt factor array
        debristhickness_array = np.load(debris_thickness_map)
        subdebrismeltfactors = MeltFactors(debristhickness_array,tt)
        debris_m = np.ones(Zgrid.shape)
        debris_m[nanlocs] = np.nan 
        print('edit!!!!!!!!!!!!!!!!!')
    else:
        print('Invalid debris treatment option')
else:
    debris_m = np.ones(Zgrid.shape)
    debris_m[nanlocs] = np.nan 


#open a writeout file to keep track of the step the model is on
rsl_file = open('writeout.txt','w')
print('number of simulations in total = ' + str(len(MF_p)))
# begin looping through each simulation
while sim<(len(MF_p)-1):
#while sim < 0:
    print('\rRun ' + str(sim+2) + ' started:',)
    rsl_file.write('Run ' + str(sim+2) + ' started:')
    rsl_file.write('\n')
    
    #counter for simulations with tested param combinations
    sim+=1
    #load param combinations for simulation    
    aice = aice_p[sim]
    asnow = asnow_p[sim]
    MF = MF_p[sim]
    
    aice_a = np.ones(Zgrid.shape) * aice_p[sim] * debris_m #multiplying my debris_m makes aice = 0 in debris covered cells
    asnow_a = np.ones(Zgrid.shape) * asnow_p[sim]
    MF_a = np.ones(Zgrid.shape) * MF_p[sim]
    MF_a[nanlocs] = np.nan #NEW LINE!! as above, trying to fix MB for non debris case
    
    for y in years:
        current_year = y
        print('starting year: ' + str(current_year))
        rsl_file.write('starting year: ' + str(current_year))
        rsl_file.write('\n')
        
    
#-------Inputs: Temperature, Precipitation, and geographical information-------------
        
#---------Set file names-----------------------------------------
        # may not be able to change too much due to some dependence on the downscaling script --> check again later
        if BC_T == True:
            File_temp_name = 'Temp' + glacier_ID + '_BC_' + str(current_year) + '.nc'
        else:
            File_temp_name = 'Temp' + glacier_ID + str(current_year) + '.nc'
        
        if BC_P == True:
            File_precip_name = 'Prcp' + glacier_ID + '_BC_' + str(current_year) + '.nc'
        else:
            File_precip_name = 'netSnow' + glacier_ID + str(current_year) + '.nc'
        
        File_PDCSR_name = 'Srad' + glacier_ID + str(current_year) + '.nc'
        File_snow_in = 'snowini' + str(current_year) + '.txt'
        File_CC_in = 'CCini' + str(current_year) + '.txt'
        
        
        File_temp_in = os.path.join(Temp_input_path,File_temp_name)
        File_precip_in = os.path.join(Precip_input_path,File_precip_name)
        File_PDCSR_in = os.path.join(Srad_input_path,File_PDCSR_name)
        

        #if debris == True:
        #    File_glacier_in = glacier_ID + '_deb.txt'
        #else:
            #File_glacier_in = glacier_ID + '.txt'
            #File_glacier_in = glacier_ID + '_deb.txt'
        #print('File_glacier_in = ' + str(File_glacier_in))
            
        #temp_inputfile = 'temp_' + glacier_ID + str(current_year) + '.nc'
        #accumulation_inputfile = 'netSnow_' + glacier_ID + str(current_year) + '.nc'
        #PDCSR_inputfile = 'PDCSR_' + glacier_ID + str(current_year) + '.nc'
        #snowpack_inputfile = 'snowpack_' + glacier_ID + str(current_year) + '.txt'
        #CC_inputfile = 'CC_' + glacier_ID + str(current_year) + '.txt'


#---------------Load netcdf data---------------------------------
        
        inT = Dataset(File_temp_in, "r")
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:]
        #print "Temperature in..."
        sys.stdout.flush()
        
        inP = Dataset(File_precip_in, "r")
        #P_var = 'Net snow'
        if BC_P == True:
            P_var = 'Net snow'
        else:
            P_var = 'Temperature'
        P_array = inP.variables[P_var][:]
        #print "Precipitation in..."
        sys.stdout.flush()
        
        inS = Dataset(File_PDCSR_in, "r")
        S_var = 'Temperature'
        S_array = inS.variables[S_var][:]
        #print "Potential direct clear-sky radiation in..."
        sys.stdout.flush()
        
        if Temp_shift == True:
            print('OG Temp = ' + str(T_array[100,100,100]))
            T_array = T_array + temp_shift_factor
            print('Shifted Temp = ' + str(T_array[100,100,100]))
        else:
            pass
        
        rainlocs = np.where(T_array > R2S)
        snowlocs = np.where(T_array <= R2S)
        
        Rain_array = P_array  * 1
        
        P_array[rainlocs] = 0. # set rain locations to zero in the snow precip array
        Rain_array[snowlocs] = 0.
        
        
        # removed glacier vectors from this section
        
        #aice_a = np.ones(Zgrid.shape) * aice_p[sim] * debris_m #multiplying my debris_m makes aice = 0 in debris covered cells
        #asnow_a = np.ones(Zgrid.shape) * asnow_p[sim]
        #MF_a = np.ones(Zgrid.shape) * MF_p[sim]
        #MF_a[nanlocs] = np.nan #NEW LINE!! as above, trying to fix MB for non debris case
        
        if considering_catchment == True:
            Topo_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, sfc_type)
            Office_locs = np.where(Topo_grid == 1)
            Onice_locs = np.where(Topo_grid == 0)
            Topo = np.empty(Topo_grid.shape)*np.nan
            Topo[Office_locs] = 0 #now all off ice points = 0, which is important for setting icemelt to 0 for off-glacier locations in the mass balance function
            Topo[Onice_locs] = 1
        else:
            Topo = np.ones(Zgrid.shape) #bc all cells are 'on ice' for the kaskawulsh only
        
        if current_year == years[-1]:
            pass
        else:
            CC = cold_content(current_year, years[:-1], P_array, T_array, glacier_ID, DHval, BC_P, Precip_input_path)
            
                   
#-------Set up output files---------------------------------------------
        Melt_output_name = 'NetMelt' + str(current_year) + str(sim)
        Melt_output_path = os.path.join(OUTPUT_PATH,Melt_output_name)
        
        Mb_output_name = 'MB' + str(current_year) + str(sim)
        Mb_output_path = os.path.join(OUTPUT_PATH,Mb_output_name)
        
        Acc_output_name = 'Accumulation' + str(current_year) + str(sim)
        Acc_output_path = os.path.join(OUTPUT_PATH,Acc_output_name)
        
        Refreezing_output_name = 'Refreezing' + str(current_year) + str(sim)
        Refreezing_output_path = os.path.join(OUTPUT_PATH,Refreezing_output_name)
        
        Rain_output_name = 'Rain' + str(current_year) + str(sim)
        Rain_output_path = os.path.join(OUTPUT_PATH,Rain_output_name)
        
        IceMelt_output_name = 'IceMelt' + str(current_year) + str(sim)
        IceMelt_output_path = os.path.join(OUTPUT_PATH,IceMelt_output_name)
        
        SnowMelt_output_name = 'SnowMelt' + str(current_year) + str(sim)
        SnowMelt_output_path = os.path.join(OUTPUT_PATH,SnowMelt_output_name)
        
        #Snowmelt_output_name = 'Snowmelt' + str(current_year)  + str(sim)  
        #Icemelt_output_name = 'Icemelt' + str(current_year)  + str(sim) 
        #Snowpack_output_name = 'Snowpack' + str(current_year) + str(sim) 

        File_sufix = ".nc"        

        # set up the 'leftover list' --> keeps track of the snowpack that did not melt and adds it to next years snowpack
        # same with tracking excess Cold Content through time
        if considering_catchment == True:
            if current_year == years[0]:
                Leftover_list = np.zeros(Zgrid.shape)
                CC_list = np.zeros(Zgrid.shape)
            else:
                File_snow_in = 'snowini_catchment' + str(current_year) + '.txt'
                File_CC_in = 'CCini_catchment' + str(current_year) + '.txt'
                Leftover_list = np.loadtxt(File_snow_in)
                CC_list = np.loadtxt(File_CC_in)
        else:
            Leftover_list = np.loadtxt(File_snow_in)   
            CC_list = np.loadtxt(File_CC_in)   
        
#-------Set up the TEMPORAL lists----------------------------------------        
        units = inT.variables['time'].units
        dates = netCDF4.num2date(inT.variables['time'][:], units)
        
        #date_list = []
        month_list = []
        
       # for d in dates:
        #    date_el = d.replace(hour=0, minute=0, second=0, microsecond=0)
         #   date_list.append(date_el)
        
        #date_list = np.array(date_list)
        
        
        st_day = start_day #called from namelist
        
        # fruit = the current date
        fruit = dates[(st_day*8) - 8]
        sys.stdout.flush()
        
        ###list to preserve timestep melt arrays
        Melthour = np.empty(T_array.shape) #Melthour is NET MELT
        Snowmelt_hour = np.empty(T_array.shape)
        Icemelt_hour = np.empty(T_array.shape)
        MBhour = np.empty(T_array.shape)
        Leftover_hour = np.empty(T_array.shape)
        Acchour = np.empty(T_array.shape)
        Refreezedhour = np.empty(T_array.shape)
        
        #set the timestep
        delta_t = time_step #from namelist --> given in units of seconds
        
        #run the model for each timestep (for each year, for each 'sim')
        while (fruit < (dates[-1] + dt.timedelta(seconds=delta_t))):
            #print(fruit)
            rsl_file.write(str(fruit))
            rsl_file.write('\n') 
            
            month_list.append(fruit.month)
        
        
            #add cold content in late october to snow pack to prevent winter melt events
            if Refreezing == True:
                if fruit.timetuple().tm_yday == 274:
                    CC_list = CC_list + CC/8.
            else:
                pass

            #determine current temporal position
            pos = np.where(dates == fruit)[0][0]
                    
            curT = T_array[pos,:,:]        
            curS = S_array[pos,:,:]
            ###units conversion to m.w.e. for precipitation done in DS step
            curP = P_array[pos,:,:]
            ###update snowpack with acumulation
            Leftover_list = Leftover_list + curP
            

            ###calculate mass balance
            MB, Melt, Leftover_list_out, Icemelt, Snowmelt, CC_out = MB_vectorized_discreteSnI(curT, curP, curS, Leftover_list, asnow_a, aice_a, MF_a, Topo, CC_list) #NEW LINE!! MF changed to MF_a
                
            deltaCC = CC_list - CC_out #calculate the change in CC over one time step
            #NEW LINE BELOW!!
            if considering_catchment == True:
                Icemelt[Office_locs] = 0 #make sure no ice melt is happening off the glacier
            # set a limit for snow melt on the off ice grid cells --> ie. the max snow melt that can occur is = to the depth of the snowpack
            else:
                pass
            
            ###update leftoverlists
            Leftover_list = Leftover_list_out  
            CC_list = CC_out       
            
            ###setup next timestep and save timestep outputs100                                                       
            move = dt.timedelta(seconds=delta_t)
            fruit = fruit + move
            Melthour[pos] = Melt
            Snowmelt_hour[pos] = Snowmelt
            Icemelt_hour[pos] = Icemelt
            MBhour[pos] = MB
            Acchour[pos] = curP
            Leftover_hour[pos] = Leftover_list 
            Refreezedhour[pos] = deltaCC

        ###generate snow content for next year
        if considering_catchment == True:
            np.savetxt('snowini_catchment' + str(current_year+1) + '.txt', Leftover_list)
            np.savetxt('CCini_catchment' + str(current_year+1) + '.txt', CC_list)
        else:
            np.savetxt('snowini' + str(current_year+1) + '.txt', Leftover_list)
            np.savetxt('CCini' + str(current_year+1) + '.txt', CC_list)
        
        ###insert outputs into .nc containers using funktion
        netcdf_container_gen(Melthour, 'Melt', Melt_output_path, File_sufix, ybounds, xbounds, current_year)
        netcdf_container_gen(MBhour, 'MB', Mb_output_path, File_sufix, ybounds, xbounds, current_year)
        netcdf_container_gen(Acchour, 'Accumulation', Acc_output_path, File_sufix, ybounds, xbounds, current_year)
        netcdf_container_gen(Refreezedhour, 'Refreezing', Refreezing_output_path, File_sufix, ybounds, xbounds, current_year)
        netcdf_container_gen(Rain_array, 'Rain', Rain_output_path, File_sufix, ybounds, xbounds, current_year)
        netcdf_container_gen(Icemelt_hour, 'IceMelt', IceMelt_output_path, File_sufix, ybounds, xbounds, current_year)
        netcdf_container_gen(Snowmelt_hour, 'SnowMelt', SnowMelt_output_path, File_sufix, ybounds, xbounds, current_year)

        inT.close()
        inP.close()
        inS.close()

rsl_file.write('ALL RUNS COMPLETE')
rsl_file.close()

print('ALL RUNS COMPLETE!')