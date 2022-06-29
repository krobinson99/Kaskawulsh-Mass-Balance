# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 13:27:44 2022

tuning script revamped!! Now with spatially distributed ELA :) 

@author: agribbon
"""

# new features to include in the tuning script: 
# 1. tuning script should take outputs from the mass balance model itself, since I will be making lots of
#       edits to the model soon and don;t want to have to also change the tuning script
# WHAT are the tuning criteria:
    #a) geodetic mass balance from EITHER 1977-2007 or 2007-2018
    #b) spatially + temporally distruted ELAs
    
# procedure:
#    1. add option to MBM to run as TUNING or run as regular
#           - this will decide if the model uses specified input parameters OR randomly generates 1000
#           - if tuning is specified in MBM run: save ONLY mass balance outputs as .npy array not NETCDF: to save time! (saves time maybe, but .npy is BIGGER than .nc)
#DONE
#    2. run MBM with tuning = TRUE --> save outputs to a specific folder 
# Running test now
#    3. calculate geodetic mass balance for the calibration period for each of the 1000 simulations
#    4. only retain the sims (files saved by sim number) that pass the geodetic mass balance
#           - save txt file with list of sim numbers that should be thrown out
#           - use bash rm * commands to delete the non-passing files (echo $listval, rm file_{$listval}) (something like that)
#    5. for each time frame that an ELA is dilineated for --> compare each of the passing sims to it
#           - Q: do I need to add up to MB to that point? OR just look at the single timestep
#                - I think I'd need to add everything up to that point since MB is (change in mass)/dt
#                - add up the MB for that hydrological balance year? iestarting Oct 1 (?)
#                - is snowline the same as equilibrium line in the model? YES IF you add the MB timesteps
#                   starting on Oct 1 (beginning of hydrological year.) because then summed MB will only be < 0 where no snowpack exists (ie total abl exceeded total acc)
#                - Q for Gwenn: confirm that I should sum MB starting at hydrological year and not starting at model beginning. 
#    6. add non-passing sims to the rm list before running bash script
#    7. for all the sims that PASS both stages: save the geodetic MB and metrics about how close the ELA was to see how strict the tuning criteria are being
#           - during the testing stages, save the MB and ELAs from ALL 1000 runs, we can use them to decide how strict to make the criteria! (ie sort all results into a bar chart)
#    8. save the params that correspond the to passing sims in a new list: will be used to run the model! 
#    
#
#

#import libraries and packages    
import numpy as np
import os
import sys
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import KWtributaries

# set initial variables for tuning script
path_to_files = 'D:\TuningOutputs\Test_0628'
param_total = 3 #how many parameter combinations were run (usually 1000)

t_start = 1977
t_end = 2022

calibration_start = 2007
calibration_end = 2008

validation_start = 1979
validation_end = 2007

mb_tuning_condition = -0.46
mb_bounds = 0.17

File_glacier_in = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/downscaled2021acc/kaskonly.txt'
Catchment = False

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
nanlocs = np.where(np.isnan(Zgrid))

# calculate geodetic mass balance for EACH SIM for the calibration period and compare to GMB for validation period
allyears = np.linspace(t_start,t_end,(t_end-t_start)+1)
calyears = np.linspace(calibration_start,calibration_end,(calibration_end-calibration_start)+1)
valyears = np.linspace(validation_start,validation_end,(validation_end-validation_start)+1)
sims = np.linspace(0,param_total-1,param_total)

paramtracker = np.zeros((param_total,2)) #[sim,net mb]
paramtracker[:,0] = sims
for sim in sims:
    print(int(sim))
    mb_list = []
    for year in calyears:
        print(int(year))
        mb = np.load(os.path.join(path_to_files,'MB' + str(int(year)) + str(int(sim)) + '.npy'))
        #print(os.path.join(path_to_files,'MB' + str(int(year)) + str(int(sim)) + '.npy'))
        mb_sum = np.nansum(mb, axis = 0)
        mb_sum[nanlocs]
        mb_annual = np.nanmean(mb_sum)

        mb_list.append(mb_annual)
    
    print(mb_list)
    mb_final = np.mean(mb_list)
    paramtracker[int(sim),1] = mb_final
    
# look at the second column in the param tracker array and make a list of the corresponding sims (column 1) that pass/dont pass

        
        
