# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:50:22 2021

@author: agribbon
"""

######Mass balance model veriosn that runs everything from netcdf and limited 
######text inputs. Output directly to netcdf format is implemented. Note: 
######Temperature downscaling and netcdf coversion must be complete, same for
######precipitation variable. Use master text file for utm corrdinates, 
######elevation, surface type, and debris cover, this is regridded in function.
######This version is a blackbox, but nice for multiple runs
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import netCDF4
from pyproj import Proj
from netCDF4 import Dataset
import datetime as dt
from datetime import datetime
from Model_functions_ver4 import T_downscale_funkfest
from Model_functions_ver4 import Precip_2_Snow
from Model_functions_ver4 import MB_vectorized_discreteSnI
from Model_functions_ver4 import MB_vectorized
from Model_functions_ver4 import rainy_day_funk
from Model_functions_ver4 import closest_node
from Model_functions_ver4 import mean_snowpacks_pts
from Model_functions_ver4 import mean_postemps_pts
#from Model_functions_ver4 import cold_content
from Model_functions_ver4 import polyfit_homebrew
from Model_functions_ver4 import regridXY_something
#from Model_functions_ver4 import get_meanSP
from Model_functions_ver4 import Check_out_that_MB
from Model_functions_ver4 import netcdf_container_gen
import sys
import os
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error


##### EDITS FOR TUNING THE MODEL FOR NO ACCUMULATION BIAS CORRECTION###
T_inputs = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Kaskonly_R2S=1'
P_inputs = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'
SR_inputs = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'
R2S=1
BC_T = True
BC_P = False

#Initialize model
sim = 0
out_label = 0

aice_l = []
asnow_l = []
MF_l = []
NETMB_l = []
mb_means_l = []
mb_stds_l = []
mb_medians_l = []
N_list = []

Glacier_id = "kaskonly"

#Load rubix cube
#params = np.loadtxt('final_params.csv')
#aice_p = params[0,:]
#asnow_p = params[1,:]
#MF_p = params[2,:]
###########################Timerange############################################

#run model for ful time period
earlist = [2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]    
#earlist = [2006] 

#calculate mean snowpack for cold content before looping
DHval = get_meanSP(earlist, Glacier_id,R2S,BC_T,BC_P,T_inputs,P_inputs)

################################################################################

#while sim<len(MF_p):
allruns_MB = np.empty((75,218,328)) #dimensions (#of sims, x, y)
while sim < 75:    
    
    print('\rRun ' + str(sim) + ' started:',)
    
    #counter for simulations with tested param combinations
    
    #load param combinations for simulation    
    #aice = aice_p[sim]
    #asnow = asnow_p[sim]
    #MF = MF_p[sim]  
    
    MB_net = np.load('Mb_net.npy')
#
#    #Model parameters, force positive values and aice greater than asnow
#
    #aice = np.random.normal(0.000003396, 0.00000438) #OG distribution
    aice = np.random.normal(0.000000114, 0.00000438) #NEW shifted distribution for Non-BC accumulation
    while aice < 0:
        #aice = np.random.normal(0.000003396, 0.00000438) #OG
        aice = np.random.normal(0.000000114, 0.00000438) #NEW
    #asnow = np.random.normal(0.000001546,0.00000085) #OG
    asnow = np.random.normal(0.00000111,0.00000085) #NEW
    while asnow < 0:
        #asnow = np.random.normal(0.000001546,0.00000085) #OG
        asnow = np.random.normal(0.00000111,0.00000085) #NEW
    #MF = np.random.normal(0.0002707,0.0001632) #OG
    MF = np.random.normal(0.00011429,0.0001632) #NEW
    while MF < 0:
        #MF = np.random.normal(0.0002707,0.0001632) #OG
        MF = np.random.normal(0.00011429,0.0001632)
        
        
    # set arbitrary values for parameters
    #aice =  2.653565027403637e-06
    #asnow = 1.2413668206900794e-06
    #MF = 0.00020652892765251417
        
        
    #setup lists for arrays    
    aice_l.append(aice)
    asnow_l.append(asnow)
    MF_l.append(MF)
    
    NETMB_run = []
    mb_means_run = []
    mb_stds_run = []
    mb_medians_run = [] 

#######################Usefule time manipulation funktions######################
#these functions are homeade to create time vectors

    str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%m/%d')
        
    def date_linspace(start, end, steps):
        delta = (end - start) / steps
        increments = start + (range(0, steps) * np.array(delta))
        return increments



###########################Begin model HERE#####################################
    MB_empty_container = np.empty((len(earlist),218,328))
    for ear in earlist:

        This_year_is = ear

########Inputs: Temperature, Precipitation, and geographical information########


        ###Set file names###

        File_glacier_in = Glacier_id + '_deb.txt'
        
        #File_temp_in = 'Temp' + Glacier_id + str(This_year_is) + '.nc'
        #File_precip_in = 'netSnow' + Glacier_id + str(This_year_is) + '.nc'
        #File_PDCSR_in = 'Srad' + Glacier_id + str(This_year_is) + '.nc'
        File_temp_name = 'Temp' + Glacier_id + '_BC_' + str(This_year_is) + '.nc'
        File_precip_name = 'netSnow' + Glacier_id + str(This_year_is) + '.nc'
        File_PDCSR_name = 'Srad' + Glacier_id + str(This_year_is) + '.nc'
        File_temp_in = os.path.join(T_inputs,File_temp_name)
        File_precip_in = os.path.join(P_inputs,File_precip_name)
        File_PDCSR_in = os.path.join(SR_inputs,File_PDCSR_name)
        
        File_snow_in = 'snowini' + str(This_year_is) + '.txt'
        File_CC_in = 'CCini' + str(This_year_is) + '.txt'
        
        
        ###Load data: netcdf###
        
        inT = Dataset(File_temp_in, "r")
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:]
        #print "Temperature in..."
        sys.stdout.flush()
        
        inP = Dataset(File_precip_in, "r")
        #P_var = 'Net snow'
        P_var = 'Temperature'
        P_array = inP.variables[P_var][:]
        #print "Precipitation in..."
        sys.stdout.flush()
        
        inS = Dataset(File_PDCSR_in, "r")
        #S_var = 'Potential direct clear sky radiation'
        S_var = 'Temperature'
        S_array = inS.variables[S_var][:]
        #print "Potential direct clear-sky radiation in..."
        sys.stdout.flush()
        
        
        ###Load data: vector form###
        File_glacier_in = 'kaskonly_deb.txt'

        glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        Ix = glacier[:,3]
        Iy = glacier[:,4]
        Ih = glacier[:,2]
        debris_array = glacier[:,6]
        
        debris_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, debris_array)
        debris_m = np.zeros(debris_grid.shape)
        debris_m[np.where(debris_grid > 100)] = 0.
        debris_m[np.where(debris_grid <= 100)] = 1.


###########################Get gridded inputs here###############################

        ####need to do some wizardry to clean up the lists and array components 
        Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
        Topo = np.ones(Zgrid.shape)
        nanlocs = np.isnan(Zgrid)
        
        if This_year_is == earlist[-1]:
            pass
        else:
            CC = cold_content(This_year_is, earlist[:-1], P_array, T_array, Glacier_id, DHval,BC_P,P_inputs)
   
        aice_a = np.ones(Zgrid.shape) * aice * debris_m
        asnow_a = np.ones(Zgrid.shape) * asnow * debris_m
        MF_a = np.ones(Zgrid.shape) * MF

###########################SETUP Model domain###################################

        ###Prep time and space resolution for model, and arrays used for
        ###Operations during MB calculation        
        
        ###arrays for dumping carry-over values for snow in  each model, set initial snow 
        ###condition here: use previous years sow profile. Also generate leftover lists for
        ###tracking snowpack through time
 
        CC_list = np.loadtxt(File_CC_in) 
        if ear == 2006:
            Leftover_list = np.zeros(CC_list.shape)
        else:
            Leftover_list = np.loadtxt(File_snow_in)  
                
        ####set start day for model run
        
        st_day = 1
        it_day = 0
        
        ###Setup starting positions for each input
        
        ###...and set up timesteps through time, with days limited array for indexing
        ###perform twice for hourly values and daily precip
        units = inT.variables['time'].units
        dates = netCDF4.num2date(inT.variables['time'][:], units)
        
        date_list = []
        month_list = []
        
        for d in dates:
            date_el = d.replace(hour=0, minute=0, second=0, microsecond=0)
            date_list.append(date_el)
        
        date_list = np.array(date_list)
        
        fruit = dates[(st_day*8) - 8]
        print("Start year:") 
        print(fruit.year)
        sys.stdout.flush()
        
        ###list to preserve timestep melt arrays
        Melthour = np.empty(T_array.shape)
        MBhour = np.empty(T_array.shape)
        Leftover_hour = np.empty(T_array.shape)

        
        ###run for each timestep
        while (fruit < (dates[-1] + dt.timedelta(hours=3))):
            
            month_list.append(fruit.month)
        
            ###determine temporal position
            pos = np.where(dates == fruit)[0][0]
                    
            curT = T_array[pos,:,:]        
            curS = S_array[pos,:,:]
            ###units conversion to m.w.e. for precipitation done in DS step
            curP = P_array[pos,:,:]
            ###update snowpack with acumulation
            Leftover_list = Leftover_list + curP
            
            ###calculate mass balance
            #MB, Melt, Leftover_list_out, Icemelt, Snowmelt, CC_out  = \
             #   MB_vectorized_discreteSnI(curT, curP, curS, Leftover_list, asnow, aice, MF, Topo, CC_list)
                
            MB, Melt, Leftover_list_out, CC_out = \
                MB_vectorized(curT, curP, curS, Leftover_list, asnow_a, aice_a, MF_a, CC_list)
                
            
            ###update leftoverlist
            Leftover_list = Leftover_list_out  
            CC_list = CC_out                  
            
            ###setup next timestep and save timestep outputs100                                                       
            move = dt.timedelta(hours=3)
            fruit = fruit + move
            Melthour[pos] = Melt
            MBhour[pos] = MB
            Leftover_hour[pos] = Leftover_list
            if This_year_is > 2006.:
                MB_net = np.nansum(np.dstack((MB_net, MB)), 2)
           
            
        #generate snow content for next year
        np.savetxt('snowini' + str(This_year_is+1) + '.txt', Leftover_list)
        np.savetxt('CCini' + str(This_year_is+1) + '.txt', CC_list)
        
        #np.save('MBhour_' + str(This_year_is) +'_' + str(sim) + '.npy',MBhour)
        
        N_MB = np.nansum(MBhour, axis = 0)
        N_MB[nanlocs] = np.nan
        MB_empty_container[ear-2006,:,:] = N_MB
                
        NETMB_year = []
        mb_std_year = []
        mb_mean_year = []
        mb_median_year = []
        
        #extract metric from melthour
        for m in range(1,13):
            locs = np.where(np.asarray(month_list) == m)
            
            netMB = np.nansum(MBhour[locs], axis = 0)
            flat_netMB = np.ndarray.flatten(netMB)
            flat_z = np.ndarray.flatten(Zgrid)
            NETMB = np.nansum(netMB)

            ELA = np.nan
            MBmeans, MBstd, MBmedians = Check_out_that_MB(flat_netMB, flat_z, ELA)
            NETMB_year.append(NETMB)
            mb_mean_year.append(MBmeans)
            mb_std_year.append(MBstd)
            mb_median_year.append(MBmedians)
        
        NETMB_run.append(NETMB_year)
        mb_means_run.append(mb_mean_year)
        mb_stds_run.append(mb_std_year)
        mb_medians_run.append(mb_median_year)
        
        
        inT.close()
        inP.close()
        inS.close()

    NETMB_l.append(NETMB_run)
    mb_means_l.append(mb_means_run)
    mb_stds_l.append(mb_stds_run)
    mb_medians_l.append(mb_medians_run)
    
    #Concatenate all years and all runs
    run_mean_MB = np.nanmean(MB_empty_container,axis=0)
    run_mean_MB[nanlocs] = np.nan
    
    allruns_MB[sim,:,:] = run_mean_MB

    if sim % 1 == 0:
        np.save('mb_means' + str(out_label) + '.npy', np.asarray(mb_means_l))        
        np.save('mb_std' + str(out_label) + '.npy', np.asarray(mb_stds_l)) 
        np.save('mb_median' + str(out_label) + '.npy', np.asarray(mb_medians_l))       
        mb_means_l = []
        mb_stds_l = []
        mb_medians_l = []
        
        out_label += 1
        
    nan_locs = np.isnan(curT)
    MB_net[nan_locs] = np.nan
    
    def balanceflux(upstream):
        return (upstream*200.*200.)/1095554179.

    N_res = balanceflux(np.nansum(MB_net))
    N_list.append(N_res)
    
    sim+=1
  
np.save('NETMB_l.npy',NETMB_l)
# CALCULATE Mean MB,ELA for all params and save outputs as npy files for second stage of tuning
np.save('Aice_' +str(sim) +'.npy',aice_l)
np.save('Asnow_' +str(sim) +'.npy',asnow_l)
np.save('Mf_' +str(sim) +'.npy',MF_l)

Overall_MB_per_run = np.nanmean(allruns_MB,axis=1) #shape is now (# sims,y)
Mb_final_list = np.nanmean(Overall_MB_per_run,axis=1) #shape is now (# sims) (ie. one MB value per sim)
np.save('Mb_' +str(sim) +'.npy',Mb_final_list)

### Calculate ELA #####
ELA_l = []
#for i in range(0,len(aice_l)):
#    MB0 = np.where(allruns_MB[i,:,:] == 0)
#    ELAs = Zgrid[MB0]
#    mean_ELA = np.mean(ELAs)
#    ELA_l.append(mean_ELA)

#np.save('ELA_' +str(sim) +'.npy',ELA_l)
for i in range(0,len(aice_l)):
    file_in = 'mb_means' + str(i) + '.npy'
    means = np.load(file_in)
    for j in range(0,len(means)):
        month_ELA = np.nansum(means[j], axis=0)
        total_ELA = np.nansum(month_ELA, axis=0)
        if np.min(total_ELA[7:36]) > 0:
            ELA_l.append(np.nan)
        else:
            ELA_funk = interp1d(total_ELA[7:36],range(7,36))
            ELA_l.append(ELA_funk(0))

np.save('ELA_' + str(sim) +'.npy',ELA_l)    

print('tuning stage 1 is DONE!!!!!')