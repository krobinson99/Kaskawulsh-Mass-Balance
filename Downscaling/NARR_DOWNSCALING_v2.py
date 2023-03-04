# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:52:00 2023

Cleaned up version of NARR_DOWNSCALING.py by Erik Young

@author: katierobinson
"""

######IMPORTANT################
#                             #
# RUN DOWNSCALINGnamelist.py  # 
#    BEFORE running this      #
#          script!            #
#                             #
###############################
        
import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate
import netCDF4
from pyproj import Proj
from netCDF4 import Dataset
import datetime as dt
from Model_functions_ver3 import T_downscale_funkfest
from Model_functions_ver3 import Precip_2_Snow
from Model_functions_ver3 import MB_341_debris
from Model_functions_ver3 import MB_341_office
from Model_functions_ver3 import rainy_day_funk
from Model_functions_ver3 import closest_node
import sys
import os 

#Import parameters from config file
#from DOWNSCALINGnamelist import *
from DOWNSCALINGnamelist import start_year
from DOWNSCALINGnamelist import end_year
from DOWNSCALINGnamelist import glacier_id
from DOWNSCALINGnamelist import glacier_outline
#from DOWNSCALINGnamelist import Topo_param #not needed - remove from config file
from DOWNSCALINGnamelist import Debris
from DOWNSCALINGnamelist import UTM
#from DOWNSCALINGnamelist import ELA
from DOWNSCALINGnamelist import snow_start
from DOWNSCALINGnamelist import save_snow
from DOWNSCALINGnamelist import NARR_subregions
from DOWNSCALINGnamelist import D_T
from DOWNSCALINGnamelist import D_P
from DOWNSCALINGnamelist import delta_t
from DOWNSCALINGnamelist import cfactor
from DOWNSCALINGnamelist import snowfactor
from DOWNSCALINGnamelist import downscaled_outputs
from DOWNSCALINGnamelist import catchment
from DOWNSCALINGnamelist import rawnarr_inputs
from DOWNSCALINGnamelist import solar_in

#set path for downscaled outputs(Temp, Precip, & Radiation) to go into 
OUTPUT_PATH = downscaled_outputs

# Define years to downscale
years = np.arange(start_year,end_year+1)

for year in years:
    print(year)
     
    #load inputs:
    
    File_elev_in = os.path.join(rawnarr_inputs,'kaskhgt.' + str(year) + '.nc')
    File_temp_in = os.path.join(rawnarr_inputs,'kaskair.' + str(year) + '.nc')
    File_precip_in = os.path.join(rawnarr_inputs,'kaskapcp.' + str(year) + '.nc')
    File_CoarseElev_in = 'kaskCE.nc'
    File_glacier_in = glacier_outline + '.txt' 
    glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
    print('Catchment outline loaded')
    
    if snow_start == True:
        print ("Snow initial condition is 0 m.w.e.,")
    else:
        print ("Snow initial condition is carried over from previous year")
        File_snow_in = 'kaskonlysnow' + str(year) + '.txt'
        snow_ini = np.genfromtxt(File_snow_in, delimiter=',')

    
    inE = Dataset(File_elev_in, "r")
    E_var = 'hgt'
    print("Geopotential in...")
    sys.stdout.flush()
    E_array = inE.variables[E_var][:]
    
    inT = Dataset(File_temp_in, "r")
    T_var = 'air'
    print("Temperature in...")
    sys.stdout.flush()
    T_array = inT.variables[T_var][:]
    
    inP = Dataset(File_precip_in, "r")
    P_var = 'apcp'
    print("Precipitation in...")
    sys.stdout.flush()
    P_array = inP.variables[P_var][:]
    
    inCE = Dataset(File_CoarseElev_in, "r")
    CE_var = 'hgt'
    print("Coarse DEM in...")
    sys.stdout.flush()
    #get coarse NARR grid
    CE_array = inCE.variables[CE_var]
    elev = CE_array[0,:,:]
    
    print("Raw NARR inputs loaded....") 
    
    #Load solar data:
    if catchment == False:
        solar_prefix = 'fixed' + 'kaskonly' + '_' #GLright/left/mid, Kaskice, kasktopo #glacier solar files
        solar_suffix = 'DS.txt'
        #EDIT! lines above should be uncommented: doing a test to see if the solar files can match Ih from kaskonly_deb
        #solar_prefix = 'fixed' + 'kask' + '_' #GLright/left/mid, Kaskice, kasktopo #catchment solar files 
        #solar_suffix = 'DS.txt'
    else:
        solar_prefix = 'fixed' + 'kask' + '_' #GLright/left/mid, Kaskice, kasktopo #catchment solar files 
        solar_suffix = 'DS.txt'
        
    #print('Solar inputs loaded')
            
                        #######OUTPUTS#######
    
    #Prep output file and variables
    print('Setting up files for downscaled outputs')
    
    #TEMPERATURE:
    Temp_outfile_name = 'Temp' + glacier_id + str(year) + '.txt'
    Temp_output_path = os.path.join(OUTPUT_PATH,Temp_outfile_name)
    Temp_out = open(Temp_output_path, 'w')
    
    #SOLAR RADIATION:
    Srad_outfile_name = 'Srad' + glacier_id + str(year) + '.txt'
    Srad_output_path = os.path.join(OUTPUT_PATH,Srad_outfile_name)
    SR_out = open(Srad_output_path, 'w')
    
    #ACCUMULATION:
    netSnow_outfile_name = 'netSnow' + glacier_id + str(year) + '.txt'
    netSnow_output_path = os.path.join(OUTPUT_PATH,netSnow_outfile_name)
    Precip_out = open(netSnow_output_path, 'w')
    
    # edited up to here Mar 3 2023: everything below is not edited yet
    
                #######SETUP Model domain#######
    
    ###Prep time and space resolution for model, and arrays used for
    ###Operations during MB calculations
    
    #set up time steps in a day...
    time_delta = delta_t
    t_st = 0
    t_nd = 24
    timeinterv = np.linspace(t_st, t_nd, t_nd/time_delta, endpoint=False)
    
    #...and set up timesteps through time, with days limited array for indexing
        #perform twice for hourly values and daily precip
    units = inE.variables['time'].units
    dates = netCDF4.num2date(inE.variables['time'][:], units)
    dates_p = netCDF4.num2date(inP.variables['time'][:], units)
    date_list = []
    date_list_p = []
    for d in dates:
        date_el = d.replace(hour=0, minute=0, second=0, microsecond=0)
        date_list.append(date_el)
    for d in dates_p:
        date_el = d.replace(hour=0, minute=0, second=0, microsecond=0)
        date_list_p.append(date_el)
    date_list = np.array(date_list)
    date_list_p = np.array(date_list_p)
    
    #get pressure level dimensions for pl variables
    plinterv = inE.variables['level'][:]
    
    #Get interv for solar data, and set t = -1 day solar values for PCT-UTC gap
    solhour_interv = 0.5
    ini_DS = np.zeros(3)
    
    #Extract grid geometry, correcting for dem/shapefile non-overlaps
    if catchment == False:
        IH = glacier[:,2]
        #Ih[Ih==-9999]=np.nan
        inval_loc = np.where(IH == -9999)
        Ih =np.delete(IH, inval_loc)
        Ix = glacier[:,3]
        Iy = glacier[:,4]
        #TOPO = glacier[:,8] #no sfc_type column in kaskonly_deb.txt
        #ELA = glacier[:,10] #no ELA column in kaskonly_deb.txt
        TOPO = np.zeros(Ih.shape) #TOPO is just an array with same shape as Ih, all vals are 0 (since all cells are on-ice)
        ELA = np.zeros(Ih.shape) #ELA is also an array with same shape as Ih, all vals are 0 
    else:
        IH = glacier[:,6]
        #Ih[Ih==-9999]=np.nan
        inval_loc = np.where(IH == -9999)
        Ih = np.delete(IH, inval_loc)
        Ix = glacier[:,4]
        Iy = glacier[:,5]
        TOPO = glacier[:,8] #topo is 0 if on-ice, 1 if off-ice
        ELA = glacier[:,10] #ELA is also an array with same shape as Ih, all vals are 0 

        
    
    #Create ELA array
    ELAs = []
    for j in range(0, len(ELA)):
        #0 means ice
        if TOPO[j] == 0:
            #when greater than ELA, 1 means add snow
            if Ih[j] > ELA[j]:
                ELAs.append(int(1))
            else:
                ELAs.append(int(0))
        else:
            ELAs.append(int(1))
    
        
    #setup factors and params for MB calculations
    cfac = cfactor
    snowfac = snowfactor
    MFi = (0.0062/len(timeinterv)) * cfac 
    MFs = (0.0054/len(timeinterv)) * cfac
    MF = (0.0026/len(timeinterv)) * cfac
    asnow = 0.00000056 * time_delta * cfac 
    aice = 0.00000086 * time_delta * cfac
    
    #arrays for dumping carry-over values for snow in  each model, set initial snow 
    #condition here: use previous years sow profile
    
    if year == years[0]:
        snow_ini_a = np.zeros(Ih.shape) # set snow profile intitial condition to zero for first year
    else:
        if snow_start == True:
            snow_ini_a = np.zeros(Ih.shape)
        else:    
            snow_ini_a = snow_ini
    
    Leftover_list = np.zeros(np.size(Ih)) + snow_ini_a
    Leftover_listSR_arc = np.zeros(np.size(Ih)) + snow_ini_a
    
    #Control for iffy debris cover input, courtesy of esri
    if Topo == 0:
        if Debris == 1:
            debris_mask[:,1], debris_mask[:,5] = zip(*sorted(zip(debris_mask[:,1], debris_mask[:,5])))
            debris_raw = debris_mask[:,5]
            debris_m = np.zeros(debris_raw.shape) + debris_raw
            debris_m[np.where(debris_raw > 100)] = 0.
            debris_m[np.where(debris_raw <= 100)] = 1.
        else:
            debris_m = np.zeros(Ih.shape) + 1.
    else:
        pass
        
    #number of counts per day for inputs (might be irrelevant in netcdf design)
    #dayT = 32
    #dayP = 2
    
    #set start day for model run
    st_day = 1
    it_day = 0
    
    #Get UTM values for each reanalysis grid point to cycle between nodes on larger 
    #modelled glaciers
    lons = inE.variables['lon'][:]
    lats = inE.variables['lat'][:]
    
    myProj = Proj('+proj=utm +zone=' + UTM_zone + ', +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
    UTMy, UTMx = myProj(lons, lats)
    UTMx[UTMx >100000000000]=np.nan
    UTMy[UTMy >100000000000]=np.nan
    #create list of array positions ungridded
    UTMx_list = UTMx.ravel()
    UTMy_list = UTMy.ravel()
    grid_pts = np.stack((UTMx_list, UTMy_list))
    
    ###Setup starting positions for each input
    fruit = dates[(st_day*8) - 8]
    print("Start year:") 
    print(fruit.year)
    sys.stdout.flush()
    
    
    


                      #######MODEL#######

    ###Begin MB calculations, where dates are used to toggle while loop which
    ###iterates through every day in the record, then through every hour in the
    ###record, in order to obtain both daily and hourly (as per timestep) MB
    ###values across the study area. 
    
    
    while (fruit < dates[-1]):
        
        
        ###Get values for current iteration  
                
        #Current position in time of interest          
        cur_day = fruit.day
        cur_month = fruit.month
        cur_year = fruit.year
        today = np.where(date_list == np.datetime64(dt.datetime(cur_year, cur_month, cur_day)))
        today_p = np.where(date_list_p == np.datetime64(dt.datetime(cur_year, cur_month, cur_day)))
    
        #select values for this day 
        selectedE = E_array[today][:][:][:]
        selectedT = T_array[today][:][:][:]
        selectedP = P_array[today_p][:][:][:]
    
        #get daily rad values from arc files
        today_num = it_day + st_day    
        #UNCOMMENT
        curSR_arc = np.genfromtxt(os.path.join(solar_in,solar_prefix + str(today_num) + solar_suffix), delimiter=' ')
        
        
        for i in range(0, len(timeinterv)):
            
            ###Get values for current hour, indexing is yucky, but straighforward
            ###Note that units for NARR precip are mm we, conversion is done here  
            
            curE = selectedE[i][:][:][:]        
            curT = selectedT[i][:][:][:]
            #units conversion
            curP = selectedP[0][:][:]/1000
            
            #all sorts of empty arrays for calculations
            Thour = []
            Phour = []
            SRhour_arc = []
            SRhour_re = []
            Dhour = []
            precipList = []
            
      
            #Apply funkfest funktion to get inversion tracking downscaled T 
            #functions for every reanalysis grid point
            
            xi_list, yi_list, xi_list_inver, yi_list_inver, funclist, \
                funclist_inver, inversion_list, y0func, y0func_inver, Lfunc, \
                    Lfunc_inver = \
                        T_downscale_funkfest(curT, curE, UTMx_list, UTMy_list) 
                
            #Apply precip funk to get multivariate based downscaling of precip values
            if i == 0:
                W1, Pfunk, I0 = rainy_day_funk(elev.ravel()[NARR_subregions], \
                    curP.ravel()[NARR_subregions], UTMx_list[NARR_subregions], \
                        UTMy_list[NARR_subregions]) 
            else:
                pass           
            
            ###For large Kask type grid use if loop to assign reanalysis grid position                             
            #Get T value for every point on glacier
            #use appropriate interp functon and extrap values
            for z in range(0, len(Ih)):
                
                #Get closest NARR grid point for appropriate downscaling T values
                grid_pt = np.asarray(([Iy[z]], [Ix[z]]))
                #print(grid_pt)
                grid_loc = closest_node(grid_pt, grid_pts)
                #print(grid_loc)
            
                #use index to get nearest grid point in u, w notation
                u = int(np.where(UTMx == grid_pts[0][grid_loc])[0])
                w = int(np.where(UTMy == grid_pts[1][grid_loc])[1])
                
                #find downscaled value at every grid point. If outside bounds, 
                #extrapolate, if within bounds... New functionuses interpolated
                #Lapse rate and y0, meaning no extrapolation routine is needed 
                
                #implement bias correction
                if DT == True:
                    dToday = DTval(today_num)
                else:
                    dToday = 0.
                
                                            
                if inversion_list[u][w] == 0:
                    k = interpolate.bisplev(Iy[z], Ix[z], Lfunc) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func)                         
                        
                else:
                    
                    if Ih[z] < yi_list[u][w][0]:
                        k = interpolate.bisplev(Iy[z], Ix[z], Lfunc_inver) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func_inver) 
                    else:
                        k = interpolate.bisplev(Iy[z], Ix[z], Lfunc) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func)  
                        
                K = (k - 273.3) + dToday  
                
                #Snag date for output for lapse rate, not currently implemented         
                D = timeinterv[i]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                    
                
                #get precip values at start of day
                if i == 0:
                    P1 = (I0 + (Pfunk[0] * Ix[z]) + (Pfunk[1] * Iy[z]) + (Pfunk[2] \
                        * Ix[z] * Iy[z]) + (Pfunk[3] * Ix[z]**2) + (Pfunk[4] * \
                            Iy[z]**2) + (Pfunk[5] * Ih[z]))*W1 
                    
                    #account for negative precip values in areas where statistical
                    #model predicts less than zero value on mx + b regression curve
                    if P1 < 0:
                        P1 = 0.
                    else:
                        pass
                    
                    P2 = curP[u][w]*(1-W1)
                    Pdval = P1 + P2
                    ##vals: Z = Ih[z], Y = Iy[z], X = Ix[z]
                    S = Precip_2_Snow(Pdval, K, snowfac)
                else:
                    S = 0.
                    Pdval = 0.
                
                #if statement to calculate the Precip bias correction factor:
                if DP == True:
                    deltaC = DPval(Ih[z])
                else:
                    deltaC = 1
                
                #Calculate SR/GR enhancement for both reanalysis dataset
                #add sum function when using hourly arc SR. NOTE: code exists for
                #averaged SR scheme in previous verison of this model
                #UNCOMMENT
                SR_arc = curSR_arc[z , :]
                #UNCOMMENT
                SR = SR_arc[int(timeinterv[i]/solhour_interv)]
                
                #reset late solar values
                #UNCOMMENT
                ini_DS = SR_arc[-4:-1]
                
                
                #Update lists  
                Thour.append(K)
                Dhour.append(D)
                Phour.append(S * deltaC)
                #UNCOMMENT
                SRhour_arc.append(SR)
                precipList.append(Pdval)
                
            
            ###Calculate and append hourly MB
            
            #get snow cover from previous timestep, but preserve snow above ELA        
            Pmelt = np.asarray(Phour) + np.asarray(Leftover_list) + np.asarray(ELAs)
            Pmelt[np.where(np.asarray(ELAs) == 1.)] = 1.
            Pmelt_SRarc = np.asarray(Phour) + np.asarray(Leftover_listSR_arc) + np.asarray(ELAs)
            Pmelt_SRarc[np.where(np.asarray(ELAs) == 1.)] = 1.
            
            #get Srad updated to units per hour
            #UNCOMMENT
            SRhour_ARC = np.asarray(SRhour_arc)/solhour_interv
    
            if TOPO[z] == 0:
                MBhour, MBhour_SRarc, Leftover_list, \
                    Leftover_listSR_arc, Melt_list, \
                            Melt_listSR_arc = \
                        MB_341_debris(Thour, Phour, Pmelt, Pmelt_SRarc, \
                            SRhour_ARC, MFi, MFs, asnow, aice, MF, debris_m)
            else:
                MBhour, MBhour_SRarc, Leftover_list, \
                    Leftover_listSR_arc, Melt_list, \
                        Melt_listSR_arc = \
                    MB_341_office(Thour, Phour, Pmelt, Pmelt_SRarc, SRhour_ARC, \
                        MFi, MFs, asnow, aice, MF)      
            
            #Get hourly MB saved
            
            MBwrite = Thour
            for im in MBwrite:
                Temp_out.write("%s," %im) #prepare for next row/day by next line 
            Temp_out.write("\n")
            

            #UNCOMMENT
            MBwrite = SRhour_ARC
            for im in MBwrite:
                SR_out.write("%s," %im) #prepare for next row/day by next line 
            SR_out.write("\n")
            
            
            MBwrite = Phour
            for im in MBwrite:
                Precip_out.write("%s," %im) #prepare for next row/day by next line 
            Precip_out.write("\n")
    
        
        #Move to next day in netcdf files
        print(fruit)
        sys.stdout.flush()
        
        move = dt.timedelta(days=1)
        fruit = fruit + move
        it_day = it_day + 1
        
    
    

    Temp_out.close()
    SR_out.close()
    Precip_out.close()
    inT.close()
    inE.close()
    inP.close()
    
    if save_snow == True:
        np.savetxt(glacier_id + 'snow' + str(current_year+1) + '.txt', Leftover_listSR_arc)
        print('Snow profile is saved as:' + glacier_id + 'snow' + str(current_year+1) + '.txt')
        
    print(glacier_id + 'Complete')
    

    
    
    
    
    
    
    
    