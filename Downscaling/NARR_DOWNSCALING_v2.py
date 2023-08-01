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
# note to self: check all KR_note comments before commiting to git.
        
import numpy as np
#from scipy.interpolate import interp1d
from scipy import interpolate
import netCDF4
from pyproj import Proj
from netCDF4 import Dataset

from Model_functions_ver3 import T_downscale_funkfest
from Model_functions_ver3 import rainy_day_funk
from Model_functions_ver3 import closest_node
import sys
import os 
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import write_config_file

#Import parameters from config file
from DOWNSCALINGnamelist import start_year, end_year
from DOWNSCALINGnamelist import Glacier_ID
from DOWNSCALINGnamelist import glacier_outline
from DOWNSCALINGnamelist import static_surface
from DOWNSCALINGnamelist import UTM
from DOWNSCALINGnamelist import NARR_subregions
from DOWNSCALINGnamelist import normalized_XYZ
from DOWNSCALINGnamelist import delta_t
from DOWNSCALINGnamelist import OUTPUT_PATH
from DOWNSCALINGnamelist import catchment
from DOWNSCALINGnamelist import Climate_inputs, Coarse_DEM_input, Easting_grid, Northing_grid

# Save configuration file for this run to output directory:
write_config_file(OUTPUT_PATH,"DOWNSCALINGnamelist.py")

# Load coarse NARR grid:
# =============================================================================
NARR_DEM = Dataset(Coarse_DEM_input, "r")
coarse_elev =  NARR_DEM.variables['hgt'][0,:,:] # surface geopotential height (time invariant), units = m a.s.l.

lons = NARR_DEM.variables['lon'][:]
lats = NARR_DEM.variables['lat'][:]
sys.stdout.flush()

Projection = Proj('+proj=utm +zone=' + UTM + ', +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
UTMy, UTMx = Projection(lons, lats)  # converts lat/lons of coarse NARR grid to easting, northing on WGS84 projection.

#create list of array positions ungridded
UTMx_list = UTMx.ravel()
UTMy_list = UTMy.ravel()
grid_pts = np.stack((UTMx_list, UTMy_list)) #KR_note: not needed, can be replaced with UTM_xlist/ylist further down in code
print("Coarse NARR grid loaded")
# =============================================================================


# Load grid to downscale onto:
# =============================================================================
Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
print('Model grid loaded')
# =============================================================================

# Load model domain:
File_glacier_in = glacier_outline + '.txt' 
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')

#Extract grid geometry, correcting for dem/shapefile non-overlaps
# KR_note: delete all this stuff once downscaling is set up for arrays
if static_surface == True:
    if catchment == False:
        IH = glacier[:,2]
        inval_loc = np.where(IH == -9999)
        Ih =np.delete(IH, inval_loc)
        Ix = glacier[:,3]
        Iy = glacier[:,4]
    else:
        IH = glacier[:,6]
        inval_loc = np.where(IH == -9999)
        Ih = np.delete(IH, inval_loc)
        Ix = glacier[:,4]
        Iy = glacier[:,5]
    print('Static Surface Loaded')
else:
    pass

# Setting-up time scale and timesteps
# =============================================================================
years = np.arange(start_year,end_year+1)
time_steps = np.arange(0,24,delta_t)
print('Downscaling years:',list(years))
# =============================================================================


# Begin looping through years: 
# =============================================================================
for year in years:
    print(year)
    # KR_note: delete once downscaling is set yp for arrays, but keep option to have static vs dynamic surface
    if static_surface == False:
        Ih_file = 'Ih_' + str(year) + '.txt'
        Ih = np.loadtxt(Ih_file)
        Ix = glacier[:,4]
        Iy = glacier[:,5]
        print('Dynamic Surface Loaded')
    else:
        pass
    
    # Load DEM: 
    Zgrid = np.loadtxt('Ih_' + str(year) + '.txt') #KR_note: replace this with 2d DEM array for each year: use consistent naming convention
    
    File_elev_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_hgt.' + str(year) + '.nc')
    File_temp_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_air.' + str(year) + '.nc')
    File_precip_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_apcp.' + str(year) + '.nc')
    
    # Load NARR inputs:
    inE = Dataset(File_elev_in, "r")
    sys.stdout.flush()
    E_array = inE.variables['hgt'][:]
    
    inT = Dataset(File_temp_in, "r")
    sys.stdout.flush()
    T_array = inT.variables['air'][:]
    
    inP = Dataset(File_precip_in, "r")
    sys.stdout.flush()
    P_array = inP.variables['apcp'][:]
    print("NARR Geopotential, Temperatue, Precipitation inputs loaded....") 
                   
                        #######OUTPUTS#######
    
    #Prep output file and variables
    print('Setting up files for downscaled outputs')
    
    #TEMPERATURE:
    Temp_outfile_name = 'Temp_' + str(year) + '.txt'
    Temp_output_path = os.path.join(OUTPUT_PATH,Temp_outfile_name)
    Temp_out = open(Temp_output_path, 'w')
    
    #PRECIPITATION:
    netSnow_outfile_name = 'Precip' + Glacier_ID + str(year) + '.txt'
    netSnow_output_path = os.path.join(OUTPUT_PATH,netSnow_outfile_name)
    Precip_out = open(netSnow_output_path, 'w')
    
    
    units = inE.variables['time'].units
    dt_daily = netCDF4.num2date(inP.variables['time'][:], units) # 365 timesteps (daily for 1 year)
    
    dt_3hourly = netCDF4.num2date(inE.variables['time'][:], units)   # 2920 timesteps (3 hourly for 1 year)
    date_list = []
    for timestep in dt_3hourly:
        d = timestep.replace(hour=0, minute=0, second=0, microsecond=0) # sets hour to zero for all dates (still 2920 timesteps)
        date_list.append(d)
    
    dt_3hourly_dayonly = np.array(date_list)
    
                      #######MODEL#######

    ###Begin MB calculations, where dates are used to toggle while loop which
    ###iterates through every day in the record, then through every hour in the
    ###record, in order to obtain both daily and hourly (as per timestep) MB
    ###values across the study area. #KR_note: rephrase this after edits are done. 
    
    for date in dt_daily:
        print(date)
       
        ##Get values for current iteration         
        hourly_indices = np.where(dt_3hourly_dayonly == date)
        daily_index = np.where(dt_daily == date)
        DOY = daily_index[0][0]+1  # Day Of Year index (e.g. Jan 1 = day 1)

        dailyE = E_array[hourly_indices]              #Selects all pressure levels and gridcells for current day. Shape=(8,29,6,6)
        dailyT = T_array[hourly_indices]
        dailyP = P_array[daily_index]
    
        for i in range(0, len(time_steps)): #Looping through each timestep (8 per day)
            ###Get values for current hour, indexing is yucky, but straighforward
            ###Note that units for NARR precip are mm we, conversion is done here  
            
            hourlyE = dailyE[i]      
            hourlyT = dailyT[i]
            hourlyP = dailyP[0]/1000 # Units for NARR precip are mm we, conversion to meters done here.  
            
            # Set up lists to contain downscaled values of T/P/SR
            Thour = []
            Phour = []
            SRhour_arc = []

            #KR_note: next steps:
            # identify exactly where T, P, and SR are downscaled, forget everything else!
            # starting with T.
            
            ########################### GET COEFFICIENTS FOR DOWNSCALING T&P ###########################
            # Young et al. (2021) Original code below
            # Apply T_downscale_funkfest function to get inversion tracking downscaled T functions for every reanalysis grid point
            xi_list, yi_list, xi_list_inver, yi_list_inver, funclist, funclist_inver, inversion_list, y0func, y0func_inver, Lfunc, Lfunc_inver = T_downscale_funkfest(hourlyT, hourlyE, UTMx_list, UTMy_list) 
               
            #Apply rainy_day_funk function to get multivariate based downscaling of precip values:
            # This function generates the b0 through b6 and r_beta squared coefficients used in the precip downscaling.
            # See Young et al. (2021) supplementary material for more info.
            
            if i == 0: #Only need to do once per day
                r_beta2, b_coeffs, b0 = rainy_day_funk(coarse_elev.ravel()[NARR_subregions], hourlyP.ravel()[NARR_subregions], UTMx_list[NARR_subregions], UTMy_list[NARR_subregions]) 
                if normalized_XYZ == True:
                    Xmax = np.max(UTMx_list[NARR_subregions])
                    Ymax = np.max(UTMy_list[NARR_subregions])
                    Zmax = np.max(coarse_elev.ravel()[NARR_subregions])
                    r_beta2, b_coeffs, b0 = rainy_day_funk(coarse_elev.ravel()[NARR_subregions]/Zmax, dailyP.ravel()[NARR_subregions], UTMx_list[NARR_subregions]/Xmax, UTMy_list[NARR_subregions]/Ymax)     
                else:
                    pass
            else:
                pass    
            
            #KR_note: move these comments somewhere else? (e.g. where these steps actually happen)
            ###For large Kask type grid use if loop to assign reanalysis grid position                             
            #Get T value for every point on glacier
            #use appropriate interp functon and extrap values
            
            for z in range(0, len(Ih)):
                
                #Get closest NARR grid point for appropriate downscaling T values
                downscaled_cell = np.asarray(([Iy[z]], [Ix[z]]))
                NARR_cell = closest_node(downscaled_cell, grid_pts) #Finds which NARR gridcell (36 total) is closest to the gridcell being downscaled.
            
                #use index to get nearest grid point in u, w notation
                u = int(np.where(UTMx == grid_pts[0][NARR_cell])[0])
                w = int(np.where(UTMy == grid_pts[1][NARR_cell])[1])
                
                #KR_note: edit/move these comments when done.
                #find downscaled value at every grid point. If outside bounds, 
                #extrapolate, if within bounds... New function uses interpolated
                #Lapse rate and y0, meaning no extrapolation routine is needed 
                
                if inversion_list[u][w] == 0:
                    k = interpolate.bisplev(Iy[z], Ix[z], Lfunc) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func)                         
                        
                else:
                    
                    if Ih[z] < yi_list[u][w][0]:
                        k = interpolate.bisplev(Iy[z], Ix[z], Lfunc_inver) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func_inver) 
                    else:
                        k = interpolate.bisplev(Iy[z], Ix[z], Lfunc) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func)  
                        
                K = (k - 273.15)
                
                if i == 0:  #Only need to do once per day
                    Plocal = (b0 + (b_coeffs[0] * Ix[z]) + (b_coeffs[1] * Iy[z]) + (b_coeffs[2] \
                        * Ix[z] * Iy[z]) + (b_coeffs[3] * Ix[z]**2) + (b_coeffs[4] * \
                            Iy[z]**2) + (b_coeffs[5] * Ih[z]))*r_beta2 
                
                    if normalized_XYZ == True:
                        Plocal = (b0 + (b_coeffs[0] * (Ix[z]/Xmax)) + (b_coeffs[1] * (Iy[z]/Ymax)) + (b_coeffs[2] \
                        * (Ix[z]/Xmax) * (Iy[z]/Ymax)) + (b_coeffs[3] * (Ix[z]/Xmax)**2) + (b_coeffs[4] * \
                            (Iy[z]/Ymax)**2) + (b_coeffs[5] * (Ih[z]/Zmax)))*r_beta2 
                        

                    # Correct for negative precip values in areas where statistical model predicts less than zero value on mx + b regression curve
                    if Plocal < 0:
                        Plocal = 0.
                    else:
                        pass
                    
                    Pregional = hourlyP[u][w]*(1-r_beta2)
                    Pdownscaled = Plocal + Pregional
                else:
                    Pdownscaled = 0.

                #Update lists  
                Thour.append(K)
                Phour.append(Pdownscaled)

            

            # Write downscaled variable to the output files:
            MBwrite = Thour
            for im in MBwrite:
                Temp_out.write("%s," %im) #prepare for next row/day by next line 
            Temp_out.write("\n")
                                    
            MBwrite = Phour
            for im in MBwrite:
                Precip_out.write("%s," %im) #prepare for next row/day by next line 
            Precip_out.write("\n")
    

    Temp_out.close()
    Precip_out.close()
    inT.close()
    inE.close()
    inP.close()
       
    print(Glacier_ID + 'Downscaling Complete')
    
    # KR_note: edited up to here Mar 7 2023: everything below is not edited yet