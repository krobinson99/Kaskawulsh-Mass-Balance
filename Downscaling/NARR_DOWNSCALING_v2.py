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
from Model_functions_ver4 import model_domain

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
units = NARR_DEM.variables['time'].units
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
print('Model coordinates loaded.')
# =============================================================================


# Setting-up time period and timesteps
# =============================================================================
years = np.arange(start_year,end_year+1)
time_steps = np.arange(0,24,delta_t)
print('Downscaling years:',list(years))
# =============================================================================


# Begin looping through years: 
# =============================================================================
for year in years:
    print('Starting downscaling for:',year)

    
    # Load inputs for downscaling (DEM, Temp, Precip, Geopotential Height) 
    # =========================================================================
    #Zgrid = np.loadtxt('DEM_' + str(Glacier_ID) + str(year) + '.txt') #KR_note: replace this with 2d DEM array for each year: use consistent naming convention
    Zgrid, Xgrid2, Ygrid2, xbounds2, ybounds2, Sfc2 = model_domain(catchment=True)
    nanlocs = np.where(np.isnan(Zgrid))
    print("Model DEM loaded.") 
    
    File_elev_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_hgt.' + str(year) + '.nc')
    File_temp_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_air.' + str(year) + '.nc')
    File_precip_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_apcp.' + str(year) + '.nc')
    
    # Geopotential Height:
    inH = Dataset(os.path.join(Climate_inputs,str(Glacier_ID) + '_hgt.' + str(year) + '.nc'), "r")
    sys.stdout.flush()
    H_array = inH.variables['hgt'][:]
    
    # Temperature:
    inT = Dataset(File_temp_in, "r")
    sys.stdout.flush()
    T_array = inT.variables['air'][:]
    
    # Precipitation:
    inP = Dataset(File_precip_in, "r")
    sys.stdout.flush()
    P_array = inP.variables['apcp'][:]
    print("NARR Geopotential, Temperature, Precipitation inputs loaded.") 
    # =========================================================================               
    
    
    # Prep output variables
    # ========================================================================= 
    Downscaled_T = np.empty((T_array.shape[0],Xgrid.shape[0],Xgrid.shape[1]))
    Downscaled_P = np.empty((T_array.shape[0],Xgrid.shape[0],Xgrid.shape[1]))

    # Get timestamps for T (3 hourly) and P (daily)
    dt_daily = netCDF4.num2date(inP.variables['time'][:], units) # 365 timesteps (daily for 1 year)
    dt_3hourly = netCDF4.num2date(inH.variables['time'][:], units)   # 2920 timesteps (3 hourly for 1 year)
    date_list = []
    for timestep in dt_3hourly:
        date_list.append(timestep.replace(hour=0, minute=0, second=0, microsecond=0)) # sets hour to zero for all dates (still 2920 timesteps)
    dt_3hourly_dayonly = np.array(date_list)
    # ========================================================================= 
    
    
    # Start looping through timesteps:
    # ========================================================================= 
    for date in dt_daily:
        print(date)
       
        # Get values for current iteration         
        hourly_indices = np.where(dt_3hourly_dayonly == date)
        daily_index = np.where(dt_daily == date)

        dailyH = H_array[hourly_indices] #Selects all pressure levels and gridcells for current day. Shape=(8,29,6,6)
        dailyT = T_array[hourly_indices] # (8,29,6,6)
        dailyP = P_array[daily_index][0]/1000  # (6,6)  # Convert daily precipitation from mm w.e. to m w.e.
    
        # PRECIP DOWNSCALING: 
        r_beta2, b_coeffs, b0 = rainy_day_funk(coarse_elev.ravel()[NARR_subregions], dailyP.ravel()[NARR_subregions], UTMx_list[NARR_subregions], UTMy_list[NARR_subregions]) 
        Plocal = (b0 + (b_coeffs[0] * Xgrid) + (b_coeffs[1] * Ygrid) + (b_coeffs[2] * (Xgrid * Ygrid)) + (b_coeffs[3] * (Xgrid**2)) + (b_coeffs[4] * (Ygrid**2)) + (b_coeffs[5] * Zgrid))*r_beta2 
        Plocal[np.where(Plocal<0)] = 0       # Correct for negative precip values in areas where statistical model predicts less than zero value on mx + b regression curve
        # KR_note: still need to add P_regional to get total downscaled Precip.
    
        Pregional = np.empty(Plocal.shape)
        for i in range(0, len(time_steps)): #Looping through each timestep (8 per day) 
            
            # Get Geopotential height and temperature vals for current timestep
            hourlyH = dailyH[i]      
            hourlyT = dailyT[i]
            
            # Apply T_downscale_funkfest function to get inversion tracking downscaled T functions for every reanalysis grid point
            xi_list, yi_list, xi_list_inver, yi_list_inver, funclist, funclist_inver, inversion_list, y0func, y0func_inver, Lfunc, Lfunc_inver = T_downscale_funkfest(hourlyT, hourlyH, UTMx_list, UTMy_list) 
            
            
            # Loop over every gridcell in downscaled reanalysis, get downscaled P and T
            for cell in range(0, len(np.where(np.isfinite(Zgrid))[0])):
                x = np.where(np.isfinite(Zgrid))[0][cell]
                y =  np.where(np.isfinite(Zgrid))[1][cell]
                
                #Get closest NARR grid point for appropriate downscaling T values
                downscaled_cell = np.asarray(([Ygrid[x,y]], [Xgrid[x,y]]))
                NARR_cell = closest_node(downscaled_cell, grid_pts) #Finds which NARR gridcell (36 total) is closest to the gridcell being downscaled.
            
                #use index to get nearest grid point in u, w notation
                u = int(np.where(UTMx == grid_pts[0][NARR_cell])[0])
                w = int(np.where(UTMy == grid_pts[1][NARR_cell])[1])
                
                #KR_note: edit/move these comments when done.
                #find downscaled value at every grid point. If outside bounds, 
                #extrapolate, if within bounds... New function uses interpolated
                #Lapse rate and y0, meaning no extrapolation routine is needed 
                
                if inversion_list[u][w] == 0:
                    # Interpolated lapse rate*elev + interpolated sea level temp
                    k = interpolate.bisplev(Ygrid[x,y],Xgrid[x,y], Lfunc) * Zgrid[x,y] + interpolate.bisplev(Ygrid[x,y],Xgrid[x,y], y0func)                         
                        
                else:
                    
                    if Zgrid[x,y] < yi_list[u][w][0]:
                        k = interpolate.bisplev(Ygrid[x,y],Xgrid[x,y], Lfunc_inver) * Zgrid[x,y] + interpolate.bisplev(Ygrid[x,y], Xgrid[x,y], y0func_inver) 
                    else:
                        k = interpolate.bisplev(Ygrid[x,y], Xgrid[x,y], Lfunc) * Zgrid[x,y] + interpolate.bisplev(Ygrid[x,y], Xgrid[x,y], y0func)  
                        
                K = (k - 273.15)
                Downscaled_T[hourly_indices[0][i],x,y] = K
                
                # Now get P_regional:
                if i == 0: #calculate for first timestep only
                    Pregional_cell = dailyP[u][w]*(1-r_beta2)
                    Pregional[x,y] = Pregional_cell
                else:
                    pass
            
            if i == 0:
                Pdownscaled = (Plocal + Pregional)/8
                Pdownscaled[nanlocs] = np.nan
            else:
                pass
            
            Downscaled_P[hourly_indices[0][i],:,:] = Pdownscaled
        
    np.save(os.path.join(OUTPUT_PATH,'DownscaledPtest' + str(year) + '.npy'),Downscaled_P)
    np.save(os.path.join(OUTPUT_PATH,'DownscaledTtest' + str(year) + '.npy'),Downscaled_T)    
    print(Glacier_ID + 'Downscaling Complete')
