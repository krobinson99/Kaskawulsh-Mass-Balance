# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 17:37:41 2023

Post-processing downscaled precip outputs

@author: katierobinson
"""

# Identify downscaled precip with anomalous values
# How to identify? Plot the b0 through b6 precip downscaling coeffs for each year.
# Flag any vals that are more than 1 std dev. from the mean

# Import packages:
import numpy as np
import netCDF4
from pyproj import Proj
from netCDF4 import Dataset
import sys, os
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

# Import parameters from config file
from DOWNSCALINGnamelist import start_year, end_year, Glacier_ID, UTM, NARR_subregions, time_step
from DOWNSCALINGnamelist import Climate_inputs, Coarse_DEM_input, Easting_grid, Northing_grid, Elev_inputs, OUTPUT_PATH, Model_functions

# Import functions for model:
sys.path.insert(1,Model_functions)
from Model_functions_ver4 import write_config_file, save_to_netcdf, closest_node
from Model_functions_ver4 import precip_downscaling, T_downscale_funkfest
from Model_functions_ver4 import model_domain

# Save configuration file for this run to output directory:
write_config_file(OUTPUT_PATH,"DOWNSCALINGnamelist.py")

# Load (time-invariant) coarse NARR grid:
# =============================================================================
NARR_DEM = Dataset(Coarse_DEM_input, "r")
coarse_elev =  NARR_DEM.variables['hgt'][0,:,:] # surface geopotential height (time invariant), units = m a.s.l.

lons = NARR_DEM.variables['lon'][:]
lats = NARR_DEM.variables['lat'][:]
units = NARR_DEM.variables['time'].units
sys.stdout.flush()

Projection = Proj('+proj=utm +zone=' +str(UTM) + ' +ellps=WGS84', preserve_units=False)
UTMx, UTMy = Projection(lons, lats)  # converts lat/lons of coarse NARR grid to easting, northing on WGS84 projection.

#create list of array positions ungridded
UTMx_list = UTMx.ravel()
UTMy_list = UTMy.ravel()
print("Coarse NARR grid loaded")
# =============================================================================


# Load grid to downscale onto:
# =============================================================================
Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
print('Model coordinates loaded.')
# =============================================================================


# Set-up time period and timesteps for downscaling:
# =============================================================================
years = np.arange(start_year,end_year+1)
time_steps = np.arange(0,24,time_step)
print('Downscaling years:',list(years))
# =============================================================================


# Prep output variables
# ========================================================================= 
r2 = []
b_0 = []
b1 = []
b2 = []
b3 = []
b4 = []
b5 = []
b6 = []
Plocal_mean = []
Plocals = []

# Begin looping through years: 
# =============================================================================
for year in years:
    print('Starting downscaling for:',year)

    
    # Load inputs for downscaling (DEM, Temp, Precip, Geopotential Height) 
    # =========================================================================
    Zgrid = np.loadtxt(os.path.join(Elev_inputs,'DEM_' + str(Glacier_ID) + '_' + str(year) + '.txt'))
    nanlocs = np.where(np.isnan(Zgrid))
    print("Model DEM loaded.") 
    
    File_precip_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_apcp.' + str(year) + '.nc')
    
    # Precipitation:
    inP = Dataset(File_precip_in, "r")
    sys.stdout.flush()
    P_array = inP.variables['apcp'][:]
    print("NARR Geopotential, Temperature, Precipitation inputs loaded.") 
    sys.stdout.flush()
    # =========================================================================               

    # Get timestamps for T (3 hourly) and P (daily)
    dt_daily = netCDF4.num2date(inP.variables['time'][:], units) # 365 timesteps (daily for 1 year)

    
    # Start looping through timesteps:
    # ========================================================================= 
    for date in dt_daily:
        print(date)
        sys.stdout.flush()
       
        # Get values for current iteration         
        daily_index = np.where(dt_daily == date)
        dailyP = P_array[daily_index][0]/1000  # (6,6)  # Convert daily precipitation from mm w.e. to m w.e.
    
# =============================================================================
#         PRECIP DOWNSCALING: 
# =============================================================================
        # Get coefficients from NARR Precip
        r_beta2, b_coeffs, b0 = precip_downscaling(coarse_elev.ravel()[NARR_subregions], dailyP.ravel()[NARR_subregions], UTMx_list[NARR_subregions], UTMy_list[NARR_subregions]) 
        Plocal = (b0 + (b_coeffs[0] * Xgrid) + (b_coeffs[1] * Ygrid) + (b_coeffs[2] * (Xgrid * Ygrid)) + (b_coeffs[3] * (Xgrid**2)) + (b_coeffs[4] * (Ygrid**2)) + (b_coeffs[5] * Zgrid))*r_beta2 
        # Correct for negative precip values in areas where statistical model predicts less than zero value on mx + b regression curve
        Plocal[np.where(Plocal<0)] = 0  
        
        r2.append(r_beta2)
        b_0.append(b0)
        b1.append(b_coeffs[0])
        b2.append(b_coeffs[1])
        b3.append(b_coeffs[2])        
        b4.append(b_coeffs[3])        
        b5.append(b_coeffs[4])        
        b6.append(b_coeffs[5])   
        Plocal_mean.append(np.nanmean(Plocal))
        Plocals.append(Plocal)
        
np.savetxt(os.path.join(OUTPUT_PATH,'r2.txt'),np.array(r2))
np.savetxt(os.path.join(OUTPUT_PATH,'b0.txt'),np.array(b_0))
np.savetxt(os.path.join(OUTPUT_PATH,'b1.txt'),np.array(b1))
np.savetxt(os.path.join(OUTPUT_PATH,'b2.txt'),np.array(b2))
np.savetxt(os.path.join(OUTPUT_PATH,'b3.txt'),np.array(b3))
np.savetxt(os.path.join(OUTPUT_PATH,'b4.txt'),np.array(b4))
np.savetxt(os.path.join(OUTPUT_PATH,'b5.txt'),np.array(b5))
np.savetxt(os.path.join(OUTPUT_PATH,'b6.txt'),np.array(b6))
np.savetxt(os.path.join(OUTPUT_PATH,'Plocal.txt'),np.array(Plocal_mean))
        
# Plot mean Plocal from each day:
dates = pd.date_range(start= str(start_year) + '-01-01 00:00:00',end= str(end_year) + '-12-31 21:00:00',freq='24H')

plt.figure(figsize=(9,4))
#plt.title('Mean Daily Downscaled P$_{local}$',fontsize=14)
plt.plot(dates,Plocal_mean,c='mediumblue')
plt.ylim(0,0.1)
plt.ylabel('P$_{local}$ (m w.e. day$^{-1}$)',fontsize=14)
plt.margins(x=0)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.hlines(y=np.mean(Plocal_mean+(2*np.std(Plocal_mean))),xmin=dates[0],xmax=dates[-1],linestyle='--',linewidth=1.3,label='Mean P$_{local}$ + 2$\sigma$',color='k')
plt.hlines(y=np.mean(Plocal_mean+(3*np.std(Plocal_mean))),xmin=dates[0],xmax=dates[-1],linestyle='--',linewidth=1.3,label='Mean P$_{local}$ + 3$\sigma$',color='red')
plt.legend(fontsize=14)
plt.tight_layout()
#plt.savefig(os.path.join(OUTPUT_PATH,'MeanDailyPlocal.pdf'),bbox_inches='tight')
       
# Plot where daily downscaled snow is > mean daily snow + 2 std dev.
anomalous_snow = np.where(Plocal_mean>(np.mean(Plocal_mean)+(2*np.std(Plocal_mean))))[0]

wierd_snowdays = []
wsd_dates = []
for i in np.array([9121, 9823,14686, 13002, 15945]):
    plt.figure()
    plt.contourf(Xgrid,Ygrid,Plocals[i],cmap='BuPu')
    wierd_snowdays.append(Plocals[i])
    plt.axis('equal')
    plt.colorbar()
    plt.title(dates[i])
    wsd_dates.append(dates[i])
    #plt.savefig(os.path.join(OUTPUT_PATH,'AnomalousPrecip' + str(dates[i])[0:10] + '.pdf'),bbox_inches='tight')
 
Zgrid, Xgrid, Ygrid, xbounds, ybounds, Sfc = model_domain(catchment=True)
Catchmentoutline = np.array(Sfc)
Catchmentoutline[np.where(np.isfinite(Sfc))] = 0
Catchmentoutline[np.where(np.isnan(Sfc))] = 1
    
plotletters = ['a) ','b) ','c) ','d) ','e) ']
plt.figure(figsize=(20,8))
for i in range(0,len(wierd_snowdays)):
    print(i) 
    plt.subplot(2,3,i+1)
    if i <= 2:
        upperlim,levelss = 0.06,13
    elif i == 3:
        upperlim,levelss = 5.6,8
    elif i == 4:
        upperlim,levelss = 0.24,9
    plt.contourf(Xgrid,np.flipud(Ygrid),wierd_snowdays[i],cmap='BuPu',levels=np.linspace(0,upperlim,levelss))
    plt.axis('equal')
    legend = plt.colorbar()
    legend.ax.set_ylabel('Accumulation (m w.e. day$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0,colors='k',linewidths=0.3,alpha=0.5)
    plt.contour(Xgrid,np.flipud(Ygrid),Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.title(plotletters[i] + str(wsd_dates[i])[:10],fontsize=16)
    plt.axis('off')

plt.tight_layout()
#plt.savefig(os.path.join(OUTPUT_PATH,'AnomalousPrecip.pdf'),bbox_inches='tight')
 
    
    
# Choose indices to redownscale
redownscale_days = dates[np.where(Plocal_mean>(np.mean(Plocal_mean)+(3*np.std(Plocal_mean))))]

def redownscale_one_day(date,NARR_subregions):
    year = int(str(date)[0:4])
    
    Zgrid = np.loadtxt(os.path.join(Elev_inputs,'DEM_' + str(Glacier_ID) + '_' + str(year) + '.txt'))
    nanlocs = np.where(np.isnan(Zgrid))
    print("Model DEM loaded.") 
    
    File_precip_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_apcp.' + str(year) + '.nc')

    # Precipitation:
    inP = Dataset(File_precip_in, "r")
    sys.stdout.flush()
    P_array = inP.variables['apcp'][:]
    print("NARR Geopotential, Temperature, Precipitation inputs loaded.") 
    sys.stdout.flush()
    # =========================================================================               

    # Get timestamps for T (3 hourly) and P (daily)
    dt_daily = netCDF4.num2date(inP.variables['time'][:], units) # 365 timesteps (daily for 1 year)

    # Start looping through timesteps:
    # ========================================================================= 
    print(date)
    sys.stdout.flush()
   
    # Get values for current iteration         
    daily_index = np.where(dt_daily == date)
    dailyP = P_array[daily_index][0]/1000  # (6,6)  # Convert daily precipitation from mm w.e. to m w.e.

# =============================================================================
#         PRECIP DOWNSCALING: 
# =============================================================================
    # Get coefficients from NARR Precip
    r_beta2, b_coeffs, b0 = precip_downscaling(coarse_elev.ravel()[NARR_subregions], dailyP.ravel()[NARR_subregions], UTMx_list[NARR_subregions], UTMy_list[NARR_subregions]) 
    Plocal = (b0 + (b_coeffs[0] * Xgrid) + (b_coeffs[1] * Ygrid) + (b_coeffs[2] * (Xgrid * Ygrid)) + (b_coeffs[3] * (Xgrid**2)) + (b_coeffs[4] * (Ygrid**2)) + (b_coeffs[5] * Zgrid))*r_beta2 
    # Correct for negative precip values in areas where statistical model predicts less than zero value on mx + b regression curve
    Plocal[np.where(Plocal<0)] = 0  

    Pregional = np.empty(Plocal.shape)
    for cell in range(0, len(np.where(np.isfinite(Zgrid))[0])):
        x = np.where(np.isfinite(Zgrid))[0][cell]
        y =  np.where(np.isfinite(Zgrid))[1][cell]
                
        #Get closest NARR grid point for appropriate downscaling T values
        downscaled_cell = np.asarray(([Xgrid[x,y]], [Ygrid[x,y]]))
        NARR_cell = closest_node(downscaled_cell, UTMx_list, UTMy_list) #Finds which NARR gridcell (36 total) is closest to the gridcell being downscaled.
            
        #use index to get nearest grid point in u, w notation
        u = int(np.where(UTMx == UTMx_list[NARR_cell])[0])
        w = int(np.where(UTMy == UTMy_list[NARR_cell])[1])

        # Now get P_regional (Precip for nearest NARR node):
        Pregional_cell = dailyP[u][w]*(1-r_beta2)
        Pregional[x,y] = Pregional_cell

    Pdownscaled = (Plocal + Pregional) # Split daily precip throughout day
    Pdownscaled[nanlocs] = np.nan
  
    return Pdownscaled

NARRsubregions_original = [8,9,10,14,15,16,20,21,22,26,27,28]

Redo2014 = redownscale_one_day(redownscale_days[0],NARRsubregions_original)
Redo2022 = redownscale_one_day(redownscale_days[1],NARRsubregions_original)

# Compare downscaled P from orographic vs original nodes
P2014_orog = Plocals[np.where(Plocal_mean>(np.mean(Plocal_mean)+(3*np.std(Plocal_mean))))[0][0]]
P2022_orog = Plocals[np.where(Plocal_mean>(np.mean(Plocal_mean)+(3*np.std(Plocal_mean))))[0][1]]

plt.figure()
plt.contourf(Xgrid,Ygrid,Redo2014,cmap='BuPu',levels=np.linspace(0,0.0012,11))
plt.axis('equal')
plt.colorbar()
plt.title(redownscale_days[0])

plt.figure()
plt.contourf(Xgrid,Ygrid,Redo2022,cmap='BuPu',levels=np.linspace(0,0.0012,11))
plt.axis('equal')
plt.colorbar()
plt.title(redownscale_days[1])

plotletters = ['a) ','b) ']
redo_snow = [Redo2014,Redo2022]
redo_dates = [dates[13002],dates[15945]]
fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(10,4))
i = 0
for ax in axes.flat:
    #im = ax.imshow(np.random.random((10,10)), vmin=0, vmax=1)
        im = ax.contourf(Xgrid,np.flipud(Ygrid),redo_snow[i],cmap='BuPu',levels=np.linspace(0,0.0012,11))
        ax.axis('equal')
        #ax.set_title(plotletter[i] + ' ' + str(binsizes[i])+ '\n C$_{mean}$ = ' + str(np.round(np.nanmean(biascorrected_snow[i]),2)),fontsize=16)
        ax.set_title(plotletters[i] + str(redo_dates[i])[:10],fontsize=16)
        #ax.text(578482,6752585,plotletters[i],fontsize=16,fontweight='bold')
        ax.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0,colors='k',linewidths=0.3,alpha=0.5)
        ax.contour(Xgrid,np.flipud(Ygrid),Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
        ax.axis('off')
        print(i)
        i+=1

#fig.subplots_adjust(right=0)
cbar_ax = fig.add_axes([0.25, 0.07, 0.5, 0.05])
legend = fig.colorbar(im, cax=cbar_ax,orientation="horizontal")
legend.formatter.set_powerlimits((0, 0))
legend.formatter.set_useMathText(True)
legend.update_ticks()
legend.ax.set_xlabel('Accumulation (m w.e. day$^{-1}$)', rotation=0,fontsize=16,labelpad=15)
legend.ax.tick_params(labelsize=16)
plt.savefig(os.path.join(OUTPUT_PATH,'Redownscaled_AnomalousPrecip.pdf'),bbox_inches='tight')
plt.show()

# Save new netcdf files where P_downscaled is replaced
def replace_bad_precip(date,INPUT_PATH,NewPrecip):
    year = int(str(date)[0:4])
    
    # Load original netcdf file
    inP = Dataset(os.path.join(INPUT_PATH,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    P_array = inP.variables['Precipitation'][:]
    sys.stdout.flush()
    inP.close()
    
    # Make a copy of the original precip array
    New_P_array = np.array(P_array)
    
    # Get indice for where the bad precip is
    dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='3H')
    index = np.where(dates==date)[0][0]
    
    for i in range(0,8):
        New_P_array[index+i,:,:] = NewPrecip/8
       
    Corrected_P_file = os.path.join(OUTPUT_PATH,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc')
    save_to_netcdf(New_P_array, 'Precipitation', Corrected_P_file, year, Xgrid, Ygrid) 
   
    return New_P_array
 
# The block of code below replaces the precip irregularities with the downscaled
# precip using the Young et al. (2021) grid cells    
# =============================================================================
# inputpath = 'D:\Downscaled_files\KRH'
# new2014_P = replace_bad_precip(redownscale_days[0],inputpath,Redo2014)
# new2022_P = replace_bad_precip(redownscale_days[1],inputpath,Redo2022)
# =============================================================================
    
# The block of code below replaces the precip irregularities with 0 precip.
# =============================================================================
Redo2014_0precip = np.zeros(Redo2014.shape)
Redo2014_0precip[np.where(np.isnan(Redo2014))] = np.nan

Redo2022_0precip = np.zeros(Redo2022.shape)
Redo2022_0precip[np.where(np.isnan(Redo2022))] = np.nan

inputpath = 'D:\Downscaled_files\KRH'
new2014_P = replace_bad_precip(redownscale_days[0],inputpath,Redo2014_0precip)
new2022_P = replace_bad_precip(redownscale_days[1],inputpath,Redo2022_0precip)
# =============================================================================


