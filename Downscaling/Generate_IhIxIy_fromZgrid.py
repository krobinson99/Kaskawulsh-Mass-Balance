# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:51:05 2023

Trying to reproduce Ih.Ix.Iy from Zgrid (inverse of the regridXYsomething function)

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something

# PLOT the final DEMs for 1977, 2018, and dhdt 1977-2018 and 

# Load the catchment grid that the model will be run on:
File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kask_catchment.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,4] 
Iy = glacier[:,5] 
Ih = glacier[:,6]  
sfc_type = glacier[:,8]
ELA = glacier[:,10]     

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
Sfc_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, sfc_type)
Ygrid = np.flipud(Ygrid)
nanlocs = np.where(np.isnan(Zgrid))
offglacier = np.where(Sfc_grid==1)
onglacier = np.where(Sfc_grid==0)

total_gridcells = Zgrid.shape[0]*Zgrid.shape[1]
total_nancells = len(nanlocs[0])
total_domaincells = total_gridcells - total_nancells
print(total_domaincells)
print(Ih.shape)

 # KR_note: the key here will be figuring out how regridXYsomething takes Ih, Ix, Iy and makes Zgrid
 # then do the reverse on the new Zgrids to make them compatible with this code. 
 # Ih, Iy, Iz are the NON-NAN gridcells only. 
 # idea: make np.zeros(Zgrid.shape), replace nanlocs with NaNs, then loop through gridcells and
 # wherever cell == 0, fill with first element of Ih, then repeat for each gridcell w a zero. (try row wise and column wise)
 # if that method works and creates the right Zgrid, then reverse that to get Ih lists from Zgrid_dynamic. 
         

# NEW IDEA!!!! CREATE ZGRID_INDICES USING REGRIDXYZ FUNCTION (REPLACE IH WITH LIST FROM 0 TO LEN(IH))
Ih_indices = np.arange(0,len(Ih))
Zgrid_indices, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih_indices)
Ygrid = np.flipud(Ygrid)

plt.figure(figsize=(8,5))
plt.title('Position of each gridcell in Ih list',fontsize=14)
plt.contourf(Xgrid,Ygrid,Zgrid_indices, cmap = 'YlGnBu', levels = np.linspace(0,42600,285))
plt.axis('equal')
legend = plt.colorbar()
legend.ax.set_ylabel('Index', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
#plt.savefig('Ih_ordering.png')

# Now try to reproduce Ih from Zgrid using the Zgrid_indices
def get_Ihlist_from_Zgrid(Zgrid,Ih=Ih,Zgrid_indices=Zgrid_indices):
    Ih_new = np.zeros(Ih.shape)

    for i in range(0,len(Zgrid)):
        for j in range(0,len(Zgrid[0])):
            if np.isnan(Zgrid[i,j]):
                pass
            else:
                ind = int(Zgrid_indices[i,j])
                Ih_new[ind] = Zgrid[i,j]
                
    return Ih_new

# this works!!!
Ih_new_test = get_Ihlist_from_Zgrid(Zgrid)
diff = Ih - Ih_new_test

plt.figure()
plt.plot(diff)

# next, turn the Zgrid from each year 1979 -- 2022 into an Ih list
# (maybe one column per year in a txt file or csv that can be read into the downscaling script easily?)

# for now - just get the 1977  and 2018 surfaces (from DEMs) 
# one txt file per year (prob need to make a data frame to save as txt?):
            
def get_array_from_tif(file):
    image = Image.open(file)
    array = np.array(image)
    array[np.where(array == -9999)] = np.nan
    
    return array

# 1977:
DEM1977_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/1977DEM_coregistered_regridded_average.tif') 
DEM1977_resampled[nanlocs] = np.nan

Ih_1977 = get_Ihlist_from_Zgrid(DEM1977_resampled)
# test if I can exactly reproduce DEM1977_resampled
Zgrid1977, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih_1977)
# it works! Zgrid1977 is an exact copy of DEM1977_resampled.

#Next: make a txt file from Zgrid1977
#saved in
#np.savetxt('Ih_1977.txt',Ih_1977)

DEM2018_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/2018DEM_coregistered_regridded_average.tif')
DEM2018_resampled[nanlocs] = np.nan

Ih_2018 = get_Ihlist_from_Zgrid(DEM2018_resampled)
Zgrid2018, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih_2018)

#np.savetxt('Ih_2018.txt',Ih_2018)

dhdt = np.load('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/1977-2018_smoothed_dhdt.npy')
plt.figure(figsize=(10,6))
#plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
plt.contourf(Xgrid,np.flipud(Ygrid),dhdt,cmap = 'RdYlBu',levels=np.linspace(-2,2,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
#plt.title('2007-2018 dh/dt')
plt.text(570000,6710000,'Off-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt[offglacier]),2)) + ' m')
plt.text(570000,6705000,'On-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt[onglacier]),2)) + ' m')
plt.tight_layout()

# get text file for every year - then incorporate in the downscaling script
years = np.arange(2022,1978,-1)
for year in years:
    if year >= 2018:
        #print(year)
        Ih = get_Ihlist_from_Zgrid(Zgrid2018)
    else:
        #print(2018-year)
        Zgrid_new = Zgrid2018 - ((2018-year)*dhdt)
        Ih = get_Ihlist_from_Zgrid(Zgrid_new)
        
    fname = 'Ih_' + str(year) + '.txt'
    print(fname)
    #np.savetxt(fname,Ih)




