# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:39:48 2022

Read coordinate info from CoordinateTable.csv and map Rounce (2021) debris map
to the 200 m mass balance model domain

@author: katierobinson
"""

# steps to generate CoordinateTable.csv in ArcMap:
# 1. upload .tiff file containing debris thickness info
# 2. use Arc Toolbox Raster to Point tool
# 3. Open attribute table in ArcMap, add fields Lat / Lon
# 4. Open editor, select "calculate geometry"
# 5. calculate x,y coordinates for each column
# 6. export shapefile as .csv

import numpy as np
import os
import sys
sys.path.insert(1,'D:\Katie\Mass Balance Model\MassBalanceModel_KatiesVersion\RunModel')
from Model_functions_ver4 import regridXY_something
import utm
import matplotlib.pyplot as plt

# READ IN COORDINATE TABLE (coords for Rounce et al. (2021) debris map)
coordinate_file = 'CoordinateTable.csv'
coordstable = np.loadtxt(coordinate_file,delimiter=",",skiprows=1) #namelist!!
debthickness = coordstable[:,1]
lat = coordstable[:,2]
lon = coordstable[:,3]
x_meters = coordstable[:,4]
y_meters = coordstable[:,5]

# GET COORDINATES FOR MODEL DOMAIN
Path2Model = 'D:\Katie\Mass Balance Model\MassBalanceModel_KatiesVersion\RunModel'
File_glacier_in = os.path.join(Path2Model,'Kaskonly_deb.txt')
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,2] 
debris_array = glacier[:,6]

#-------Turn vectors into 3D gridded inputs--------------------

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs = np.where(np.isnan(Zgrid))

#Setup debris mask for use in radiation parameters
debris_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, debris_array)
debris_m = np.zeros(debris_grid.shape)
debris_m[np.where(debris_grid > 100)] = 1
debris_m[np.where(debris_grid <= 100)] = np.nan
debris_m[nanlocs] = np.nan

Ygrid_flipped = np.flipud(Ygrid) # Ygrid is originally upside down (values decreasing south to north)
easting_EY = []
northing_EY = []
for i in range(0,len(debris_m)):
    for j in range(0,len(debris_m[1])):
        if debris_m[i,j] == 1:
            easting_EY.append(Xgrid[i,j])
            northing_EY.append(Ygrid_flipped[i,j])
        else:
            pass
        
# CONVERT EASTING/NORTHING TO LAT/LON using python module utm
easting_DR = []
northing_DR =[]
for i in range(0,len(lat)):
    x = utm.from_latlon(lat[i],lon[i])
    easting_DR.append(x[0])
    northing_DR.append(x[1])
    
# PLOT THE TWO DIFFERENT GRIDS! COMPARE
plt.figure(figsize=(8,5))
plt.scatter(easting_EY,northing_EY,color='red')
plt.scatter(easting_DR,northing_DR,color='dodgerblue')
plt.legend(['EY Debris Cells','DR Debris Cells'],fontsize=12)
plt.xlabel('Easting (m)',fontsize=12)
plt.ylabel('Northing (m)',fontsize=12)
#plt.savefig('EY_DR_debriscomparison.png')

DR_cellarea = 0.1*0.1
EY_cellarea = 0.2*0.2

DR_debrisarea = len(easting_DR)*DR_cellarea
EY_debrisarea = len(easting_EY)*EY_cellarea

print('DR debris area = ' + str(np.round(DR_debrisarea,2)) + ' km2')
print('EY debris area = ' + str(np.round(EY_debrisarea,2)) + ' km2')


#Problem: There are debris cells in DR's projection that are OUTSIDE the model domain from EY
# ie. np.max(easting_EY) - np.max(Xgrid) = 223.4 m
    



