# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 14:20:48 2023

The purpose of this script is to resample the 1977,2007,and 2018 DEMs


@author: katierobinson
"""

import numpy as np
import os
import sys
import utm
import matplotlib.pyplot as plt
import PIL
from PIL import Image
import cmocean
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something

#open TIF files from Etienne Berthier for 1977 DEM and 1977-2007 dh/dt

DEM1977_tiffile = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/1977DEM/StElias_Map_DEMsub.tif'
DEM1977_image = Image.open(DEM1977_tiffile)
DEM1977 = np.array(DEM1977_image)
DEM1977.shape

plt.figure(figsize=(6,5))
plt.contourf(DEM1977, cmap = 'Greys', levels = [-2,-1])
legend = plt.colorbar(ticks=[-1])
legend.ax.set_ylabel('Raster Value', rotation=270,fontsize=14,labelpad=2)
plt.title('StElias_Map_DEMsub.tif')
plt.savefig('StElias_Map_DEMsub_tif.png')

dhdt_77_07_tiffile = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/1977DEM/StElias_dh_rate_sub.tif'
dhdt_77_07_image = Image.open(dhdt_77_07_tiffile)
dhdt_77_07 = np.array(dhdt_77_07_image)
dhdt_77_07.shape

plt.figure(figsize=(6,5))
plt.contourf(dhdt_77_07, cmap = 'Greys', levels = np.linspace(-2100000000,2100000000,22))
legend = plt.colorbar()
legend.ax.set_ylabel('Raster Value', rotation=270,fontsize=14,labelpad=20)
plt.title('StElias_dh_rate_sub.tif')
plt.savefig('StElias_dh_rate_sub_tif.png')

satdate_tiffile = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/1977DEM/Satellite_DATE_on_StElias_sub.tif'
satdate_image = Image.open(satdate_tiffile)
satdate = np.array(satdate_image)
satdate.shape

mapdate_tiffile = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/1977DEM/Map_DATE_on_StElias_sub.tif'
mapdate_image = Image.open(mapdate_tiffile)
mapdate = np.array(mapdate_image)
mapdate.shape

plt.figure(figsize=(12.5,5))
plt.subplot(1,2,1)
plt.contourf(satdate, cmap = 'Greys', levels = np.linspace(np.min(satdate),np.max(satdate),20))
legend = plt.colorbar()
legend.ax.set_ylabel('Raster Value', rotation=270,fontsize=14,labelpad=20)
plt.title('Satellite_DATE_on_StElias_sub.tif')
plt.subplot(1,2,2)
plt.contourf(mapdate, cmap = 'Greys', levels = np.linspace(np.min(mapdate),np.max(mapdate),11))
legend = plt.colorbar()
legend.ax.set_ylabel('Raster Value', rotation=270,fontsize=14,labelpad=20)
plt.title('Map_DATE_on_StElias_sub.tif')
plt.savefig('SatelliteDate_and_MapDate_tifs.png')
