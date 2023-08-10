# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:17:17 2023

This script will make corrections to the 2018 DEM, 
generate a smoothed dh/dt map, and use them to generate
DEMs for the catchment from 1979-2022.

These DEMs will then be used to downscaled NARR Temperature
and Precipitation, and calculate potential direct solar radiation
inputs for the mass-balance model.

@author: katierobinson
"""

import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import utm
from dbfread import DBF
import PIL
from PIL import Image
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import shiftedColorMap
from Model_functions_ver4 import model_domain


Zgrid, Xgrid, Ygrid, xbounds, ybounds, Sfc = model_domain(catchment=True)
Ygrid = np.flipud(Ygrid)
nanlocs = np.where(np.isnan(Zgrid))
catchment_tribarray = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Tributaries.txt')


def get_array_from_tif(file):
    image = Image.open(file)
    array = np.array(image)
    #array[np.where(array == -9999)] = np.nan
    
    return array

DEM1977_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/1977DEM_coregistered_regridded_average.tif') 
DEM1977_resampled[nanlocs] = np.nan
np.where(DEM1977_resampled == -9999) # No holes
np.where(np.isfinite(DEM1977_resampled))[0].shape # Same num of non-nan cells as Zgrid
np.where(np.isfinite(Zgrid))[0].shape

DEM2018_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/2018DEM_coregistered_regridded_average.tif')
np.where(DEM2018_resampled == -9999) # Has some (26) holes that need to be filled
#DEM2018_resampled[np.where(DEM2018_resampled == -9999)] = np.nan


# =============================================================================
# FILL HOLES (-9999s) IN 2018 DEM
# =============================================================================

# Make the edge value (1,193) a hole (-9999) so that it gets replaced.
#DEM2018_resampled[1,193] = -9999

def gapfill(gappy_raster):
    gapfree_raster = np.array(gappy_raster)     # make a copy of the original raster
    
    # change holes to nans to avoid skewing the average elev of surrounding cells:
    # For each hole in the raster:
    for i in range (0,len(np.where(gappy_raster == -9999)[0])):
        Xgrid_copy = np.array(Xgrid)
        Ygrid_copy = np.array(Ygrid)
        DEM_copy = np.array(gappy_raster)
        DEM_copy[np.where(gappy_raster == -9999)] = np.nan
        
        row = np.where(gappy_raster == -9999)[0][i]
        col = np.where(gappy_raster == -9999)[1][i]
    
        # Get the coordinates of the hole cell
        x = Xgrid[row,col]
        y = Ygrid[row,col]
        
        # set the corresponding cell in Xgrid/Ygrid copy arrays to NaN
        Xgrid_copy[row,col] = np.nan
        Ygrid_copy[row,col] = np.nan
        
        # Find the distance of every surrounding cell to the NaN cell
        x_dist = Xgrid_copy - x
        y_dist = Ygrid_copy - y
        distance = np.sqrt((x_dist**2)+(y_dist**2))
        distance[row,col] = np.nan
       
        # Get the 4 closest cells (up, down, left, right)
        closest_cells = np.where(distance == np.nanmin(distance)) # should be 4 each time

        # If the closest cells are all NaNs, set find the NEXT closest cells that are non-NaN
        while np.isnan(np.nanmean(DEM_copy[closest_cells])):
            distance[closest_cells] = np.nan
            closest_cells = np.where(distance == np.nanmin(distance))

        # Get the mean elevation of the surrounding cells, and replace the gap in the original raster.    
        new_elev = np.nanmean(DEM_copy[closest_cells])
        gapfree_raster[row,col] = new_elev
        #gapfree_raster[row,col] = 888
        print('NaN replaced with new elev:' + str(new_elev))
        gapfree_raster[nanlocs] = np.nan
    
    return gapfree_raster
  
print('Gap-filling the 2018 DEM')    
DEM2018_gapfilled = gapfill(DEM2018_resampled)
print('number of NaN cells in new raster:',str(np.where(np.isnan(DEM2018_gapfilled))[0].shape))
np.where(np.isfinite(DEM2018_gapfilled))[0].shape # Same num of non-nan cells as Zgrid

# The array DEM2018_gapfilled is ready to be tested with the downscaling.
np.savetxt('DEM2018_gapfilled.txt',DEM2018_gapfilled)

# =============================================================================
# GENERATE PREVIOUS YEAR SURFACES FROM 2018 DEM AND DH/DT MAP
# =============================================================================

# Load smoothed dh/dt map
dhdt = np.load('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/1977-2018_smoothed_dhdt.npy')

def generate_DEM(year):
    dt = year - 2018
    if year > 2018:
        dt = 0
    print(year, dt)
    DEM = DEM2018_gapfilled + (dhdt*dt)
    np.savetxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Zgrids/DEM_KRH_' + str(year) + '.txt',DEM)
    #return DEM
    
#for year in range(1979,2022+1):
#    generate_DEM(year)
    
















# =============================================================================
# REPLACE ANOMALOUS ELEV VALUES IN 2018 AND 1977 DEMS
# =============================================================================
# =============================================================================
# 
# # First, identify which elev vals are "out of place" or anomalous (ie. don't match surrounding terrain)
# # Can use temperature values from a previous downscaling to achieve this
# # (the previous downscaling used the non-corrected DEM, so there are temp values that are out of place)
# 
# meanT_1979 = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Climate/meanT1979_downscaledwithgappyDEM.txt')
# x_anomalies = Xgrid[np.where(meanT_1979>5)] # from Young et al. (2021), mean annual downscaled temp maxes at 2.5
# y_anomalies = np.flipud(Ygrid)[np.where(meanT_1979>5)]
# 
# # also one cold spot:
# np.where(meanT_1979==meanT_1979[np.where(Zgrid<2000)][876]) # 27,238
# coldx = Xgrid[27,238]
# coldy = np.flipud(Ygrid)[27,238]
# 
# # Some plots to help diagnose the problem:
# plt.figure()
# plt.title('2018 DEM',fontsize=14)
# plt.contourf(Xgrid,np.flipud(Ygrid),DEM2018_gapfilled,cmap='PuRd',levels=np.linspace(700,3600,30))
# legend = plt.colorbar()
# plt.xlabel('Easting',fontsize=14)
# plt.ylabel('Northing',fontsize=14)
# plt.axis('equal')
# legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
# #plt.scatter(x_anomalies,y_anomalies,s=200,facecolors='none',edgecolor='k')
# plt.scatter(coldx,coldy,s=200,facecolors='none',edgecolor='k')
# 
# 
# plt.figure()
# plt.scatter(meanT_1979,DEM2018_gapfilled)
# plt.xlabel('Mean Annual Temperature  ($\degree$C)',fontsize=14)
# plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
# plt.tight_layout()
# 
# plt.figure()
# plt.title('Mean Annual Temperature (1979)')
# plt.contourf(Xgrid,np.flipud(Ygrid),meanT_1979,cmap='YlOrRd',levels=np.linspace(-17,10,28))
# legend = plt.colorbar()
# plt.xlabel('Easting',fontsize=14)
# plt.ylabel('Northing',fontsize=14)
# plt.axis('equal')
# legend.ax.set_ylabel('Temperature ($\degree$C)', rotation=270,fontsize=14,labelpad=20)
# plt.scatter(x_anomalies[0],y_anomalies[0],s=200,facecolors='none',edgecolor='k')
# plt.scatter(coldx,coldy,s=200,facecolors='none',edgecolor='k')
# 
# # Check for big differences in the 1977 and 2018 surfaces
# plt.figure()
# plt.title('2018 DEM minus TanDEM-X',fontsize=14)
# plt.contourf(Xgrid,np.flipud(Ygrid),DEM2018_gapfilled-Zgrid,cmap='RdBu',levels=np.linspace(-400,400,17))
# legend = plt.colorbar()
# legend.ax.set_ylabel('Elevation Difference (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
# plt.xlabel('Easting',fontsize=14)
# plt.ylabel('Northing',fontsize=14)
# plt.axis('equal')
# plt.scatter(x_anomalies[:5],y_anomalies[:5],s=300,facecolors='none',edgecolor='k')
# 
# off_trunk_anomalies = np.where(DEM2018_gapfilled-DEM1977_resampled>200) # located both points on the long off-trunk tributary, matches up w the anomalies
# (DEM2018_gapfilled-Zgrid)[off_trunk_anomalies]
# plt.scatter(Xgrid[off_trunk_anomalies],np.flipud(Ygrid)[off_trunk_anomalies],s=200,facecolors='none',edgecolor='k')
# # these points are 1000 and 320 m above where they are in the Zgrid DEM
# 
# a_anomalies = np.where(DEM2018_gapfilled-DEM1977_resampled<-200)
# (DEM2018_gapfilled-Zgrid)[a_anomalies]
# plt.scatter(Xgrid[a_anomalies],np.flipud(Ygrid)[a_anomalies],s=200,facecolors='none',edgecolor='k')
# 
# # LEFT OFF HERE ^ 07-25
# # looking for the locations of elevation anomalies, then replace them with the equivalent val from Zgrid
# # Then run downscaling for just 1979 and see if that fixed the temperature anomalies.
# 
# DEM2018_gapfilled[off_trunk_anomalies] = Zgrid[off_trunk_anomalies]
# # the line above gets rid of the high topogrpahy points on the small trunk tributary, array = 27, 238
# 
# 
# # Next, find out the difference between the elev. in the anomalous gridcells and the surrounding gridcells
# Zgrid[np.where(meanT_1979>5)] - DEM2018_gapfilled[np.where(meanT_1979>5)]
# 
# DEM2018_fixed = np.array((DEM2018_gapfilled))
# anomalousT = []
# elevdiff = []
# for gridcell in range(0,len([np.where(meanT_1979>5)][0][0])):
#     row = [np.where(meanT_1979>5)][0][0][gridcell]
#     col = [np.where(meanT_1979>5)][0][1][gridcell]
#     #print(row,col)
#     current_z = DEM2018_gapfilled[row,col]
#     North_z = DEM2018_gapfilled[row+1,col]
#     South_z = DEM2018_gapfilled[row-1,col]
#     East_z = DEM2018_gapfilled[row,col+1]
#     West_z = DEM2018_gapfilled[row,col-1]
#     mean_surrounding_z = np.nanmean([North_z,South_z,East_z,West_z])
#     #print(current_z-mean_surrounding_z)
#     print(current_z,current_z-mean_surrounding_z,North_z,South_z,East_z,West_z)
#     elevdiff.append(current_z-mean_surrounding_z)
#     
#     current_T = meanT_1979[row,col]
#     anomalousT.append(current_T)
#     North_T = meanT_1979[row+1,col]
#     South_T = meanT_1979[row-1,col]
#     East_T = meanT_1979[row,col+1]
#     West_T = meanT_1979[row,col-1]
#     mean_surrounding_T = np.nanmean([North_T,South_T,East_T,West_T])
#     print(current_T,current_T-mean_surrounding_T,North_T,South_T,East_T,West_T)
#     
#     print(current_T-mean_surrounding_T)
#     DEM2018_fixed[row,col] = mean_surrounding_z
# 
# plt.figure()
# plt.contourf(Xgrid,np.flipud(Ygrid),DEM2018_fixed,cmap='PuRd')
# legend = plt.colorbar()
# 
# plt.figure(figsize=(5,7))
# plt.scatter(anomalousT,elevdiff)
# plt.xlabel('Temperature (of the hot/cold spots)',fontsize=14)
# plt.ylabel('Elev. of hot spot minus mean elev. of surrounding gridcells',fontsize=14)
# 
# 
# plt.figure(figsize=(4,4))
# plt.scatter(Zgrid,DEM2018_gapfilled)
# plt.xlabel('TanDEM-X (2011-2015) elev (m a.s.l.)')
# plt.ylabel('2018 DEM elev (m a.s.l.)')
# plt.plot([0,4000],[0,4000],linestyle='--',linewidth=0.5,c='k')
# 
# plt.figure(figsize=(4,4))
# plt.scatter(Zgrid,DEM1977_resampled)
# plt.xlabel('Zgrid elev (m a.s.l.)')
# plt.ylabel('1977 DEM elev (m a.s.l.)')
# plt.plot([0,4000],[0,4000],linestyle='--',linewidth=0.5,c='k')
# 
# 
# 
# 
# 
# plt.figure()
# plt.subplot(2,2,1)
# plt.contourf(Xgrid,np.flipud(Ygrid),DEM1977_resampled,cmap='RdPu',levels=np.linspace(700,3600,30))
# plt.colorbar()
# plt.subplot(2,2,2)
# plt.contourf(Xgrid,np.flipud(Ygrid),DEM2018_gapfilled,cmap='RdPu',levels=np.linspace(700,3600,30))
# plt.colorbar()
# plt.subplot(2,2,3)
# plt.contourf(Xgrid,np.flipud(Ygrid),DEM2018_gapfilled-DEM1977_resampled,cmap='RdPu')
# plt.colorbar()
# 
# 
# # Check out the slope:
# slope = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Solar/Hock_model_inputs/KW_Slope_1979.txt',skiprows=6)
# 
# plt.figure()
# plt.contourf(Xgrid,np.flipud(Ygrid),slope,levels=np.linspace(0,60,21),cmap='YlGnBu')
# legend = plt.colorbar()
# legend.ax.set_ylabel('Slope ($\circ$)', rotation=270,fontsize=14,labelpad=20)
# plt.axis('equal') 
# plt.title('Slope (' + str(1979) + ')',fontsize=14)
# plt.xlabel('Easting',fontsize=14)
# plt.ylabel('Northing',fontsize=14)
# plt.scatter(x_anomalies[:5],y_anomalies[:5],s=300,facecolors='none',edgecolor='k')
# 
# plt.figure()
# plt.title('1979')
# plt.scatter(meanT_1979[np.where(meanT_1979>5)],slope[np.where(meanT_1979>5)])
# plt.scatter(meanT_1979[27,238],slope[27,238])
# plt.xlabel('Mean Annual Temperature ($\degree$C)',fontsize=14)
# plt.ylabel('Slope ($\degree$)',fontsize=14)
# plt.ylim(0,60)
# 
# slope = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Solar/Hock_model_inputs/KW_Slope_2022.txt',skiprows=6)
# plt.figure()
# plt.title('2018')
# plt.scatter(meanT_1979[np.where(meanT_1979>5)],slope[np.where(meanT_1979>5)])
# plt.scatter(meanT_1979[27,238],slope[27,238])
# plt.xlabel('Mean Annual Temperature ($\degree$C)',fontsize=14)
# plt.ylabel('Slope ($\degree$)',fontsize=14)
# plt.ylim(0,60)
# 
# asp = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Solar/Hock_model_inputs/KW_Aspect_2022.txt',skiprows=6)
# plt.figure()
# plt.title('2018')
# plt.scatter(meanT_1979[np.where(meanT_1979>5)],asp[np.where(meanT_1979>5)])
# plt.scatter(meanT_1979[27,238],asp[27,238])
# plt.xlabel('Mean Annual Temperature ($\degree$C)',fontsize=14)
# plt.ylabel('Aspect ($\degree$)',fontsize=14)
# plt.ylim(0,360)
# 
# asp = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Solar/Hock_model_inputs/KW_Aspect_1979.txt',skiprows=6)
# plt.figure()
# plt.title('1979')
# plt.scatter(meanT_1979[np.where(meanT_1979>5)],asp[np.where(meanT_1979>5)])
# plt.scatter(meanT_1979[27,238],asp[27,238])
# plt.xlabel('Mean Annual Temperature ($\degree$C)',fontsize=14)
# plt.ylabel('Aspect ($\degree$)',fontsize=14)
# plt.ylim(0,360)
# =============================================================================
