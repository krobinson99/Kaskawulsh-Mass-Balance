# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 12:33:01 2023

A script to look at the coregistered kaskawulsh DEMs
and construct surface elevation from 1979-2022

@author: katierobinson
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import utm
from dbfread import DBF
import PIL
from PIL import Image
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

plt.figure(figsize=(8,5))
plt.title('Outline from kask_catchment.txt',fontsize=14)
plt.contourf(np.flipud(Zgrid), cmap = 'Greys_r', levels = np.linspace(700,3700,16))
#plt.contourf(np.flipud(Zgrid2), cmap = 'bone', levels = np.linspace(700,3700,16))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
#plt.contour(np.flipud(Sfc_grid), colors = 'black', levels = 0, linewidth = 0.2)
plt.tight_layout()

plt.figure(figsize=(8,5))
plt.title('Sfc_type',fontsize=14)
plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
#plt.contourf(np.flipud(Zgrid2), cmap = 'bone', levels = np.linspace(700,3700,16))
legend = plt.colorbar()
legend.ax.set_ylabel('Surface Type ( 0 = ice, 1 = off-glacier)', rotation=270,fontsize=14,labelpad=20)
#plt.contour(np.flipud(Sfc_grid), colors = 'black', levels = 0, linewidth = 0.2)
plt.tight_layout()
#plt.savefig('Sfc_type.png')


# 1977 DEM, 2018 DEM, and 2007-2018 dhdt were regridded to the model domain using GDALWARP -r average:
# load the regridded files here: replace off-glacier terrain with Zgrid (TanDEMx 2012 elevs)
nanlocs = np.where(np.isnan(Zgrid))
offglacier = np.where(Sfc_grid==1)
onglacier = np.where(Sfc_grid==0)


def get_array_from_tif(file):
    image = Image.open(file)
    array = np.array(image)
    array[np.where(array == -9999)] = np.nan
    
    return array

DEM1977_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/1977DEM_coregistered_regridded_average.tif') 
DEM1977_resampled[nanlocs] = np.nan

plt.figure(figsize=(10,6))
plt.contourf(Xgrid,Ygrid,DEM1977_resampled,cmap = 'viridis',levels=np.linspace(600,5000,85))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.title('1977 DEM')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.savefig('1977DEM_resampledd_coregistered.png')

DEM2018_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/2018DEM_coregistered_regridded_average.tif')
DEM2018_resampled[nanlocs] = np.nan

plt.figure(figsize=(10,6))
plt.contourf(Xgrid,Ygrid,DEM2018_resampled,cmap = 'Greys_r',levels=np.linspace(600,5000,45))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.title('2018 DEM (Coregistered)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
#plt.savefig('2018DEM_resampled_coregistered.png')


dhdt1977_2018_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/dhdt_1977-2018_coregistered_regridded_average.tif') 
dhdt1977_2018_resampled[nanlocs] = np.nan

plt.figure(figsize=(10,6))
plt.contourf(Xgrid,Ygrid,dhdt1977_2018_resampled,cmap = 'RdYlBu',levels=np.linspace(-4,4,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('1977-2018 dh/dt (Coregistered)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Mean off-glacier dh/dt = ' + str(np.round(np.nanmean((dhdt1977_2018_resampled)[offglacier]),2)) + ' m a$^{-1}$')
plt.text(570000,6705000,'Mean on-glacier dh/dt = ' + str(np.round(np.nanmean((dhdt1977_2018_resampled)[onglacier]),2)) +' m a$^{-1}$')
plt.savefig('1977-2018_dhdt_resampled_coregistered.png')

print(np.nanmean(dhdt1977_2018_resampled[onglacier]))
# -0.51
# could this be bc of the extra glaciers (non KW) that has a more negative MB?
# test this by isolating KW and calculating the dh/dt
print(np.nanmean(dhdt1977_2018_resampled[offglacier]))
# -0.012 (1.2 cm per yr)
###############################################################################
# Check that off glacier terrain is approx. unchanged b/w 1977 to 2018 DEMs

DEMdifference_2018to1977 = (DEM2018_resampled - DEM1977_resampled)/42

plt.figure(figsize=(16,5))
plt.subplot(1,2,1)
plt.contourf(Xgrid,Ygrid,dhdt1977_2018_resampled,cmap = 'RdYlBu',levels=np.linspace(-4,4,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('A. Coregistered 1977-2018 dh/dt (from E.B.)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off glacier dh/dt = ' + str(np.round(np.nanmean(dhdt1977_2018_resampled[offglacier]),2)) + ' m a$^{-1}$')
plt.text(570000,6705000,'On glacier dh/dt = ' + str(np.round(np.nanmean(dhdt1977_2018_resampled[onglacier]),2)) + ' m a$^{-1}$')
plt.subplot(1,2,2)
plt.contourf(Xgrid,Ygrid,DEMdifference_2018to1977,cmap = 'RdYlBu',levels=np.linspace(-4,4,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('B. 2018 DEM (coregistered from E.B.) - 1977 DEM (coregistered from E.B.)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off glacier dh/dt = ' + str(np.round(np.nanmean(DEMdifference_2018to1977[offglacier]),2)) + ' m a$^{-1}$')
plt.text(570000,6705000,'On glacier dh/dt = ' + str(np.round(np.nanmean(DEMdifference_2018to1977[onglacier]),2)) + ' m a$^{-1}$')
plt.tight_layout()
#plt.savefig('1977-2018_dhdt_2methods.png')


methods_difference = dhdt1977_2018_resampled - DEMdifference_2018to1977

plt.figure(figsize=(10,6))
plt.contourf(Xgrid,Ygrid,methods_difference,cmap = 'RdYlBu',levels=np.linspace(-0.75,0.75,151))
#plt.pcolormesh(Xgrid,Ygrid,methods_difference,cmap = 'RdYlBu',vmin=-0.75,vmax=0.75)
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Difference between dh/dt maps (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('A minus B')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off glacier difference = ' + str(np.round(np.nanmean(methods_difference[offglacier]),2)) + ' m a$^{-1}$')
plt.text(570000,6705000,'On glacier difference = ' + str(np.round(np.nanmean(methods_difference[onglacier]),2)) + ' m a$^{-1}$')
plt.tight_layout()
#plt.savefig('differencebw_dhdt_methods.png')


plt.figure(figsize=(10,6))
plt.contourf(Xgrid,Ygrid,DEMdifference_2018to1977*42,cmap = 'RdYlBu',levels=np.linspace(-120,120,41))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m)', rotation=270,fontsize=14,labelpad=20)
plt.title('Total elevation change 1977-2018')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off-glacier avg change = ' + str(np.round(np.nanmean((DEMdifference_2018to1977*42)[offglacier]),2)) + ' m')
plt.text(570000,6705000,'On-glacier avg change = ' + str(np.round(np.nanmean((DEMdifference_2018to1977*42)[onglacier]),2)) + ' m')
plt.tight_layout()
#plt.savefig('Total_dh_1977-2018.png')

dhdt0718_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/2007/dhdt_2007-2018_regridded_average.tif') #has No NaNs - completely gap-filled
#gapsdhdt = np.where(np.isnan(dhdt0718_resampled))
#dhdt0718_resampled[gapsdhdt] = 0 # need to fill with something else..
dhdt0718_resampled[nanlocs] = np.nan

plt.figure(figsize=(10,6))
#plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
plt.contourf(Xgrid,Ygrid,dhdt0718_resampled,cmap = 'RdYlBu',levels=np.linspace(-3,3,21))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007-2018 dh/dt')
plt.text(570000,6710000,'Off-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[offglacier]),2)) + ' m')
plt.text(570000,6705000,'On-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[onglacier]),2)) + ' m')
plt.tight_layout()
plt.savefig('2017-2018dhdt_resampled.png')

dhdt0718_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/2007/dhdt_2007-2018_regridded_average.tif') #has No NaNs - completely gap-filled
gapsdhdt = np.where(np.isnan(dhdt0718_resampled))
dhdt0718_resampled[gapsdhdt] = 0 # need to fill with something else..
dhdt0718_resampled[nanlocs] = np.nan

plt.figure(figsize=(10,6))
#plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
plt.contourf(Xgrid,Ygrid,dhdt0718_resampled,cmap = 'RdYlBu',levels=np.linspace(-3,3,21))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007-2018 dh/dt')
plt.text(570000,6710000,'Off-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[offglacier]),2)) + ' m')
plt.text(570000,6705000,'On-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[onglacier]),2)) + ' m')
plt.tight_layout()


###############################################################################
###############################################################################
# construct the surface elevation for each year from 1979 to 2022

# trust the 2018 DEM and the 1977-2018 dh/dt (most recently sent by Etienne)
#subtract the dh/dt map from 2018 to get 2017 and so on
def dynamic_surface_linear():
    years = np.arange(2022,1978,-1)
    Zgrid_dynamic = np.empty((len(years),Zgrid.shape[0],Zgrid.shape[1]))
    
    yy = 1979
    for i in range (len(years)):
        Zgrid_dynamic[i,0,0] = yy
        yy+=1
    
    i = -1
    for y in years:
        #print(y)
        #print(i)
        if y >= 2018:
            Zgrid_dynamic[i,:,:] = DEM2018_resampled
        else:
            Zgrid_dynamic[i,:,:] = Zgrid_dynamic[i+1,:,:] - dhdt1977_2018_resampled
        i += -1
        
    for i in range(0,len(years)):
        Zgrid_dynamic[i,:,:][offglacier] = Zgrid[offglacier]
        
    return Zgrid_dynamic

Zgrid_dynamic_linear1977to2018 = dynamic_surface_linear()
  
#np.save('Zgrid_dynamic_linear1977to2018.npy',Zgrid_dynamic_linear1977to2018)

def dynamic_surface_2phases():
    years = np.arange(2022,1978,-1)
    Zgrid_dynamic = np.empty((len(years),Zgrid.shape[0],Zgrid.shape[1]))
    
    yy = 1979
    for i in range (len(years)):
        Zgrid_dynamic[i,0,0] = yy
        yy+=1
    
    i = -1
    for y in years:
        #print(y)
        #print(i)
        if y >= 2018:
            Zgrid_dynamic[i,:,:] = DEM2018_resampled
        elif y >= 2007 and y < 2018:
            #print(y)
            Zgrid_dynamic[i,:,:] = Zgrid_dynamic[i+1,:,:] - dhdt0718_resampled
        elif y == 2006:
            dhdt1977_2007 = (Zgrid_dynamic[i+1,:,:] - DEM1977_resampled)/((2007-1977)+1)
            Zgrid_dynamic[i,:,:] = Zgrid_dynamic[i+1,:,:] - dhdt1977_2007
        elif y < 2006:
            Zgrid_dynamic[i,:,:] = Zgrid_dynamic[i+1,:,:] - dhdt1977_2007
            
        i += -1
        
    for i in range(0,len(years)):
        Zgrid_dynamic[i,:,:][offglacier] = Zgrid[offglacier]
        
    return Zgrid_dynamic, dhdt1977_2007


Zgrid_dynamic_2phases, dhdt1977_2007_inferred = dynamic_surface_2phases()

plt.figure(figsize=(16,5))
plt.subplot(1,2,1)
plt.contourf(Xgrid,Ygrid,dhdt0718_resampled,cmap = 'RdYlBu',levels=np.linspace(-4,4,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007-2018 dh/dt')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[offglacier]),2)) + ' m a$^{-1}$')
plt.text(570000,6705000,'On glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[onglacier]),2)) + ' m a$^{-1}$')
plt.subplot(1,2,2)
plt.contourf(Xgrid,Ygrid,dhdt1977_2007_inferred,cmap = 'RdYlBu',levels=np.linspace(-4,4,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('1977-2007 dh/dt')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off glacier dh/dt = ' + str(np.round(np.nanmean(dhdt1977_2007_inferred[offglacier]),2)) + ' m a$^{-1}$')
plt.text(570000,6705000,'On glacier dh/dt = ' + str(np.round(np.nanmean(dhdt1977_2007_inferred[onglacier]),2)) + ' m a$^{-1}$')
plt.tight_layout()
#plt.savefig('1977-2007-2018_dhdtmaps.png')

plt.figure(figsize=(16,5))
plt.subplot(1,2,1)
plt.contourf(Xgrid,Ygrid,Zgrid_dynamic_linear1977to2018[0],cmap = 'Purples',levels=np.linspace(700,4000,34))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('1979 Surface (derived from 1977-2018 dh/dt)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.subplot(1,2,2)
plt.contourf(Xgrid,Ygrid,Zgrid_dynamic_2phases[0],cmap = 'Purples',levels=np.linspace(700,4000,34))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('1979 Surface (derived from 2007-2018 dh/dt and 1977-2007 dh/dt)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.tight_layout()
#plt.savefig('1979Zgrid_2methods.png')

difference_1979Zgrids = Zgrid_dynamic_linear1977to2018[0] - Zgrid_dynamic_2phases[0]

plt.figure(figsize=(10,6))
#plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
plt.contourf(Xgrid,Ygrid,difference_1979Zgrids,cmap = 'RdBu',levels=np.linspace(-10,10,101))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('1979 Surface Elevation difference for 1977-2018 dh/dt minus 2007-2018 and 1977-2007 dh/dt')
plt.text(570000,6710000,'Mean off-glacier difference = ' + str(np.round(np.nanmean(difference_1979Zgrids[offglacier]),2)) + ' m')
plt.text(570000,6705000,'Mean on-glacier difference = ' + str(np.round(np.nanmean(difference_1979Zgrids[onglacier]),2)) + ' m')
plt.tight_layout()
plt.savefig('difference_1979_Zgrids.png')
