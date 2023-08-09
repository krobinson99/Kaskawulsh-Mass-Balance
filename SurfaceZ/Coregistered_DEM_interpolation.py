# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 12:33:01 2023

A script to look at the coregistered kaskawulsh DEMs
and construct surface elevation from 1979-2022

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
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import shiftedColorMap


catchment_tribarray = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Tributaries.txt')

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
plt.contourf(Xgrid,Ygrid,DEM1977_resampled,cmap = 'Greys_r',levels=np.linspace(600,5000,85))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.title('1977 DEM')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
#plt.savefig('1977DEM_resampled_coregistered.png')

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
#plt.savefig('1977-2018_dhdt_resampled_coregistered.png')

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
#dhdt0718_resampled[np.where(dhdt0718_resampled==0)] = np.nan

plt.figure(figsize=(10,6))
#plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
plt.contourf(Xgrid,Ygrid,dhdt0718_resampled,cmap = 'RdYlBu',levels=np.linspace(-3,3,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007-2018 dh/dt')
plt.text(570000,6710000,'Off-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[offglacier]),2)) + ' m')
plt.text(570000,6705000,'On-glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[onglacier]),2)) + ' m')
plt.tight_layout()
#plt.savefig('2007-2018dhdt_resampled.png')

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
#np.save('Zgrid_dynamic_2phases.npy',Zgrid_dynamic_2phases)

plt.figure(figsize=(16,5))
plt.subplot(1,2,2)
plt.contourf(Xgrid,Ygrid,dhdt0718_resampled,cmap = 'RdYlBu',levels=np.linspace(-4,4,81))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007-2018 dh/dt')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[offglacier]),2)) + ' m a$^{-1}$')
plt.text(570000,6705000,'On glacier dh/dt = ' + str(np.round(np.nanmean(dhdt0718_resampled[onglacier]),2)) + ' m a$^{-1}$')
plt.subplot(1,2,1)
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
#plt.savefig('1977-2007-2018_dhdtmaps.pdf')

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
plt.contourf(Xgrid,Ygrid,difference_1979Zgrids,cmap = 'RdBu',levels=np.linspace(-30,30,61))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m)', rotation=270,fontsize=14,labelpad=20)
plt.title('1979 Surface Elevation difference for 1977-2018 dh/dt minus 2007-2018 and 1977-2007 dh/dt')
plt.text(570000,6710000,'Mean off-glacier difference = ' + str(np.round(np.nanmean(difference_1979Zgrids[offglacier]),2)) + ' m')
plt.text(570000,6705000,'Mean on-glacier difference = ' + str(np.round(np.nanmean(difference_1979Zgrids[onglacier]),2)) + ' m')
plt.tight_layout()
plt.savefig('difference_1979_Zgrids.png')

plt.figure(figsize=(16,5))
plt.subplot(1,2,1)
plt.contourf(Xgrid,Ygrid,Zgrid_dynamic_linear1977to2018[28],cmap = 'Purples',levels=np.linspace(700,4000,34))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007 Surface (derived from 1977-2018 dh/dt)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.subplot(1,2,2)
plt.contourf(Xgrid,Ygrid,Zgrid_dynamic_2phases[28],cmap = 'Purples',levels=np.linspace(700,4000,34))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007 Surface (derived from 2007-2018 dh/dt and 1977-2007 dh/dt)')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.tight_layout()
#plt.savefig('2007Zgrid_2methods.png')

difference_2007Zgrids = Zgrid_dynamic_linear1977to2018[28] - Zgrid_dynamic_2phases[28]

plt.figure(figsize=(10,6))
#plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
plt.contourf(Xgrid,Ygrid,difference_2007Zgrids,cmap = 'RdBu',levels=np.linspace(-30,30,61))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m)', rotation=270,fontsize=14,labelpad=20)
plt.title('2007 Surface Elevation difference for 1 dh/dt vs 2 dh/dt methods')
plt.text(570000,6710000,'Mean off-glacier difference = ' + str(np.round(np.nanmean(difference_2007Zgrids[offglacier]),2)) + ' m')
plt.text(570000,6705000,'Mean on-glacier difference = ' + str(np.round(np.nanmean(difference_2007Zgrids[onglacier]),2)) + ' m')
plt.tight_layout()
plt.savefig('difference_2007_Zgrids.png')

#plot area vs elevation
def area_vs_elev_curve(Zgrid,binnum):
    counts, bins = np.histogram(Zgrid[onglacier],bins=binnum,range=(600,4000))
    #with 12 bins the bin width is 300m
    area = counts*(0.2*0.2) #multiplied by area per gridcell in units km^2
    elevs = bins[1:] #bins are the upper bound (corresponds to total area under that elev.)
    return area,elevs

plt.figure(figsize=(5,7))
for i in range(0,len(Zgrid_dynamic_2phases)):  
    area, elev = area_vs_elev_curve(Zgrid_dynamic_2phases[i],24)
    plt.plot(area,elev)
plt.xlabel('Area (km$^2$)')
plt.ylabel('Elevation (m a.s.l.)')

plt.figure(figsize=(4,8))
area_i, elev_i = area_vs_elev_curve(Zgrid_dynamic_2phases[0],24)
plt.plot(area_i,elev_i,label='1979')
area_f, elev_f = area_vs_elev_curve(Zgrid_dynamic_2phases[-1],24)
plt.plot(area_f,elev_f,label='2022')
plt.legend()
plt.xlabel('Area (km$^2$)')
plt.ylabel('Elevation (m a.s.l.)')

area_i, elev_i = area_vs_elev_curve(Zgrid_dynamic_2phases[0],17)
area_f, elev_f = area_vs_elev_curve(Zgrid_dynamic_2phases[-1],17)
x = np.arange(len(elev_i))
width = 0.4
plt.figure(figsize=(6,8))
plt.barh(x+0.2, area_i, width,color='deeppink')
plt.barh(x-0.2, area_f, width,color='cornflowerblue')
plt.legend(['1979','2022'],fontsize=14)
plt.yticks(x,elev_i, fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlabel('Area (km$^2$)',fontsize=14)
plt.tight_layout()
#plt.savefig('Hypsometry_1979vs2022.png')

area_lin, elev_lin = area_vs_elev_curve(Zgrid_dynamic_linear1977to2018[28],17)
area_2p, elev_2p = area_vs_elev_curve(Zgrid_dynamic_2phases[28],17)
x = np.arange(len(elev_lin))
width = 0.4
plt.figure(figsize=(6,8))
plt.barh(x+0.2, area_lin, width,color='gold')
plt.barh(x-0.2, area_2p, width,color='blueviolet')
plt.legend(['2007 Surface\n(1977-2018 dh/dt)','2007 Surface\n(1977-2007 & 2007-2018 dh/dt)'],fontsize=14)
plt.yticks(x,elev_i, fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlabel('Area (km$^2$)',fontsize=14)
plt.tight_layout()
#plt.savefig('Hypsometry_2007.png')

plt.figure(figsize=(10,6))
plt.contourf(Xgrid,Ygrid,Zgrid_dynamic_2phases[-1]-Zgrid_dynamic_2phases[0],cmap = 'RdBu',levels=np.linspace(-120,120,41))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m)', rotation=270,fontsize=14,labelpad=20)
plt.title('Total elevation change 1977-2018')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.text(570000,6710000,'Off-glacier avg change = ' + str(np.round(np.nanmean((Zgrid_dynamic_2phases[-1]-Zgrid_dynamic_2phases[0])[offglacier]),2)) + ' m')
plt.text(570000,6705000,'On-glacier avg change = ' + str(np.round(np.nanmean((Zgrid_dynamic_2phases[-1]-Zgrid_dynamic_2phases[0])[onglacier]),2)) + ' m')
plt.tight_layout()
plt.savefig('Total_dh_1977-2018_COREGISTERED.png')

coregistered_0718 = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/dhdt_2018-2007.tif') #has No NaNs - completely gap-filled

plt.figure(figsize=(7,6))
plt.axis('equal') 
plt.contourf(np.flipud(coregistered_0718), cmap = 'RdYlBu',levels=np.linspace(-3,3,25))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation Change (m a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.title('dh/dt 2007--2018, no fill')


dhdt0718_resampled = get_array_from_tif('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/2007/dhdt_2007-2018_regridded_average.tif') #has No NaNs - completely gap-filled
#gapsdhdt = np.where(np.isnan(dhdt0718_resampled))
#dhdt0718_resampled[gapsdhdt] = 0 # need to fill with something else..
dhdt0718_resampled[nanlocs] = np.nan

elev = []
dhdt = []
for i in range(0,len(dhdt0718_resampled)):
    for j in range(0,len(dhdt0718_resampled[0])):
        if np.isnan(dhdt0718_resampled[i,j]):
            pass
        elif Sfc_grid[i,j] == 1:
            pass
        else:
            elev.append(Zgrid[i,j])
            dhdt.append(dhdt0718_resampled[i,j])
            
plt.figure()
plt.scatter(dhdt,elev)

def dhdt_vs_elev_curve(dhdtmap,Zgrid=Zgrid,binnum=32):
    dhdt = []
    counts, bins = np.histogram(Zgrid[onglacier],bins=binnum,range=(600,3800))
    #with 12 bins the bin width is 300m
    for i in range(0,len(bins)-1):
        bin_l = bins[i]
        bin_u = bins[i+1]
        #print(bin_l,bin_u)
        targetelevs = (Zgrid[onglacier] >= (bin_l)) & (Zgrid[onglacier] <= (bin_u))
        dhdt.append(np.nanmean(dhdtmap[onglacier][targetelevs]))
    
    binsize = bins[1] - bins[0]
    elevs = bins[:-1] + (0.5*binsize) #elevs are the midpoint value of each bin (ie. 150m is the bin for 100-200m)
    return dhdt,elevs

dhdt,elev = dhdt_vs_elev_curve(dhdt0718_resampled)
x = np.arange(len(elev))
width = 0.7
plt.figure(figsize=(5,7))
plt.title('2007-2018 dh/dt',fontsize=14)
plt.barh(x, dhdt, width,color='cornflowerblue')
#plt.barh(x-0.2, area_f, width,color='cornflowerblue')
#plt.legend(['1979','2022'],fontsize=14)
plt.yticks(x,elev, fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
plt.tight_layout()
#plt.savefig('2007-2018dhdt_elev_histogram.png')

dhdt1977_2018_resampled[np.where(np.isnan(dhdt0718_resampled))] = np.nan
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2018_resampled)
x = np.arange(len(elev))
width = 0.7
plt.figure(figsize=(5,7))
plt.title('1977-2018 dh/dt',fontsize=14)
plt.barh(x, dhdt, width,color='cornflowerblue')
#plt.barh(x-0.2, area_f, width,color='cornflowerblue')
#plt.legend(['1979','2022'],fontsize=14)
plt.yticks(x,elev, fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
plt.tight_layout()  
#plt.savefig('1977-2018dhdt_elev_histogram.png')   

dhdt1977_2007_inferred[np.where(np.isnan(dhdt0718_resampled))] = np.nan
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred)
x = np.arange(len(elev))
width = 0.7
plt.figure(figsize=(5,7))
plt.title('1977-2007 dh/dt',fontsize=14)
plt.barh(x, dhdt, width,color='cornflowerblue')
#plt.barh(x-0.2, area_f, width,color='cornflowerblue')
#plt.legend(['1979','2022'],fontsize=14)
plt.yticks(x,elev, fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
plt.tight_layout()  
#plt.savefig('1977-2007dhdt_elev_histogram.png')   

x = np.arange(len(elev))
width = 0.3
plt.figure(figsize=(6,8))
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2018_resampled)
plt.barh(x-0.3, dhdt, width,color='cornflowerblue')
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred)
plt.barh(x, dhdt, width,color='gold')
dhdt,elev = dhdt_vs_elev_curve(dhdt0718_resampled)
plt.barh(x+0.3, dhdt, width,color='deeppink')
plt.legend(['1977-2018 dh/dt','1977-2007 dh/dt','2007-2018 dh/dt'],fontsize=14)
plt.yticks(x,elev, fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
#plt.plot(dhdt,elev)
plt.tight_layout()  
#plt.savefig('dhdt_elev_histogram.png')   

# use np.polyfit to get dh/dt curve
def deg1_interp(dhdt,elev):
    m, b = np.polyfit(np.array(dhdt)[~np.isnan(dhdt)],np.array(elev)[~np.isnan(dhdt)],deg=1)
    x = np.linspace(-3.5,1,31)
    y = m*x + b
    return x,y

def deg2_interp(dhdt,elev):
    a, b, c = np.polyfit(np.array(dhdt)[~np.isnan(dhdt)],np.array(elev)[~np.isnan(dhdt)],deg=2)
    x = np.linspace(-3.5,1,31)
    y = a*(x**2) + b*x + c
    return x,y

def deg3_interp(dhdt,elev):
    a, b, c, d = np.polyfit(np.array(dhdt)[~np.isnan(dhdt)],np.array(elev)[~np.isnan(dhdt)],deg=3)
    x = np.linspace(-3.5,1,31)
    y = a*(x**3) + b*(x**2) + c*x + d
    return x,y

plt.figure(figsize=(12,7))
plt.subplot(1,2,1)
x = np.arange(650,3751,100)
width = 30
#plt.figure(figsize=(6,8))
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2018_resampled)
plt.barh(x-30, dhdt, width,color='cornflowerblue',label='1977-2018 dh/dt')
xl, yl = deg1_interp(dhdt,elev)
plt.plot(xl,yl,color='cornflowerblue',linewidth=3)
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred)
plt.barh(x, dhdt, width,color='gold',label='1977-2007 dh/dt')
xl, yl = deg1_interp(dhdt,elev)
plt.plot(xl,yl,color='gold',linewidth=3)
dhdt,elev = dhdt_vs_elev_curve(dhdt0718_resampled)
plt.barh(x+30, dhdt, width,color='deeppink',label='2007-2018 dh/dt')
xl, yl = deg1_interp(dhdt,elev)
plt.plot(xl,yl,color='deeppink',linewidth=3)
plt.legend(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
#plt.plot(dhdt,elev)
plt.xlim(-3.5,0.5)
plt.ylim(700,3700)
#plt.tight_layout()
  
plt.subplot(1,2,2)
x = np.arange(650,3751,100)
width = 30
#plt.figure(figsize=(6,8))
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2018_resampled)
plt.barh(x-30, dhdt, width,color='cornflowerblue',label='1977-2018 dh/dt')
xl, yl = deg2_interp(dhdt,elev)
plt.plot(xl,yl,color='cornflowerblue',linewidth=3)
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred)
plt.barh(x, dhdt, width,color='gold',label='1977-2007 dh/dt')
xl, yl = deg2_interp(dhdt,elev)
plt.plot(xl,yl,color='gold',linewidth=3)
dhdt,elev = dhdt_vs_elev_curve(dhdt0718_resampled)
plt.barh(x+30, dhdt, width,color='deeppink',label='2007-2018 dh/dt')
xl, yl = deg2_interp(dhdt,elev)
plt.plot(xl,yl,color='deeppink',linewidth=3)
plt.legend(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
#plt.plot(dhdt,elev)
plt.xlim(-3.5,0.5)
plt.ylim(700,3700)
plt.tight_layout()  
#plt.savefig('dhdt_elev_histogram_curvefit.png',bbox_inches='tight') 

plt.figure(figsize=(5,7))
x = np.arange(650,3751,100)
width = 30
#plt.figure(figsize=(6,8))
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2018_resampled)
plt.barh(x-30, dhdt, width,color='cornflowerblue',label='1977-2018 dh/dt')
xl, yl = deg3_interp(dhdt,elev)
plt.plot(xl,yl,color='cornflowerblue',linewidth=3)
dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred)
plt.barh(x, dhdt, width,color='gold',label='1977-2007 dh/dt')
xl, yl = deg3_interp(dhdt,elev)
plt.plot(xl,yl,color='gold',linewidth=3)
dhdt,elev = dhdt_vs_elev_curve(dhdt0718_resampled)
plt.barh(x+30, dhdt, width,color='deeppink',label='2007-2018 dh/dt')
xl, yl = deg3_interp(dhdt,elev)
plt.plot(xl,yl,color='deeppink',linewidth=3)
plt.legend(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
#plt.plot(dhdt,elev)
plt.xlim(-3.5,0.5)
plt.ylim(700,3700)
plt.tight_layout() 
#plt.savefig('dhdt_elev_histogram_curvefit_deg3.png',bbox_inches='tight')  

width = 0.3
plt.figure(figsize=(6,8))
dhdt1,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred,binnum=32)
x = np.arange(len(elev))
plt.barh(x, dhdt1, width,color='gold')
dhdt2,elev = dhdt_vs_elev_curve(dhdt0718_resampled,binnum=32)
plt.barh(x+0.3, dhdt2, width,color='deeppink')
#dhdt,elev = dhdt_vs_elev_curve(dhdt1977_2018_resampled,binnum=16)
plt.barh(x-0.3, ((np.array(dhdt1) + np.array(dhdt2))/2), width,color='darkblue')
plt.legend(['1977-2007 dh/dt','2007-2018 dh/dt','Mean dh/dt'],fontsize=14)
plt.yticks(x,elev, fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
#plt.plot(dhdt,elev)
plt.tight_layout()  
#plt.savefig('dhdt_elev_histogram_meandhdt_100mbins.png',bbox_inches='tight')   


plt.figure(figsize=(6,8))
x = np.arange(700,3701,200)
width = 60
dhdt1,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred,binnum=16)
plt.barh(x, dhdt1, width,color='gold',label='1977-2007 dh/dt')
dhdt2,elev = dhdt_vs_elev_curve(dhdt0718_resampled,binnum=16)
plt.barh(x+60, dhdt2, width,color='deeppink',label='2007-2018 dh/dt')
plt.barh(x-60, ((np.array(dhdt1) + np.array(dhdt2))/2), width,color='darkblue',label='Mean dh/dt')
x1, y1 = deg1_interp(((np.array(dhdt1) + np.array(dhdt2))/2),elev)
plt.plot(x1,y1,color='darkblue',linewidth=2)
x2, y2 = deg2_interp(((np.array(dhdt1) + np.array(dhdt2))/2),elev)
#plt.plot(x2,y2,color='darkblue',linewidth=2)
x3, y3 = deg3_interp(((np.array(dhdt1) + np.array(dhdt2))/2),elev)
#plt.plot(x3,y3,color='darkblue',linewidth=2)
plt.legend(fontsize=14)
#plt.legend(['1977-2007 dh/dt','2007-2018 dh/dt','Mean dh/dt'],fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
#plt.plot(dhdt,elev)
plt.xlim(-3.5,0.5)
plt.ylim(700,3700)
plt.tight_layout() 
#plt.savefig('dhdt_elev_histogram_meandhdt_deg1.png',bbox_inches='tight') 

# Revisions suggested by GF:
# 1. DONE. get rid of the bars in 1977--2007 200m bin that are influenced by the thinning at high elevs (feature of Map DEMs - not real)
# 2. DONE Calculate and plot time weighted mean bw 1977-2007 (30 yrs) and 2007-2018 (12 yrs)
# 3. Fit 2 or 3 degree polynomial to the 77-07, 07-18 and time weighted mean
# 4. Calculate and plot resulting dh/dt map from all 3 curves
# 5. Adjust curves as necessary to get correct dh/dt (from EMY?) - maybe gdalwarp kaskonly map to catchment map.

# find where exactly in the 1977-2007 map is the red blob (elevation band)
zl = np.where(Zgrid < 2400)
zu = np.where(Zgrid > 2700)

dhdt1977_2007_inferred[zl] = np.nan
dhdt1977_2007_inferred[zu] = np.nan

plt.figure()
plt.contourf(Xgrid,Ygrid,dhdt1977_2007_inferred,cmap = 'RdYlBu',levels=np.linspace(-4,4,81))
legend = plt.colorbar()

dhdt1,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred,binnum=16)
elev[np.where(~np.isnan(dhdt1))]
# shows that elev bands 2600 and 2800 are where dhdt is too negative 

# now remake figure with the incorrect bars replaced w those from the 2007-2018 dhdt map
Zgrid_dynamic_2phases, dhdt1977_2007_inferred = dynamic_surface_2phases()

plt.figure(figsize=(6,8))
x = np.arange(700,3701,200)
width = 60
dhdt1,elev = dhdt_vs_elev_curve(dhdt1977_2007_inferred,binnum=16)
dhdt2,elev = dhdt_vs_elev_curve(dhdt0718_resampled,binnum=16)
#replace dhdt1 vals at 2600 and 2800m elev with the vals from dhdt2
dhdt1[9:11] = dhdt2[9:11]
# calculate and plot time weighted mean over 42 yrs
meandhdt = ((30/42)*np.array(dhdt1)) + ((12/42)*np.array(dhdt2))
plt.barh(x, dhdt1, width,color='gold',label='1977-2007 dh/dt')
plt.barh(x+60, dhdt2, width,color='deeppink',label='2007-2018 dh/dt')
plt.barh(x-60,meandhdt, width,color='darkblue',label='Time-weighted mean dh/dt')
x1, y1 = deg2_interp(dhdt1,elev)
y1min=np.where(y1==np.min(y1))[0][0]
x2, y2 = deg2_interp(dhdt2,elev)
y2min=np.where(y2==np.min(y2))[0][0]
x3, y3 = deg2_interp(meandhdt,elev)
y3min=np.where(y3==np.min(y3))[0][0]
plt.plot(x1[y1min:],y1[y1min:]-(np.min(y1)-np.nanmin(DEM2018_resampled)),color='gold',linewidth=2,linestyle='--')
plt.plot(x2[y2min:],y2[y2min:]-(np.min(y2)-np.nanmin(DEM2018_resampled)),color='deeppink',linewidth=2,linestyle='--')
plt.plot(x3[y3min:],y3[y3min:]-(np.min(y3)-np.nanmin(DEM2018_resampled)),color='darkblue',linewidth=2,linestyle='--')
plt.legend(fontsize=14)
#plt.legend(['1977-2007 dh/dt','2007-2018 dh/dt','Time-weighted mean dh/dt'],fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Elevation Band (m a.s.l.)',fontsize=14)
plt.xlabel('dh/dt (m a$^{-1}$)',fontsize=14)
plt.xlim(-3.5,0.5)
#plt.plot(dhdt,elev)
plt.xlim(-3.5,0.5)
plt.ylim(600,3700)
plt.grid()
plt.tight_layout()
#plt.savefig('dhdt_elev_histogram_timeweightedmeandhdt.png',bbox_inches='tight')

# calculate distributed dh/dt from curves:
# 1. get the correct interpolated curve so that
def get_distributed_dhdt_from_curve(x,y,Zgrid):
    miny = np.where(y==np.min(y))[0][0]
    DPval = interp1d(y[miny:],x[miny:], kind = 'linear')    
    
    dhdt_smoothed = np.zeros(Zgrid.shape)
    for i in range(0,len(Zgrid)):
        for j in range(0,len(Zgrid[0])):
            if np.isnan(Zgrid[i,j]):
                dhdt_smoothed[i,j] = np.nan
            else:
                z = Zgrid[i,j]
                dhdt_smoothed[i,j] = DPval(z)
    dhdt_smoothed[offglacier] = 0
    return dhdt_smoothed

shiftedColorMap(matplotlib.cm.RdBu,start=0,midpoint=0.921,stop=0.8,name='customcmap')

Sfc_grid[np.where(np.isfinite(Sfc_grid))] = 0
Sfc_grid[np.where(np.isnan(Sfc_grid))] = 1

# Get Gapfilled 2018 DEM:
DEM2018_gapfilled = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/DEM2018_gapfilled.txt')

# Double check that this actually has a -0.46 m w.e. KASKAWULSH wide glacier mass balance. 
mean_dhdt2 = (get_distributed_dhdt_from_curve(x3-0.07,y3-(np.min(y3)-np.nanmin(DEM2018_gapfilled)),DEM2018_gapfilled))*0.9 # convert from m to m w.e. using ice density = 900 kg/m3
plt.figure(figsize=(10,6))
#plt.contourf(np.flipud(Sfc_grid), cmap = 'GnBu')
plt.contourf(Xgrid,Ygrid,mean_dhdt2,cmap = 'customcmap',levels=np.linspace(-1.75,0.15,39))
legend = plt.colorbar()
plt.axis('equal') 
legend.ax.set_ylabel('Elevation change (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.contour(Xgrid,Ygrid,Sfc_grid,levels=0,colors='k',linewidths=1)
plt.text(570000,6707000,'Off-glacier dh/dt = ' + str(np.round(np.nanmean(mean_dhdt2[offglacier]),4)) + ' m w.e.')
plt.text(570000,6703000,'Kaskawulsh dh/dt = ' + str(np.round(np.nanmean(mean_dhdt2[np.where(np.isfinite(catchment_tribarray))]),4)) + ' m w.e.')
plt.tight_layout()
#plt.savefig('dhdt_smoothed.png',bbox_inches='tight')   
#np.save('1977-2018_smoothed_dhdt.npy',mean_dhdt2)

