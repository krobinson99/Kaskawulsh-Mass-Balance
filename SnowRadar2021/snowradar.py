# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:08:00 2022

Download and plot data from NASA Icebride snow radar flyover of the Kaskawulsh May 3 2021
Data is available at https://data.cresis.ku.edu/data/snow/2021_Alaska_SO

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import utm
import os
import sys
sys.path.insert(1,'F:\Mass Balance Model\RawNARR_Catchment')
from Model_functions_ver4 import regridXY_something


# load .csv file with radar data
radardata = np.loadtxt('F:/Mass Balance Model/SnowRadar2021/Data_20210510_03.csv',delimiter=',',skiprows=1)

# load individual variables from the .csv file, eliminate last line since it is a Null value
lat = radardata[:,0][:-1]
lon = radardata[:,1][:-1]
elevation = radardata[:,4][:-1]
depth = radardata[:,6][:-1]/100
frame = radardata[:,5][:-1]

# convert lat/lon to easting/northing:
easting = []
northing =[]
for i in range(0,len(lat)):
    x = utm.from_latlon(lat[i],lon[i])
    easting.append(x[0])
    northing.append(x[1])

# get easting, northing, elevation, and accumulation data from Kask model domain
#load easting, northing, elevation:
Path2KWoutline = 'F:\Mass Balance Model\RawNARR_Catchment'
File_glacier_in = os.path.join(Path2KWoutline,'kaskonly_deb.txt')
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,2] 
debris_array = glacier[:,6]

#-------Turn vectors into 3D gridded inputs--------------------

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs = np.where(np.isnan(Zgrid))
#Xgrid[nanlocs] = np.nan
#Ygrid[nanlocs] = np.nan
Ygrid_flipped = np.flipud(Ygrid)

#Load final model accumulation field (averaged over all years, sims)
KW_accumulation = 'F:/Mass Balance Model/OUTPUTS/Plots/BaselineDebris/allrunAccumulations_BaselineDebris.npy'
uncorrected_accumulation_file = 'F:/Mass Balance Model/Kaskonly_Downscaled_NoBC/UncorrectedAccumulation.npy'

# plotting code from the PlotOutputs.py file
overall_Accumulation = np.load(KW_accumulation)
distributed_accumulation = np.nanmean(overall_Accumulation, axis = 0)
#distributed_accumulation_flipped = np.flipud(distributed_accumulation)

uncorrected_acc = np.load(uncorrected_accumulation_file)


#PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
plt.figure(figsize=(12,7))
plt.contourf(np.flipud(distributed_accumulation), cmap = 'BuPu', levels = np.linspace(0,10,41))
#plt.contourf((distributed_accumulation_flipped), cmap = 'BuPu', levels = np.linspace(0,10,41))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.xticks(ticks=[0,50,100,150,200,250,300],labels=[int(Xgrid[0][0]),int(Xgrid[0][50]),int(Xgrid[0][100]),int(Xgrid[0][150]),int(Xgrid[0][200]),int(Xgrid[0][250]),int(Xgrid[0][300])],fontsize=10,rotation=20)
plt.yticks(ticks=[0,50,100,150,200],labels=[int(Ygrid[:,0][0]),int(Ygrid[:,0][50]),int(Ygrid[:,0][100]),int(Ygrid[:,0][150]),int(Ygrid[:,0][200])],fontsize=10,rotation=20)
plt.title('Downscaled & Bias Corrected NARR Accumulation (2007-2018)',fontsize=14)
#plt.savefig('NARRAccumulation.png',bbox_inches = 'tight')

plt.figure(figsize=(12,7))
plt.scatter(easting,northing,c=depth, cmap="BuPu") #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('CReSIS Snow Radar (2021)',fontsize=14)
#plt.savefig('CReSISAccumulation.png',bbox_inches = 'tight')

NARR_elevation = []
NARR_accumulation = []
NARR_uncorrected_acc = []
NARR_easting = []
NARR_northing = []
for i in range(0,len(Xgrid)):
    for j in range(0,len(Xgrid[1])):
        if np.isnan(Zgrid[i,j]):
            pass
        else:
            NARR_elevation.append(Zgrid[i,j])
            NARR_accumulation.append(distributed_accumulation[i,j])
            NARR_uncorrected_acc.append(uncorrected_acc[i,j])
            NARR_easting.append(Xgrid[i,j])
            NARR_northing.append(Ygrid[i,j])
            #print(Ygrid[i,j])
            #what is wrong with this code!!
    
#----------------------------------------------------------------------------#
# REMOVE WEIRD RADAR POINTS
# METHOD 1: REMOVE THE FULL CROSS SECTION FROM NORTH TO WEST ACROSS THE DOMAIN
s=41000 #this is where wierdness starts
plt.scatter(easting[:s],northing[:s],c=depth[:s], cmap="BuPu") #levels = np.linspace(0,4.5,19))
    
e = 61000 #this is where wierdness ends
plt.scatter(easting[:e],northing[:e],c=depth[:e], cmap="BuPu") #levels = np.linspace(0,4.5,19))    
    
plt.scatter(easting[s:e],northing[s:e],c=depth[s:e], cmap="BuPu") #levels = np.linspace(0,4.5,19))        

badradar_flightpath = depth[s:e]
badradar_flightpath_E = easting[s:e]
badradar_flightpath_N = northing[s:e]    
plt.figure(figsize=(12,7))
plt.scatter(badradar_flightpath_E,badradar_flightpath_N,c=badradar_flightpath, cmap="BuPu") #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('CReSIS Snow Radar (2021) \n where it flies NW across the domain',fontsize=14)

# change points on the cross KW flight path to NaNs
depth_fixed = np.array(depth)
elevation_fixed = np.array(elevation)
easting_fixed = np.array(easting)
northing_fixed = np.array(northing)

depth_fixed[s:e] = np.nan
elevation_fixed[s:e] = np.nan
easting_fixed[s:e] = np.nan
northing_fixed[s:e] = np.nan

plt.figure(figsize=(12,7))
plt.scatter(easting_fixed,northing_fixed,c=depth_fixed, cmap="BuPu") #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('CReSIS Snow Radar (2021) \n without the NW cross-section',fontsize=14)
#plt.savefig('CReSISAccumulation.png',bbox_inches = 'tight')

plt.figure(figsize=(10,8))
plt.subplot(1,2,1)
plt.scatter(depth_fixed,elevation)

# METHOD 2: MANUALLY REMOVE FRAMES THAT LOOK WIERD IN THE ECHO: https://data.cresis.ku.edu/data/snow/2021_Alaska_SO/images/20210510_03/
bad_frames = np.array([11,12,13,14,15,25,26,27,28,29,30,31,32,33,34,35,36])
questionable_frames = [16,17,37] #check these ones with gwenn

frame0 = 2021051003000
bad_frames = bad_frames + frame0

for i in bad_frames:
    badpoints = np.where(frame == i)
    depth[badpoints] = np.nan
    elevation[badpoints] = np.nan
    
plt.figure(figsize=(12,7))
plt.scatter(easting,northing,c=depth, cmap="BuPu") #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('CReSIS Snow Radar (2021)',fontsize=14)
#plt.savefig('CReSISAccumulation.png',bbox_inches = 'tight')
    
plt.figure(figsize=(12,6))
plt.subplot(1,3,1)
plt.scatter(depth,elevation,marker='.')
plt.title('CReSIS Snow Radar (2021)',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,4.5)
plt.ylim(500,3500)
plt.subplot(1,3,2)
plt.scatter(NARR_accumulation,NARR_elevation,marker='.')
plt.title('NARR Accumulation (2007-2018) \n Downscaled + Bias Corrected',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,4.5)
plt.ylim(500,3500)
plt.subplot(1,3,3)
plt.scatter(NARR_uncorrected_acc,NARR_elevation,marker='.')
plt.title('NARR Accumulation (2007-2018) \n Uncorrected',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,4.5)
plt.ylim(500,3500)
#plt.savefig('NARR_CReSISvsElevation.png',bbox_inches = 'tight')

# PLOT THE CRESIS ACC. OVERLAYED ON THE NARR ACCUMULATION
#plt.figure(figsize=(14,7))
#plt.scatter(NARR_easting,np.flipud(NARR_northing),c=NARR_accumulation, cmap="BuPu")
#legend1 = plt.colorbar()
#legend1.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
#plt.scatter(easting,northing,c=depth, facecolor='k',marker='o',linewidth=3) #levels = np.linspace(0,4.5,19))
#legend2 = plt.colorbar()
#legend2.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
#plt.xlabel('Easting',fontsize=14)
#plt.ylabel('Northing',fontsize=14)
#plt.title('CReSIS Snow Radar Compared to NARR Accumulation',fontsize=14)
#plt.savefig('NARR_CReSIS_diffcolour.png',bbox_inches = 'tight')

#plt.figure(figsize=(8,5))
#plt.scatter(NARR_easting,np.flipud(NARR_northing),c=NARR_accumulation,cmap="BuPu")
#plt.scatter(easting,northing,c=depth, cmap="BuPu",marker='o',linewidth=3) #levels = np.linspace(0,4.5,19))
#legend = plt.colorbar()
#legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
#plt.xlabel('Easting',fontsize=14)
#plt.ylabel('Northing',fontsize=14)
#plt.title('CReSIS Snow Radar Compared to NARR Accumulation',fontsize=14)
#plt.savefig('NARR_CReSIS_samecolours.png',bbox_inches = 'tight')


plt.figure(figsize=(5,7))
plt.scatter(NARR_accumulation,NARR_elevation,marker='.',color='purple')
plt.scatter(depth,elevation,marker='.',color='turquoise')
plt.title('CReSIS Snow Radar (2021)',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.legend(['NARR Accumulation','CReSIS Snow Radar'])
plt.xlim(0,4.5)
plt.ylim(500,4000)
#plt.savefig('NARR_CReSIS_vs_elevation.png',bbox_inches = 'tight')

###########################################################################
Zgrid, EY_Xgrid, EY_Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs = np.where(np.isnan(Zgrid))

#Setup debris mask for use in radiation parameters
debris_grid, EY_Xgrid, EY_Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, debris_array)
debris_m = np.zeros(debris_grid.shape)
debris_m[np.where(debris_grid > 100)] = 1
debris_m[np.where(debris_grid <= 100)] = np.nan
debris_m[nanlocs] = np.nan

EY_Ygrid_flipped = np.flipud(EY_Ygrid) # EY_Ygrid is originally upside down (values decreasing south to north)

easting_EY = []
northing_EY = []
elevation_EY = []
acc_EY = []
for i in range(0,len(Xgrid)):
    for j in range(0,len(Xgrid[1])):
        if np.isnan(Zgrid[i,j]):
            pass
        else:
            easting_EY.append(Xgrid[i,j])
            northing_EY.append(EY_Ygrid_flipped[i,j])
            elevation_EY.append(Zgrid[i,j])
            acc_EY.append(distributed_accumulation[i,j])

plt.figure(figsize=(14,7))
plt.scatter(easting_EY,northing_EY,c=acc_EY,cmap="BuPu")
legend1 = plt.colorbar()
legend1.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.scatter(easting,northing,c=depth, facecolor='k',marker='o',linewidth=3) #levels = np.linspace(0,4.5,19))
legend2 = plt.colorbar()
#legend2.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('CReSIS Snow Radar Compared to NARR Accumulation',fontsize=14)
#plt.savefig('NARR_CReSIS_diffcolour.png',bbox_inches = 'tight')

plt.figure(figsize=(12,7))
plt.scatter(easting_EY,northing_EY,c=acc_EY,cmap="BuPu")
plt.scatter(easting,northing,c=depth, cmap="BuPu",marker='o',linewidth=3) #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('CReSIS Snow Radar Compared to NARR Accumulation',fontsize=14)
#plt.savefig('NARR_CReSIS_samecolours.png',bbox_inches = 'tight')
        

