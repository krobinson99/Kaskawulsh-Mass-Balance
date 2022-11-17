# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 19:48:33 2022

recreate the EMY accumulation bias correction

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
from netCDF4 import Dataset
import sys
import os
import scipy.io
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
sys.path.insert(1,'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel')
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import KWtributaries
sys.path.insert(2,'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021')
from OIB_picked_data import tributary_accvselev #importing from another file runs the whole file in this console
from OIB_picked_data import accumulationvselevation
from OIB_picked_data import meanaccumulation_vs_z
from OIB_picked_data import snow2021_bc_plussummersnow, snow2021_nobc_plussummersnow, kw_elevs, kw_snow_mwe
from OIB_picked_data import tribarray



    


# from Young et al. (2021): At each location, we calculate the difference between 
# measured (Cobs) and downscaled (Cds) seasonal accumulation on the date of measurement. 

# import snowpit data with dates and elevations filled out
# calculate yearly accumulation up to the date of observation
# find nearest neighbour downscaled cell
# append to list and calculate difference
# plot the difference versus elevation 
File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kask_catchment.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,4] 
Iy = glacier[:,5] 
Ih = glacier[:,6]  
sfc_type = glacier[:,8] 

#-------Turn vectors into 3D gridded inputs--------------------

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs = np.where(np.isnan(Zgrid))

#kask only gridding:
File_glacier_in_KW = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kaskonly.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier_KW = np.genfromtxt(File_glacier_in_KW, skip_header=1, delimiter=',')
        
Ix = glacier_KW[:,3] 
Iy = glacier_KW[:,4] 
Ih = glacier_KW[:,5]       
        
Zgrid_kw, Xgrid_kw, Ygrid_kw, xbounds_kw, ybounds_kw = regridXY_something(Ix, Iy, Ih)

#OG Bias Correction
DPval = interp1d(np.asarray([500, 900, 1100, 1300, 1500, 1700, 2100, 2300, 2500, \
        2700, 3100, 3700, 5000]), np.asarray([ 1.28530635,  1.26180458,  1.23360247,  1.25240269,  1.28293036,
        1.26320929,  1.24449735,  1.30082501,  2.46677552,  4.24804321,
        5.45333788,  5.45333788, 5.45333788]), kind = 'linear')

OG_corrections = np.asarray([ 1.28530635,  1.26180458,  1.23360247,  1.25240269,  1.28293036,
        1.26320929,  1.24449735,  1.30082501,  2.46677552,  4.24804321,
        5.45333788,  5.45333788, 5.45333788])
OG_elevations = np.asarray([500, 900, 1100, 1300, 1500, 1700, 2100, 2300, 2500, \
        2700, 3100, 3700, 5000])


df = pd.read_csv('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Accumulation/Snowpitdata_2007-2018.csv',encoding='cp1252')

snowpits = df['Snowpit']
sp_z = df['Elevation']
sp_easting = df['Easting']
sp_northing = df['Northing']
sp_mwe = df['Snow water eq. (m w.e.)']
sp_dates = df['Date']
sp_usealldata = df['Use all data']
sp_catchmentdata = df['Use catchment only']

sim_dates365 = pd.date_range(start="2007-01-01 00:00:00",end="2007-12-31 21:00:00",freq='3H')


def find_Cds(year1,year2,end_month,end_day,elevation):
    
    year1_file = 'F:/Mass Balance Model/Catchment/Catchment_Downscaled_NoBC/netSnowkaskonly' + str(year1) + '.nc'
    year2_file = 'F:/Mass Balance Model/Catchment/Catchment_Downscaled_NoBC/netSnowkaskonly' + str(year2) + '.nc'
    Fvar1 = 'Temperature'
    Fvar2 = 'Temperature'
    if year2 == 2022:
        year1_file = 'F:/Mass Balance Model/Downscaled_files/missing_trib/DownscaledNARR_1979-2022/Constant_Z/netSnowkaskonly' + str(year1) + '.nc'
        year2_file = 'F:/Mass Balance Model/Downscaled_files/missing_trib/DownscaledNARR_1979-2022/Constant_Z/netSnowkaskonly' + str(year2) + '.nc'
        Fvar1 = 'Precipitation'
        Fvar2 = 'Precipitation'
    
    inF1 = Dataset(year1_file,'r')
    snowyr1 = inF1.variables[Fvar1][:] # has shape (time,y,x)
    inF2 = Dataset(year2_file,'r')
    snowyr2 = inF2.variables[Fvar2][:] # has shape (time,y,x)
    
    snowstart = (np.where(sim_dates365 == np.datetime64(dt.datetime(2007,8,1))))[0][0]
    digging_index = (np.where(sim_dates365 == np.datetime64(dt.datetime(2007,end_month,end_day,12,00))))[0][0]
    
    snowaug = np.nansum(snowyr1[snowstart:,:,:],axis=0) #units are m.w.e.
    snowmay = np.nansum(snowyr2[:digging_index,:,:],axis=0)
    
    totalsnow = np.add(snowaug,snowmay)
    nanlocs = np.where(totalsnow == 0)
    totalsnow[nanlocs] = np.nan 
    if year2 == 2022:
        targetelevs = (Zgrid_kw >= (elevation - 1)) & (Zgrid_kw <= (elevation + 1)) 
    else:
        targetelevs = (Zgrid >= (elevation - 1)) & (Zgrid <= (elevation + 1))
        
    snow_at_targetz = np.nanmean(totalsnow[targetelevs])
    
    return snow_at_targetz, totalsnow    
    
#def meanaccumulation_vs_z(Zlist,snowlist,z_start,z_end,delta_z=100):
    '''
    calculate the mean accumulation in each elevation band and return a list of mean accumulation vs mean
    elevation, to be plotted on top of the scatter plots
    '''
#    Zarray = np.array(Zlist)
#    snowarray = np.array(snowlist)
#    mean_snow = []
    
#    elevation_bins = []
#    z = z_start
#    while z <= z_end:
#        elevation_bins.append(z)
#        z += delta_z
        
#    for zbin in elevation_bins:
        #print(zbin)
#        z_bottom = zbin - (0.5*delta_z)
#        z_top = zbin + (0.5*delta_z)
        #print(z_bottom, z_top)
        #zrange = Zarray[(Zarray >= z_bottom) & (Zarray < z_top)]
#        zrange = (Zarray >= z_bottom) & (Zarray < z_top)
        #elevs = Zarray[zrange]
        #print(elevs)
        #print(np.min(elevs), np.max(elevs))
#        snow_in_zrange = snowarray[zrange]
#        mean_snow.append(np.nanmean(snow_in_zrange))
        
#    return mean_snow, elevation_bins
    
C_obs = []
C_ds = []
C_z = []
pitnames = []
for i in range(0,len(snowpits)):
    #print(snowpits[i])
    if sp_usealldata[i] == 0:
        pass
    else:
        pitnames.append(snowpits[i])
        targetelev = sp_z[i]
        C_z.append(targetelev)
        C_obs.append(float(sp_mwe[i]))
        
        year = int(sp_dates[i][0:4])
        month = int(sp_dates[i][5:7])
        day = int(sp_dates[i][8:10]) 
        print(year,month,day)
        
        Cds = find_Cds(year-1,year,month,day,targetelev)[0]
        C_ds.append(Cds)
  
C_bias = np.array(C_obs) - np.array(C_ds)
C_multipfact = np.array(C_obs) / np.array(C_ds)
Cz = np.array(C_z)  

mean_Cbias = np.array(meanaccumulation_vs_z(Cz,C_bias,700,3000,200)[0])
mean_Cz = np.array(meanaccumulation_vs_z(Cz,C_bias,700,3000,200)[1])
nanmask = np.isfinite(mean_Cbias)

mean_Cmultifact = np.array(meanaccumulation_vs_z(Cz,C_multipfact,700,3000,200)[0])

cluster = np.array([1,1,99,1,1,2,2,2,99,1,1,2,2,1,1,2,2,1,1,2,2,1,1,3,2,2,0,0,4,4,4,2,13,13,13,0,0,0,100,100,100,0,0,0,0,0,0,0,0,0,0])
#99 = canada creek
# 0 = KW
# 100 = Lowell

plt.figure(figsize=(10,7))
plt.subplot(1,2,1)
plt.scatter(C_bias[cluster==1],Cz[cluster==1],c='magenta')
plt.scatter(C_bias[cluster==2],Cz[cluster==2],c='red')
plt.scatter(C_bias[cluster==3],Cz[cluster==3],c='navy')
plt.scatter(C_bias[cluster==4],Cz[cluster==4],c='cyan')
plt.scatter(C_bias[cluster==13],Cz[cluster==13],c='dodgerblue')
plt.scatter(C_bias[cluster==99],Cz[cluster==99],c='lawngreen')
plt.scatter(C_bias[cluster==0],Cz[cluster==0],c='blueviolet')
plt.scatter(C_bias[cluster==100],Cz[cluster==100],c='gold')
plt.legend(['GL1','GL2','GL3','GL4','GL13','Canada\nCreek','KW','LWL'],fontsize=14)
plt.plot(mean_Cbias[nanmask],mean_Cz[nanmask],c='k')
plt.xlabel('(C$_{obs}$ - C$_{ds}$) (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(-1.5,1.5)
plt.grid()
plt.subplot(1,2,2)
plt.scatter(C_multipfact[cluster==1],Cz[cluster==1],c='magenta')
plt.scatter(C_multipfact[cluster==2],Cz[cluster==2],c='red')
plt.scatter(C_multipfact[cluster==3],Cz[cluster==3],c='navy')
plt.scatter(C_multipfact[cluster==4],Cz[cluster==4],c='cyan')
plt.scatter(C_multipfact[cluster==13],Cz[cluster==13],c='dodgerblue')
plt.scatter(C_multipfact[cluster==99],Cz[cluster==99],c='lawngreen')
plt.scatter(C_multipfact[cluster==0],Cz[cluster==0],c='blueviolet')
plt.scatter(C_multipfact[cluster==100],Cz[cluster==100],c='gold')
plt.legend(['GL1','GL2','GL3','GL4','GL13','Canada\nCreek','KW','LWL'],fontsize=14)
plt.plot(mean_Cmultifact[nanmask],mean_Cz[nanmask],c='k')
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,6)
plt.grid()
plt.suptitle('All Snowpits \n Comparison by elevation (not location specific)',fontsize=14,y=1.03)
plt.tight_layout()
#plt.savefig('BC_allsnowpits.png',bbox_inches = 'tight')

# KW catchment pits only
  
C_obs_KW = []
C_ds_KW = []
C_z_KW = []
pitnames_KW = []
for i in range(0,len(snowpits)):
    #print(snowpits[i])
    if sp_catchmentdata[i] == 0:
        pass
    else:
        pitnames_KW.append(snowpits[i])
        targetelev = sp_z[i]
        C_z_KW.append(targetelev)
        C_obs_KW.append(float(sp_mwe[i]))
        
        year = int(sp_dates[i][0:4])
        month = int(sp_dates[i][5:7])
        day = int(sp_dates[i][8:10]) 
        print(year,month,day)
        
        Cds = find_Cds(year-1,year,month,day,targetelev)[0]
        C_ds_KW.append(Cds)
  
C_bias_KW = np.array(C_obs_KW) - np.array(C_ds_KW)
C_multipfact_KW = np.array(C_obs_KW) / np.array(C_ds_KW)
Cz_KW = np.array(C_z_KW)  

mean_Cbias_KW = np.array(meanaccumulation_vs_z(Cz_KW,C_bias_KW,1100,2700,200)[0])
mean_Cz_KW = np.array(meanaccumulation_vs_z(Cz_KW,C_bias_KW,1100,2700,200)[1])
nanmask = np.isfinite(mean_Cbias_KW)

mean_Cmultifact_KW = np.array(meanaccumulation_vs_z(Cz_KW,C_multipfact_KW,1100,2700,200)[0])


cluster_KW = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1000,4,4,4,18,18,18,22,22,22,22,22,22,22])
#99 = canada creek
# 0 = KW
# 100 = Lowell

plt.figure(figsize=(10,7))
plt.subplot(1,2,1)
plt.scatter(C_bias_KW[cluster_KW==1],Cz_KW[cluster_KW==1],c='deeppink')
plt.scatter(C_bias_KW[cluster_KW==4],Cz_KW[cluster_KW==4],c='dodgerblue')
plt.scatter(C_bias_KW[cluster_KW==18],Cz_KW[cluster_KW==18],c='gold')
plt.scatter(C_bias_KW[cluster_KW==22],Cz_KW[cluster_KW==22],c='blueviolet')
plt.scatter(C_bias_KW[cluster_KW==1000],Cz_KW[cluster_KW==1000],c='cyan')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14)
plt.plot(mean_Cbias_KW[nanmask],mean_Cz_KW[nanmask],c='k')
plt.xlabel('(C$_{obs}$ - C$_{ds}$) (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(-1.5,1.5)
plt.grid()
plt.subplot(1,2,2)
plt.scatter(C_multipfact_KW[cluster_KW==1],Cz_KW[cluster_KW==1],c='deeppink')
plt.scatter(C_multipfact_KW[cluster_KW==4],Cz_KW[cluster_KW==4],c='dodgerblue')
plt.scatter(C_multipfact_KW[cluster_KW==18],Cz_KW[cluster_KW==18],c='gold')
plt.scatter(C_multipfact_KW[cluster_KW==22],Cz_KW[cluster_KW==22],c='blueviolet')
plt.scatter(C_multipfact_KW[cluster_KW==1000],Cz_KW[cluster_KW==1000],c='cyan')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14)
plt.plot(mean_Cmultifact_KW[nanmask],mean_Cz_KW[nanmask],c='k')
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,6)
plt.grid()
plt.suptitle('Catchment Snowpits \n Comparison by elevation (not location specific)',fontsize=14,y=1.04)
plt.tight_layout()
#plt.savefig('BC_KWcatchment_snowpits.png',bbox_inches = 'tight')


###############################################################################
#COMPARE C_obs and C_ds at exact locations and times!!!!
###############################################################################

def nearestneighbour(NARRsnow,snowpitdata,snowpit_e,snowpit_n,xgrid,ygrid,Zgrid):
    '''
    NARRsnow should have shape (x,y)
    snowpitdata should be a single snow depth in m we
    snowpit_e, snowpit_n, snowpit_z is single values of snowpit coordinates
    Ygrid should be flipped with np.flipud()
    
    returns 3 lists: 
    1. distance to NN
    2. Narr accumulation for each NN cell
    3. elevation of each NN cell (all lists are same size as snowpit_e/n/z)
    '''
    distance_to_NN = []
    C_ds = []
    NARRelevs = []
    x_dist = xgrid - snowpit_e
    y_dist = ygrid - snowpit_n
    distance = np.sqrt((x_dist**2)+(y_dist**2))
    closest_cell = np.where(distance == np.nanmin(distance))
    while np.isnan(NARRsnow[closest_cell][0]):
        print('snow is nan in closest cell (' + str(np.nanmin(distance)) + '). searching for next closest cell')
        distance[closest_cell] = 999
        closest_cell = np.where(distance == np.nanmin(distance))
    distance_to_NN = np.nanmin(distance)
    C_ds = NARRsnow[closest_cell]
    NARRelevs = Zgrid[closest_cell]
        
    return distance_to_NN, C_ds, NARRelevs


NN_dists = []
C_ds_NN = []
C_ds_z = []  
for i in range(0,len(np.array(snowpits)[sp_catchmentdata==1])):
    print(np.array(snowpits)[sp_catchmentdata==1][i])
    
    #calculate the distirbuted accumulation field (ie sum of accumulation from
    # aug1 to date of snowpit digging)
    
    year = int(np.array(sp_dates)[sp_catchmentdata==1][i][0:4])
    month = int(np.array(sp_dates)[sp_catchmentdata==1][i][5:7])
    day = int(np.array(sp_dates)[sp_catchmentdata==1][i][8:10]) 
    print(year, month, day)
    
    targetelev = np.array(sp_z)[sp_catchmentdata==1][i]
    
    snow_to_date = find_Cds(year-1,year,month,day,targetelev)[1]
    
    C_obs = np.array(sp_mwe)[sp_catchmentdata==1][i]
    C_e = np.array(sp_easting)[sp_catchmentdata==1][i]
    C_n = np.array(sp_northing)[sp_catchmentdata==1][i]
    
    if year == 2022:
        NN_dist, NN_snow, NN_z = nearestneighbour(snow_to_date,C_obs,C_e,C_n,Xgrid_kw,np.flipud(Ygrid_kw),Zgrid_kw)
    else:
        NN_dist, NN_snow, NN_z = nearestneighbour(snow_to_date,C_obs,C_e,C_n,Xgrid,np.flipud(Ygrid),Zgrid)
    
    
    NN_dists.append(NN_dist)
    C_ds_NN.append(float(NN_snow))
    C_ds_z.append(float(NN_z))
    

NN_bias_KW = np.array(C_obs_KW) - np.array(C_ds_NN)
NN_multipfact_KW = np.array(C_obs_KW) / np.array(C_ds_NN)

mean_NNbias_KW= np.array(meanaccumulation_vs_z(Cz_KW,NN_bias_KW,1100,2700,200)[0])
mean_Cz_NN = np.array(meanaccumulation_vs_z(Cz_KW,NN_bias_KW,1100,2700,200)[1])
nanmask_NN = np.isfinite(mean_NNbias_KW)

mean_NNmultifact_KW = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1100,2700,200)[0])

plt.figure(figsize=(7,7))
for i in range(0,len(Cz_KW)):
    plt.scatter(Cz_KW[i],C_ds_z[i],s=70)
plt.legend(pitnames_KW,bbox_to_anchor=(1.1, 1.05))
plt.xlim(1700,2600)
plt.ylim(1700,2600)
plt.plot([1700,2600],[1700,2600],linestyle='--',color='k')
plt.xlabel('Snowpit elevations (m a.s.l.)',fontsize=14)
plt.ylabel('NARR gridcell elevations (m a.s.l.)',fontsize=14)
#plt.tight_layout()
plt.savefig('snowpitelevations.png',bbox_inches = 'tight')

cluster_KW_NN = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1000,4,4,4,18,18,18,22,22,22,22,22,22,22])

plt.figure(figsize=(10,7))
plt.subplot(1,2,1)
plt.scatter(NN_bias_KW[cluster_KW==1],Cz_KW[cluster_KW==1],c='deeppink')
plt.scatter(NN_bias_KW[cluster_KW==4],Cz_KW[cluster_KW==4],c='dodgerblue')
plt.scatter(NN_bias_KW[cluster_KW==18],Cz_KW[cluster_KW==18],c='gold')
plt.scatter(NN_bias_KW[cluster_KW==22],Cz_KW[cluster_KW==22],c='blueviolet')
plt.scatter(NN_bias_KW[cluster_KW==1000],Cz_KW[cluster_KW==1000],c='cyan')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14)
plt.plot(mean_NNbias_KW[nanmask],mean_Cz_NN[nanmask],c='k')
plt.xlabel('(C$_{obs}$ - C$_{ds}$) (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(-1.5,1.5)
plt.ylim(1200,2800)
plt.grid()
plt.subplot(1,2,2)
plt.scatter(NN_multipfact_KW[cluster_KW==1],Cz_KW[cluster_KW==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW==4],Cz_KW[cluster_KW==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW==18],Cz_KW[cluster_KW==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW==22],Cz_KW[cluster_KW==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW==1000],Cz_KW[cluster_KW==1000],c='cyan')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14)
plt.plot(mean_NNmultifact_KW[nanmask],mean_Cz_NN[nanmask],c='k')
#plt.plot(OG_corrections,OG_elevations)
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,6)
plt.ylim(1200,2800)
plt.grid()
plt.suptitle('Location-specific comparison (within catchment only)',fontsize=14,y=1.01)
plt.tight_layout()
#plt.savefig('Cobs_vs_Cds_locationspecific.png',bbox_inches = 'tight')

##################################################################################################
### HERE IS THE ELEVATION DEPENDANT BIAS CORRECTION BASED ON SNOWPITS WITHIN THE CATCHMENT ONLY:
##################################################################################################

catchmentBC_elevs = [500]
for i in mean_Cz_NN[nanmask]:
    catchmentBC_elevs.append(i)
catchmentBC_elevs.append(5000)

catchmentBC_corrfact = [mean_NNmultifact_KW[nanmask][0]]
for i in mean_NNmultifact_KW[nanmask]:
    catchmentBC_corrfact.append(i)
catchmentBC_corrfact.append(mean_NNmultifact_KW[nanmask][-1])

DPval_catchment = interp1d(np.asarray(catchmentBC_elevs), np.asarray(catchmentBC_corrfact), kind = 'linear')

plt.figure(figsize=(5,7))
#plt.plot(mean_NNmultifact_KW[nanmask],mean_Cz_KW[nanmask],c='k')
plt.plot(catchmentBC_corrfact,catchmentBC_elevs,c='k')
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
plt.legend(['Average (C$_{obs}$/C$_{ds}$)\nfrom all catchment snowpits','Original Bias Correction','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.5))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,6)
plt.ylim(1000,3000)
plt.grid()
#plt.savefig('catchmentBC_vs_originalBC',bbox_inches = 'tight')

##################################################################################################
##################################################################################################
##################################################################################################
### BIAS CORRECT NARR ACCUMULATION USING THIS NEW CATCHMENT FORM AND COMPARE WITH OIB DATA! ###
##################################################################################################
##################################################################################################
##################################################################################################

#bias correct the downscaled narr for 2020-2021 using the catchment bias correction:
def biascorrection(snow_array,Zgrid,DPvals):
    '''
    snow array should have shape (x,y) and be non bias corrected sum of accumulation for that season
    Zgrid should have same shape as snow array
    DPvals is the bias correction function
    
    returns the corrected snow array with same shape as the input snow_array
    '''
    nanlocs = np.where(np.isnan(snow_array))
    
    correctedsnow = np.zeros(snow_array.shape)
    correctedsnow[nanlocs] = np.nan
    
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            z = Zgrid[x][y] #find out the elevation at each x,y coordinate
            C_ds = snow_array[x][y]
            DeltaC = DPvals(z) # find the Bias Correction factor (DeltaC)
            correctedsnow[x][y] = C_ds*DeltaC
            
    return correctedsnow
    
snow2021_catchmentBC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DPval_catchment)

#get OIB data from each tributary (kw_snow_mwe and kw_elevs loaded from OIB_picked_data)
NAs_OIB = kw_snow_mwe[0:7845]
NAz_OIB = kw_elevs[0:7845]

SAs_OIB = kw_snow_mwe[7850:12255]
SAz_OIB = kw_elevs[7850:12255]

CAs_OIB = kw_snow_mwe[12260:-1]
CAz_OIB = kw_elevs[12260:-1]


#get average accumulations and elevation
avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],500,4000,delta_z=10)
avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],500,4000,delta_z=10)
avgsnow_catchBC,zbins_catchBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_catchmentBC)[1],accumulationvselevation(Zgrid_kw,snow2021_catchmentBC)[0],500,4000,delta_z=10)
avgsnow_OIB,zbins_OIB = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,1500,3000,delta_z=10)

plt.figure(figsize=(6,8))
plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],marker='.',color='orange')
plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta')
plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_catchmentBC)[0],accumulationvselevation(Zgrid_kw,snow2021_catchmentBC)[1],marker='.',color='deepskyblue')
plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='mediumaquamarine')
plt.legend(['NARR Accumulation: Uncorrected','NARR Accumulation w/ Original Bias Correction','NARR Accumulation w/ New Bias Correction','OIB Snow Radar (May 2021)'],fontsize=12,loc='lower right')
plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod')
plt.plot(avgsnow_BC,zbins_BC,color='indigo')
plt.plot(avgsnow_OIB,zbins_OIB,color='darkcyan')
plt.plot(avgsnow_catchBC,zbins_catchBC,color='blue')
plt.title('Kaskawulsh Accumulation May 2021',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,4.5)
plt.ylim(500,3700)
#plt.tight_layout()
#plt.savefig('KW_BiasCorrection_comparisons.png',bbox_inches = 'tight')

#calculate averages for each trib
# BC = ORIGINAL BIAS CORRECTION
# catchBC = NEW BIAS CORRECTION W CATCHMENT SNOWPITS ONLY
avgsnow_noBC_SA,zbins_noBC_SA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[0][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[0][0],500,4000,delta_z=10)
avgsnow_BC_SA,zbins_BC_SA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[0][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[0][0],500,4000,delta_z=10)
avgsnow_catchBC_SA,zbins_catchBC_SA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[0][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[0][0],500,4000,delta_z=10)
avgsnow_OIB_SA,zbins_OIB_SA = meanaccumulation_vs_z(SAz_OIB,SAs_OIB,500,4000,delta_z=10)

avgsnow_noBC_SW,zbins_noBC_SW = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[1][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[1][0],500,4000,delta_z=10)
avgsnow_BC_SW,zbins_BC_SW = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[1][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[1][0],500,4000,delta_z=10)
avgsnow_catchBC_SW,zbins_catchBC_SW = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[1][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[1][0],500,4000,delta_z=10)
#avgsnow_OIB_SA,zbins_OIB_SA = meanaccumulation_vs_z(SAz_OIB,SAs_OIB,500,4000,delta_z=10)

avgsnow_noBC_CA,zbins_noBC_CA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[2][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[2][0],500,4000,delta_z=10)
avgsnow_BC_CA,zbins_BC_CA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[2][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[2][0],500,4000,delta_z=10)
avgsnow_catchBC_CA,zbins_catchBC_CA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[2][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[2][0],500,4000,delta_z=10)
avgsnow_OIB_CA,zbins_OIB_CA = meanaccumulation_vs_z(CAz_OIB,CAs_OIB,500,4000,delta_z=10)

avgsnow_noBC_NA,zbins_noBC_NA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[3][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[3][0],500,4000,delta_z=10)
avgsnow_BC_NA,zbins_BC_NA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[3][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[3][0],500,4000,delta_z=10)
avgsnow_catchBC_NA,zbins_catchBC_NA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[3][1],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[3][0],500,4000,delta_z=10)
avgsnow_OIB_NA,zbins_OIB_NA = meanaccumulation_vs_z(NAz_OIB,NAs_OIB,500,4000,delta_z=10)

plt.figure(figsize=(10,12))
plt.subplot(2,2,1)
plt.title('South Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[0][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[0][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[0][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[0][1],color='orange')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[0][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[0][1],color='deepskyblue')
plt.scatter(SAs_OIB,SAz_OIB,color='mediumaquamarine')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR Accumulation: Uncorrected','NARR Accumulation w/ Original Bias Correction','NARR Accumulation w/ New Bias Correction','OIB Snow Radar (May 2021)'])
plt.plot(avgsnow_noBC_SA,zbins_noBC_SA,color='darkgoldenrod')
plt.plot(avgsnow_BC_SA,zbins_BC_SA,color='indigo')
plt.plot(avgsnow_catchBC_SA,zbins_catchBC_SA,color='blue')
plt.plot(avgsnow_OIB_SA,zbins_OIB_SA,color='darkcyan')
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.subplot(2,2,2)
plt.title('Stairway Glacier',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[1][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[1][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[1][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[1][1],color='orange')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[1][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[1][1],color='deepskyblue')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR Accumulation: Uncorrected','NARR Accumulation w/ Original Bias Correction','NARR Accumulation w/ New Bias Correction','OIB Snow Radar (May 2021)'])
plt.plot(avgsnow_noBC_SW,zbins_noBC_SW,color='darkgoldenrod')
plt.plot(avgsnow_BC_SW,zbins_BC_SW,color='indigo')
plt.plot(avgsnow_catchBC_SW,zbins_catchBC_SW,color='blue')
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.subplot(2,2,3)
plt.title('Central Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[2][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[2][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[2][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[2][1],color='orange')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[2][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[2][1],color='deepskyblue')
plt.scatter(CAs_OIB,CAz_OIB,color='mediumaquamarine')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR Accumulation: Uncorrected','NARR Accumulation w/ Original Bias Correction','NARR Accumulation w/ New Bias Correction','OIB Snow Radar (May 2021)'])
plt.plot(avgsnow_noBC_CA,zbins_noBC_CA,color='darkgoldenrod')
plt.plot(avgsnow_BC_CA,zbins_BC_CA,color='indigo')
plt.plot(avgsnow_catchBC_CA,zbins_catchBC_CA,color='blue')
plt.plot(avgsnow_OIB_CA,zbins_OIB_CA,color='darkcyan')
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.subplot(2,2,4)
plt.title('North Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[3][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_bc_plussummersnow)[3][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[3][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_nobc_plussummersnow)[3][1],color='orange')
plt.scatter(tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[3][0],tributary_accvselev(tribarray,Zgrid_kw,snow2021_catchmentBC)[3][1],color='deepskyblue')
plt.scatter(NAs_OIB,NAz_OIB,color='mediumaquamarine')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR Accumulation: Uncorrected','NARR Accumulation w/ Original Bias Correction','NARR Accumulation w/ New Bias Correction','OIB Snow Radar (May 2021)'])
plt.plot(avgsnow_noBC_NA,zbins_noBC_NA,color='darkgoldenrod')
plt.plot(avgsnow_BC_NA,zbins_BC_NA,color='indigo')
plt.plot(avgsnow_catchBC_NA,zbins_catchBC_NA,color='blue')
plt.plot(avgsnow_OIB_NA,zbins_OIB_NA,color='darkcyan')
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.tight_layout()
#plt.savefig('KW_BiasCorrection_Trib_comparisons.png',bbox_inches = 'tight')

################################################################################
# DEFINE A RANGE OF FUNCTIONAL FORMS TO DESCRIBE THE NEW ACCUMULATION BIAS CORRECTION
################################################################################
# define a range of functional forms describing the accumulation bias correction
   # by taking the mean and std of Cobs/Cds in each elevation band,
   # move through those stds to create a "min" and "max" (least aggressive vs most aggresive?)

#snow_bin100, z_bin100 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,100))
snow_bin100, z_bin100 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,100))
snow_bin150, z_bin150 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,150))
snow_bin200, z_bin200 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,200))
snow_bin250, z_bin250 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,250))
snow_bin300, z_bin300 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,300))
snow_bin350, z_bin350 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,350))
snow_bin400, z_bin400 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,400))
snow_bin450, z_bin450 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,450))
snow_bin500, z_bin500 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,500))
snow_bin550, z_bin550 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,550))
snow_bin600, z_bin600 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,560))
snow_bin2000, z_bin2000 = np.array(meanaccumulation_vs_z(Cz_KW,NN_multipfact_KW,1000,3000,2001))
#OIBsnow_bin250,OIBz_bin250 = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,1000,3000,250)

#make bias correction function out of each of the arrays above^
def generate_BiasCorrection(snow_vals,elev_vals):
    
    BC_elevs = [500]
    for i in elev_vals[np.isfinite(snow_vals)]:
        BC_elevs.append(i)
    BC_elevs.append(5000)

    BC_corrfact = [snow_vals[np.isfinite(snow_vals)][0]]
    for i in snow_vals[np.isfinite(snow_vals)]:
        BC_corrfact.append(i)
    BC_corrfact.append(snow_vals[np.isfinite(snow_vals)][-1])

    DPvals = interp1d(np.asarray(BC_elevs), np.asarray(BC_corrfact), kind = 'linear')
    
    return DPvals

DP_bin100 = generate_BiasCorrection(snow_bin100,z_bin100)
DP_bin150 = generate_BiasCorrection(snow_bin150,z_bin150)
DP_bin200 = generate_BiasCorrection(snow_bin200,z_bin200)
DP_bin250 = generate_BiasCorrection(snow_bin250,z_bin250)
DP_bin300 = generate_BiasCorrection(snow_bin300,z_bin300)
DP_bin350 = generate_BiasCorrection(snow_bin350,z_bin350)
DP_bin400 = generate_BiasCorrection(snow_bin400,z_bin400)
DP_bin450 = generate_BiasCorrection(snow_bin450,z_bin450)
DP_bin500 = generate_BiasCorrection(snow_bin500,z_bin500)
DP_bin550 = generate_BiasCorrection(snow_bin550,z_bin550)
DP_bin600 = generate_BiasCorrection(snow_bin600,z_bin600)
DP_bin2000 = generate_BiasCorrection(snow_bin2000,z_bin2000)

elevs = np.linspace(1000,3000,2000)

#plot piecewise linear functions with different elevation bins: (ie. 100,150,200,250), plus a linear
plt.figure(figsize=(5,7))
#plt.plot(mean_NNmultifact_KW[nanmask],mean_Cz_KW[nanmask],c='k')
#plt.plot(catchmentBC_corrfact,catchmentBC_elevs,c='k')
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
plt.plot(DP_bin150(elevs),elevs,linewidth=3,c='orange')
plt.plot(DP_bin200(elevs),elevs,linewidth=3,c='mediumblue')
plt.plot(DP_bin250(elevs),elevs,linewidth=3,c='deepskyblue')
plt.plot(DP_bin300(elevs),elevs,linewidth=3,c='red')
plt.plot(DP_bin2000(elevs),elevs,linewidth=3,c='springgreen')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
plt.legend(['Original Bias Correction','bin size = 150 m','bin size = 200 m','bin size = 250 m','bin size = 300 m','bin size = 2000 m','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
#plt.plot(NN_multipfact_KW,poly1d_fn(NN_multipfact_KW))
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid()
#plt.savefig('BC_with_diff_elevationbins.png',bbox_inches = 'tight')

####LINEAR FIT ON THE C_OBS/C_DS DATA##########################################
#linear fit to all the points
coef = np.polyfit(NN_multipfact_KW,Cz_KW,1)
poly1d_fn = np.poly1d(coef) 

#linear fit with forced intersection at the origin
def func(x, a):
    return a * x

xdata = NN_multipfact_KW
ydata = Cz_KW
popt, pcov = curve_fit(func, xdata, ydata)

plt.figure(figsize=(5,7))
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
plt.plot(np.linspace(0,5,len(NN_multipfact_KW)),poly1d_fn(np.linspace(0,5,len(NN_multipfact_KW))),linewidth=3,c='k')
plt.plot(np.linspace(0,5,len(NN_multipfact_KW)), func(np.linspace(0,5,len(NN_multipfact_KW)), popt),"b",linewidth=3)
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
plt.legend(['Original Bias Correction','Linear fit','Linear fit with y-intercept at origin','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid()
#plt.savefig('linearfitBC.png',bbox_inches = 'tight')

DP_linear = interp1d(np.array([0,func(4, popt)]),np.array([0,4]),kind = 'linear')


fig = plt.figure(figsize=(10, 10))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
sub1 = fig.add_subplot(2,3,1) # two rows, two columns, fist cell
#plt.annotate('sub1', xy = (0.5, 0.5), va = 'center', ha = 'center',  weight='bold', fontsize = 15)
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3, label='Original Bias Correction')
plt.plot(DP_bin150(elevs),elevs,linewidth=3,c='orange',label='bin size = 150 m')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
#plt.legend(['Original Bias Correction','bin size = 150 m','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.title('bin size = 150 m')
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid()
sub2 = fig.add_subplot(2,3,2) # two rows, two columns, second cell
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
plt.plot(DP_bin200(elevs),elevs,linewidth=3,c='mediumblue',label='bin size = 200 m')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
#plt.legend(['Original Bias Correction','bin size = 150 m','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.title('bin size = 200 m')
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid()
sub2 = fig.add_subplot(2,3,3) # two rows, two columns, second cell
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
plt.plot(DP_bin250(elevs),elevs,linewidth=3,c='deepskyblue',label='bin size = 250 m')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
#plt.legend(['Original Bias Correction','bin size = 150 m','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.title('bin size = 250 m')
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid() 
sub3 = fig.add_subplot(2,3,(4)) # two rows, two colums, combined third and fourth cell
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
plt.plot(DP_bin300(elevs),elevs,linewidth=3,c='red',label='bin size = 300 m')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
#plt.legend(['Original Bias Correction','bin size = 150 m','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.title('bin size = 300 m')
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid() 
sub3 = fig.add_subplot(2,3,(5)) # two rows, two colums, combined third and fourth cell
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
plt.plot(DP_bin2000(elevs),elevs,linewidth=3,c='springgreen',label='bin size = 2000 m')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink',label='GL1')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue',label='GL4')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold',label='KW 2018')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet',label='KW 2022')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan',label='Icefield')
#plt.legend(['Original Bias Correction','bin size = 150 m','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.title('bin size = 2000 m')
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid() 
fig.legend(fontsize=14,bbox_to_anchor=(0.99, 0.4))
#fig.legend(['Original Bias Correction','bin size = 150 m','bin size = 200 m','bin size = 250 m','bin size = 300 m','bin size = 2000 m','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.6))
plt.tight_layout()
#plt.savefig('C_diffbinsizes_subplots.png', dpi = 300, bbox_inches = 'tight')

binsizes = [100,150,200,250,300,350,400,450,500,550,999,2000]
i = 0
subplot = 1
fig = plt.figure(figsize=(12,11))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
for BC in [DP_bin100,DP_bin150,DP_bin200,DP_bin250,DP_bin300,DP_bin350,DP_bin400,DP_bin450,DP_bin500,DP_bin550,DP_linear,DP_bin2000]:
    fig.add_subplot(3,4,subplot)
    if subplot == 1:
        plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3, label='Original Bias Correction')
        plt.plot(BC(elevs),elevs,linewidth=3,c='k',label='New Bias Correction')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink',label='GL1')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue',label='GL4')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold',label='KW 2018')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet',label='KW 2022')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan',label='Icefield')
    else: 
        plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
        plt.plot(BC(elevs),elevs,linewidth=3,c='k')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
    plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=12)
    plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
    if binsizes[i] == 999:
        plt.title('linear fit bias correction',fontsize=11)
    else:
        plt.title('bin size = ' + str(binsizes[i]) + ' m',fontsize=11)
    plt.xlim(0,5)
    plt.ylim(1000,3000)
    plt.grid()
    i += 1
    subplot+=1
fig.legend(fontsize=10,bbox_to_anchor=(0.9, 0.025), ncol=8)
plt.tight_layout()
#plt.savefig('C_diffbinsizes_subplots_v2.png', dpi = 300, bbox_inches = 'tight')

#bias correct the 2021 snow data and plot along with the OIB data (binned with the same bin size)
snow2021_bin100BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin100)
snow2021_bin150BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin150)
snow2021_bin200BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin200)
snow2021_bin250BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin250)
snow2021_bin300BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin300)
snow2021_bin350BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin350)
snow2021_bin400BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin400)
snow2021_bin450BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin450)
snow2021_bin500BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin500)
snow2021_bin550BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin550)
snow2021_bin600BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin600)
snow2021_bin2000BC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_bin2000)
snow2021_linearBC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_linear)



binsizes = [100,150,200,250,300,350,400,450,500,550,600,2000]
i = 0
subplot = 1
fig = plt.figure(figsize=(12,11))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_bin600BC,snow2021_bin2000BC]:
    #print(str(BC), binsizes[i])
    #plt.figure(figsize=(6,8))
    fig.add_subplot(3,4,subplot)
    
    if subplot==1:
        avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],-500,5000,delta_z=binsizes[i])
        avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],-500,5000,delta_z=binsizes[i])
        avgsnow_binnedBC,avgz_binnedBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,BC)[1],accumulationvselevation(Zgrid_kw,BC)[0],-500,5000,delta_z=binsizes[i])
        avgsnow_OIB,zbins_OIB = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,500,5000,delta_z=binsizes[i])
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],marker='.',color='orange',label='Uncorrected NARR')
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta',label='Original Bias Correction')
        plt.scatter(accumulationvselevation(Zgrid_kw,BC)[0],accumulationvselevation(Zgrid_kw,BC)[1],marker='.',color='deepskyblue',label='New Bias Correction')
        plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='mediumaquamarine',label='OIB Snow Radar (May 2021)')
    else:
        avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],500,5000,delta_z=binsizes[i])
        avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],500,5000,delta_z=binsizes[i])
        avgsnow_binnedBC,avgz_binnedBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,BC)[1],accumulationvselevation(Zgrid_kw,BC)[0],500,5000,delta_z=binsizes[i])
        avgsnow_OIB,zbins_OIB = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,500,5000,delta_z=binsizes[i])
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],marker='.',color='orange')
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta')
        plt.scatter(accumulationvselevation(Zgrid_kw,BC)[0],accumulationvselevation(Zgrid_kw,BC)[1],marker='.',color='deepskyblue')
        plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='mediumaquamarine')
    
    #RMSE = np.sqrt(np.nanmean(np.square(np.subtract(np.array(avgsnow_binnedBC),np.array(avgsnow_OIB)))))
    #plt.text(1.75,1000,'RMSE = ' + str(np.round(RMSE,2)),fontsize=12)
    plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod',linewidth=3)
    plt.plot(avgsnow_BC,zbins_BC,color='indigo',linewidth=3)
    plt.plot(avgsnow_OIB,zbins_OIB,color='darkcyan',linewidth=3)
    plt.plot(avgsnow_binnedBC,avgz_binnedBC,color='blue',linewidth=3)
    plt.xlabel('Accumulation (m w.e.)',fontsize=12)
    plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
    plt.title('bin size used in\nbias correction = ' + str(binsizes[i]) + ' m',fontsize=11)
    plt.xlim(0,4.5)
    plt.ylim(500,3700)

    i+=1
    subplot+=1
fig.legend(fontsize=10,bbox_to_anchor=(0.9, 0.025), ncol=8)
plt.tight_layout()
#plt.savefig('binnedBCs_appliedtoNARR.png',bbox_inches = 'tight')

binsizes = [100,150,200,250,300,350,400,450,500,550,999,2000]
i = 0
subplot = 1
fig = plt.figure(figsize=(12,11))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_linearBC,snow2021_bin2000BC]:
    #print(str(BC), binsizes[i])
    #plt.figure(figsize=(6,8))
    fig.add_subplot(3,4,subplot)
    
    dz = 10
    avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],500,5000,delta_z=dz)
    avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],500,5000,delta_z=dz)
    avgsnow_binnedBC,avgz_binnedBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,BC)[1],accumulationvselevation(Zgrid_kw,BC)[0],500,5000,delta_z=dz)
    avgsnow_OIB,zbins_OIB = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,500,5000,delta_z=dz)
    if subplot==1:
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],marker='.',color='orange',label='Uncorrected NARR')
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta',label='Original Bias Correction')
        plt.scatter(accumulationvselevation(Zgrid_kw,BC)[0],accumulationvselevation(Zgrid_kw,BC)[1],marker='.',color='deepskyblue',label='New Bias Correction')
        plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='mediumaquamarine',label='OIB Snow Radar (May 2021)')
    else:
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],marker='.',color='orange')
        plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta')
        plt.scatter(accumulationvselevation(Zgrid_kw,BC)[0],accumulationvselevation(Zgrid_kw,BC)[1],marker='.',color='deepskyblue')
        plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='mediumaquamarine')
        
    #calculate RMSE:
    RMSE = np.sqrt(np.nanmean(np.square(np.subtract(np.array(avgsnow_binnedBC),np.array(avgsnow_OIB)))))
    
    plt.text(1.75,1000,'RMSE = ' + str(np.round(RMSE,2)),fontsize=12)
    plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod',linewidth=2)
    plt.plot(avgsnow_BC,zbins_BC,color='indigo',linewidth=2)
    plt.plot(avgsnow_OIB,zbins_OIB,color='teal',linewidth=2)
    plt.plot(avgsnow_binnedBC,avgz_binnedBC,color='blue',linewidth=2)
    plt.xlabel('Accumulation (m w.e.)',fontsize=12)
    plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
    if binsizes[i] == 999:
        plt.title('linear fit bias correction',fontsize=11)
    else:
        plt.title('bin size used in\nbias correction = ' + str(binsizes[i]) + ' m',fontsize=11)
    plt.xlim(0,4.5)
    plt.ylim(500,3700)
    
    i+=1
    subplot+=1

fig.legend(fontsize=10,bbox_to_anchor=(0.9, 0.025), ncol=8)
plt.tight_layout()
#plt.savefig('binnedBCs_appliedtoNARR_nosmoothing_dz10.png',bbox_inches = 'tight')

fig = plt.figure(figsize=(6,7))
dz = 10
BC = snow2021_linearBC
avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],500,5000,delta_z=dz)
avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],500,5000,delta_z=dz)
avgsnow_binnedBC,avgz_binnedBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,BC)[1],accumulationvselevation(Zgrid_kw,BC)[0],500,5000,delta_z=dz)
avgsnow_OIB,zbins_OIB = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,500,5000,delta_z=dz)
plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],marker='.',color='orange',label='Uncorrected NARR')
plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta',label='Original Bias Correction')
plt.scatter(accumulationvselevation(Zgrid_kw,BC)[0],accumulationvselevation(Zgrid_kw,BC)[1],marker='.',color='deepskyblue',label='New Bias Correction')
plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='mediumaquamarine',label='OIB Snow Radar (May 2021)')
#calculate RMSE:
RMSE = np.sqrt(np.nanmean(np.square(np.subtract(np.array(avgsnow_binnedBC),np.array(avgsnow_OIB)))))
plt.text(1.75,1000,'RMSE = ' + str(np.round(RMSE,2)),fontsize=12)
plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod',linewidth=2)
plt.plot(avgsnow_BC,zbins_BC,color='indigo',linewidth=2)
plt.plot(avgsnow_OIB,zbins_OIB,color='teal',linewidth=2)
plt.plot(avgsnow_binnedBC,avgz_binnedBC,color='blue',linewidth=2)
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.title('Linear bias correction',fontsize=12)
plt.xlim(0,4.5)
plt.ylim(500,3700)
plt.legend(fontsize=12,bbox_to_anchor=(0.5, 0.5), ncol=1)
plt.tight_layout()
#plt.savefig('NARRvsOIB_linearBC.png',bbox_inches = 'tight')
