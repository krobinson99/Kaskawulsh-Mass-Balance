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
sys.path.insert(1,'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel')
from Model_functions_ver4 import regridXY_something
from scipy.interpolate import interp1d

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
    
def meanaccumulation_vs_z(Zlist,snowlist,z_start,z_end,delta_z=100):
    '''
    calculate the mean accumulation in each elevation band and return a list of mean accumulation vs mean
    elevation, to be plotted on top of the scatter plots
    '''
    Zarray = np.array(Zlist)
    snowarray = np.array(snowlist)
    mean_snow = []
    
    elevation_bins = []
    z = z_start
    while z <= z_end:
        elevation_bins.append(z)
        z += delta_z
        
    for zbin in elevation_bins:
        #print(zbin)
        z_bottom = zbin - (0.5*delta_z)
        z_top = zbin + (0.5*delta_z)
        #print(z_bottom, z_top)
        #zrange = Zarray[(Zarray >= z_bottom) & (Zarray < z_top)]
        zrange = (Zarray >= z_bottom) & (Zarray < z_top)
        #elevs = Zarray[zrange]
        #print(elevs)
        #print(np.min(elevs), np.max(elevs))
        snow_in_zrange = snowarray[zrange]
        mean_snow.append(np.nanmean(snow_in_zrange))
        
    return mean_snow, elevation_bins
    
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
plt.suptitle('All Snowpits',fontsize=14,y=1.01)
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


cluster_KW = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1000,4,4,4,18,18,18,22,22,22,22,22,22,22,22,22,22])
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
plt.suptitle('Catchment Snowpits',fontsize=14,y=1.01)
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

cluster_KW_NN = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1000,4,4,4,18,18,18,22,22,22,22,22,22,22,22,22,22])

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
plt.ylim(1200,3200)
plt.grid()
#plt.savefig('catchmentBC_vs_originalBC',bbox_inches = 'tight')

##################################################################################################
### BIAS CORRECT NARR ACCUMULATION USING THIS NEW CATCHMENT FORM AND COMPARE WITH OIB DATA!
##################################################################################################






