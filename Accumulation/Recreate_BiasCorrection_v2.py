# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 19:48:33 2022

recreate the EMY accumulation bias correction

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import datetime as dt
from netCDF4 import Dataset
import sys
import os
import cmocean
import utm
import scipy.io
from PIL import Image
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.stats.stats import pearsonr  
from scipy.stats import linregress
sys.path.insert(1,'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel')
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import KWtributaries
from Model_functions_ver4 import model_domain
sys.path.insert(2,'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021')
from OIB_picked_data import tributary_accvselev #importing from another file runs the whole file in this console
from OIB_picked_data import accumulationvselevation
from OIB_picked_data import meanaccumulation_vs_z
from OIB_picked_data import snow2021_bc_plussummersnow, snow2021_nobc_plussummersnow, kw_elevs, kw_snow_mwe
from OIB_picked_data import tribarray
from OIB_picked_data import Colocated_avg_OIBsnow



plotting_accumulation_vs_elevation_RMSE = False #set false to avoid replotting all the outdated plots
  # (ie. the ones where RMSE is calculated based on elevation bins instead of co-located cells)


# from Young et al. (2021): At each location, we calculate the difference between 
# measured (Cobs) and downscaled (Cds) seasonal accumulation on the date of measurement. 

# import snowpit data with dates and elevations filled out
# calculate yearly accumulation up to the date of observation
# find nearest neighbour downscaled cell
# append to list and calculate difference
# plot the difference versus elevation 
Zgrid, Xgrid, Ygrid, xbounds, ybounds, Sfc = model_domain(catchment=True)
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


df = pd.read_csv('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Accumulation/Snowpitdata_2007-2022.csv',encoding='cp1252')

snowpits = df['Snowpit']
sp_z = df['Elevation']
sp_easting = df['Easting']
sp_northing = df['Northing']
sp_mwe = df['Snow water eq. (m w.e.)']
sp_dates = df['Date']
sp_usealldata = df['Use all data']
sp_catchmentdata = df['Use catchment only']
sp_regionalobs = df['Young et al. "Regional Obs"']
sp_youngetal = df['Young et al. Data']
sp_2022 = df['2022']

sim_dates365 = pd.date_range(start="2007-01-01 00:00:00",end="2007-12-31 21:00:00",freq='3H')


def find_Cds(year1,year2,end_month,end_day,elevation):
    
    year1_file = 'D:/Downscaled_files/Catchment/DynamicSurface/FinalDownscaling_1979-2022/netSnowkaskonly' + str(year1) + '.nc'
    year2_file = 'D:/Downscaled_files/Catchment/DynamicSurface/FinalDownscaling_1979-2022/netSnowkaskonly' + str(year2) + '.nc'
    Fvar1 = 'Precipitation'
    Fvar2 = 'Precipitation'
    
    inF1 = Dataset(year1_file,'r')
    snowyr1 = inF1.variables[Fvar1][:] # has shape (time,y,x)
    inF2 = Dataset(year2_file,'r')
    snowyr2 = inF2.variables[Fvar2][:] # has shape (time,y,x)
    
    snowstart = (np.where(sim_dates365 == np.datetime64(dt.datetime(2007,9,1))))[0][0]
    digging_index = (np.where(sim_dates365 == np.datetime64(dt.datetime(2007,end_month,end_day,12,00))))[0][0]
    
    snowaug = np.nansum(snowyr1[snowstart:,:,:],axis=0) #units are m.w.e.
    snowmay = np.nansum(snowyr2[:digging_index,:,:],axis=0)
    
    totalsnow = np.add(snowaug,snowmay)
    nanlocs = np.where(totalsnow == 0)
    totalsnow[nanlocs] = np.nan 

    targetelevs = (Zgrid >= (elevation - 1)) & (Zgrid <= (elevation + 1))
        
    snow_at_targetz = np.nanmean(totalsnow[targetelevs])
    
    return snow_at_targetz, totalsnow    
    
    
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

cluster = np.array([1,1,99,1,1,2,2,2,99,1,1,2,2,1,1,2,2,1,1,2,2,1,1,3,2,2,0,0,4,4,4,2,13,13,13,0,0,0,100,100,100,0,0,0,0,0,0,0,0,0,0,2018])
#99 = canada creek
# 0 = KW
# 100 = Lowell
# 2018 = Ochwat Core

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
plt.scatter(C_bias[cluster==2018],Cz[cluster==2018],c='springgreen')
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
plt.scatter(C_multipfact[cluster==2018],Cz[cluster==2018],c='springgreen')
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


cluster_KW = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1000,4,4,4,18,18,18,22,22,22,22,22,22,22,2018])
#99 = canada creek
# 0 = KW
# 100 = Lowell
#2018 = Ochwat core

plt.figure(figsize=(10,7))
plt.subplot(1,2,1)
plt.scatter(C_bias_KW[cluster_KW==1],Cz_KW[cluster_KW==1],c='deeppink')
plt.scatter(C_bias_KW[cluster_KW==4],Cz_KW[cluster_KW==4],c='dodgerblue')
plt.scatter(C_bias_KW[cluster_KW==18],Cz_KW[cluster_KW==18],c='gold')
plt.scatter(C_bias_KW[cluster_KW==22],Cz_KW[cluster_KW==22],c='blueviolet')
plt.scatter(C_bias_KW[cluster_KW==1000],Cz_KW[cluster_KW==1000],c='cyan')
plt.scatter(C_bias_KW[cluster_KW==2018],Cz_KW[cluster_KW==2018],c='springgreen')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield','KW Core (Ochwat et al. 2018)'],fontsize=14)
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
plt.scatter(C_multipfact_KW[cluster_KW==2018],Cz_KW[cluster_KW==2018],c='springgreen')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield','KW Core (Ochwat et al. 2018)'],fontsize=14)
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
#for i in range(0,len(Cz_KW)):
#    plt.scatter(Cz_KW[i],C_ds_z[i],s=100)
colors = plt.cm.viridis(np.linspace(0,1,len(Cz_KW)))
for i in range(0,len(Cz_KW)):
    plt.scatter(sorted(zip(Cz_KW,C_ds_z))[i][0],sorted(zip(Cz_KW,C_ds_z))[i][1],s=100,color=colors[i])
pitnames_KW_sorted = []
for i in sorted(zip(Cz_KW,pitnames_KW)):
        #print(i[1])
    pitnames_KW_sorted.append(i[1])
plt.legend(pitnames_KW_sorted,bbox_to_anchor=(1.1, 1.05),ncol=2,fontsize=14)
plt.xlim(1100,2800)
plt.ylim(1100,2800)
plt.plot([1100,2800],[1100,2800],linestyle='--',color='k')
plt.xlabel('Snowpit elevations (m a.s.l.)',fontsize=14)
plt.ylabel('NARR gridcell elevations (m a.s.l.)',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.tight_layout()
#plt.savefig('snowpitelevations.pdf',bbox_inches = 'tight')

cluster_KW_NN = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1000,4,4,4,18,18,18,22,22,22,22,22,22,22,2018])

plt.figure(figsize=(10,7))
plt.subplot(1,2,1)
plt.scatter(NN_bias_KW[cluster_KW==1],Cz_KW[cluster_KW==1],c='deeppink')
plt.scatter(NN_bias_KW[cluster_KW==4],Cz_KW[cluster_KW==4],c='dodgerblue')
plt.scatter(NN_bias_KW[cluster_KW==18],Cz_KW[cluster_KW==18],c='gold')
plt.scatter(NN_bias_KW[cluster_KW==22],Cz_KW[cluster_KW==22],c='blueviolet')
plt.scatter(NN_bias_KW[cluster_KW==1000],Cz_KW[cluster_KW==1000],c='cyan')
plt.scatter(NN_bias_KW[cluster_KW==2018],Cz_KW[cluster_KW==2018],c='springgreen')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield','KW Core (Ochwat et al. 2018)'],fontsize=14)
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
plt.scatter(NN_multipfact_KW[cluster_KW==2018],Cz_KW[cluster_KW==2018],c='springgreen')
plt.legend(['GL1','GL4','KW 2018','KW 2022','Icefield','KW Core (Ochwat et al. 2018)'],fontsize=14)
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
plt.scatter(NN_multipfact_KW[cluster_KW_NN==2018],Cz_KW[cluster_KW_NN==2018],c='springgreen')
plt.legend(['Average (C$_{obs}$/C$_{ds}$)\nfrom all catchment snowpits','Original Bias Correction','GL1','GL4','KW 2018','KW 2022','Icefield','KW Core (Ochwat et al. 2018)'],fontsize=14,bbox_to_anchor=(0.5, 0.5))
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
plt.scatter(NN_multipfact_KW[cluster_KW_NN==2018],Cz_KW[cluster_KW_NN==2018],c='springgreen')
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
def lin(x, a):
    return a * x

xdata = NN_multipfact_KW
ydata = Cz_KW
popt, pcov = curve_fit(lin, xdata, ydata)
DP_linear = interp1d(np.array([0,lin(4, popt)]),np.array([0,4]),kind = 'linear')

##############################################################################
# linear fit through point (1,0)


xdata = NN_multipfact_KW - 1
ydata = Cz_KW - 0
popt2, pcov2 = curve_fit(lin, xdata, ydata)
DP_linear_SL_v1 = interp1d(np.array([0,lin(4, popt2)]),np.array([1,5]),kind = 'linear')

lm = LinearRegression(fit_intercept = False)

# center data on x_0, y_0
x0 = 1
y0 = 0
x = NN_multipfact_KW
y = Cz_KW
x2 = NN_multipfact_KW - x0
y2 = Cz_KW - y0

# fit model
lm.fit(x2.reshape(-1, 1), y2)
# predict line
preds = lm.predict(np.arange(0, 4, 0.1).reshape(-1,1))
# plot on original scale
plt.figure(figsize=(5,7))
plt.plot(x, y, "o")
#plt.plot(x0, y0, "o")
# add x_0 and y_0 back to the predictions
#plt.plot(x2, y2, "o")
#plt.plot(np.arange(0, 3, 0.1), preds)
plt.plot(np.arange(0, 4, 0.1) + x0, preds  + y0)
plt.xlim(0,5)
plt.ylim(0,8000)
plt.grid()

DP_linear_SL_v2 = interp1d(np.array([0,lm.predict(np.array([3]).reshape(-1,1))]),np.array([1,4]),kind = 'linear')

###what if I just try adding like 1000 points of (1,0) to the arrays and then fit from there
x = NN_multipfact_KW
y = Cz_KW
i = 0
while i < 10000:
    x = np.append(x,1)
    y = np.append(y,0)
    i+=1
    
coef_SL = np.polyfit(x,y,1)
poly1d_fn_SL = np.poly1d(coef_SL) #same answer as curve_fit AND lm.predict

#hand code a linear fit to the data. then make sure it worked by comparing with the poly1d_fn fit
#linear fit through the origin (to get the slope) --> use data shifted so that the point
# we want the line to go through is centered (ie move (1,0) to (0,0):
def linearfit_origin(x,y):
    b0 = 0
    b1 = (np.mean(y)-b0)/np.mean(x)
    
    yfit = b0 + (b1*x)
    return b1, yfit

b1 = linearfit_origin(NN_multipfact_KW-1,Cz_KW)[0]

#use the slope (b1) calculated for the transformed data to calculate the intercept
# for the un-transformed data:
def linearfit_sealevel(b1,x):
    x0 = 1
    y0 = 0
    b0 = y0 - b1*x0
    yfit = b0 + (b1*x)
    return yfit

bestfit_linear = linearfit_sealevel(b1,NN_multipfact_KW)
#FINAL AND BEST LINEAR BC:
DP_linear_SL = interp1d(bestfit_linear,NN_multipfact_KW,kind = 'linear')

plt.figure(figsize=(5,7))
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
#plt.plot(np.linspace(0,5,len(NN_multipfact_KW)),poly1d_fn(np.linspace(0,5,len(NN_multipfact_KW))),linewidth=3,c='k')
#plt.plot(np.linspace(0,5,len(NN_multipfact_KW)), lin(np.linspace(0,5,len(NN_multipfact_KW)), popt),"b",linewidth=3)
plt.plot(DP_linear_SL_v2(elevs),elevs,linewidth=8,c='orange')
plt.plot(DP_linear_SL_v1(elevs),elevs,linewidth=4,c='blueviolet')
plt.plot(np.linspace(0,5,100),poly1d_fn_SL(np.linspace(0,5,100)),linewidth=2,c='turquoise',linestyle='--')
#plt.plot(NN_multipfact_KW,bestfit_linear,linewidth=3,c='crimson')
plt.plot(DP_linear_SL(elevs),elevs,linewidth=3,c='crimson')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==2018],Cz_KW[cluster_KW_NN==2018],c='springgreen')
plt.legend(['Original Bias Correction','Linear fit through (1,0) (sklearn.linear_model)','Linear fit through (1,0) (scipy.optimize)','Linear fit through (1,0) (numpy.polyfit)','Linear fit through (1,0) (hand code)','GL1','GL4','KW 2018','KW 2022','Icefield'],fontsize=14,bbox_to_anchor=(0.5, 0.7))
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,5)
plt.ylim(1000,3000)
plt.grid()
#plt.savefig('manylinearfit_BCs.png',bbox_inches = 'tight')

binsizes = [100,150,200,250,300,350,400,450,500,550,2000,999]
i = 0
subplot = 1
fig = plt.figure(figsize=(12,11))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
for BC in [DP_bin100,DP_bin150,DP_bin200,DP_bin250,DP_bin300,DP_bin350,DP_bin400,DP_bin450,DP_bin500,DP_bin550,DP_bin2000,DP_linear_SL]:
    fig.add_subplot(3,4,subplot)
    if subplot == 1:
        plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3, label='Original Bias Correction')
        plt.plot(BC(elevs),elevs,linewidth=3,c='k',label='New Bias Correction')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink',label='GL1')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue',label='GL4')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold',label='KW 2018')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet',label='KW 2022')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan',label='Icefield')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==2018],Cz_KW[cluster_KW_NN==2018],c='springgreen',label='KW Core (Ochwat et al. 2021)')
    else: 
        plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3)
        plt.plot(BC(elevs),elevs,linewidth=3,c='k')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan')
        plt.scatter(NN_multipfact_KW[cluster_KW_NN==2018],Cz_KW[cluster_KW_NN==2018],c='springgreen')
    plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=14)
    plt.xticks(np.arange(0,5,1))
    plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
    
    RMSE = np.sqrt(np.nanmean(np.square(np.subtract(NN_multipfact_KW,BC(Cz_KW)))))
    #plt.text(1.8,1250,'RMSE = ' + str(np.round(RMSE,2)),fontsize=12)
    
    if binsizes[i] == 999:
        plt.title('Linear regression',fontsize=14)
    else:
        plt.title('Bin size = ' + str(binsizes[i]) + ' m',fontsize=14)
    plt.xlim(0,4.5)
    plt.ylim(1000,3000)
    plt.grid(color='gainsboro',zorder=-10)
    i += 1
    subplot+=1
fig.legend(fontsize=14,bbox_to_anchor=(1.3,0.95), ncol=1, borderaxespad=0.19)
plt.tight_layout()
plt.savefig('C_diffbinsizes_subplots_v2.pdf', dpi = 300, bbox_inches = 'tight')

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
snow2021_linearBC = biascorrection(snow2021_nobc_plussummersnow,Zgrid_kw,DP_linear_SL)

if plotting_accumulation_vs_elevation_RMSE == True:
    binsizes = [100,150,200,250,300,350,400,450,500,550,10,2000]
    i = 0
    subplot = 1
    fig = plt.figure(figsize=(12,11))
    plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
    for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_linearBC,snow2021_bin2000BC]:
        #print(str(BC), binsizes[i])
        #plt.figure(figsize=(6,8))
        fig.add_subplot(3,4,subplot)
        
        if subplot==1:
            avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],500,5000,delta_z=binsizes[i])
            avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],500,5000,delta_z=binsizes[i])
            avgsnow_binnedBC,avgz_binnedBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,BC)[1],accumulationvselevation(Zgrid_kw,BC)[0],500,5000,delta_z=binsizes[i])
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
        
        RMSE = np.sqrt(np.nanmean(np.square(np.subtract(np.array(avgsnow_binnedBC),np.array(avgsnow_OIB)))))
        plt.text(1.75,1000,'RMSE = ' + str(np.round(RMSE,2)),fontsize=12)
        plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod',linewidth=3)
        plt.plot(avgsnow_BC,zbins_BC,color='indigo',linewidth=3)
        plt.plot(avgsnow_OIB,zbins_OIB,color='darkcyan',linewidth=3)
        plt.plot(avgsnow_binnedBC,avgz_binnedBC,color='blue',linewidth=3)
        plt.xlabel('Accumulation (m w.e.)',fontsize=12)
        plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
        if binsizes[i] == 10:
            plt.title('linear bias correction',fontsize = 11)
        else:
            plt.title('bin size used in\nbias correction = ' + str(binsizes[i]) + ' m',fontsize=11)
        plt.xlim(0,4.5)
        plt.ylim(500,3700)
    
        i+=1
        subplot+=1
    fig.legend(fontsize=10,bbox_to_anchor=(0.9, 0.025), ncol=8)
    plt.tight_layout()
    plt.savefig('binnedBCs_appliedtoNARR.png',bbox_inches = 'tight')
    
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
    plt.savefig('binnedBCs_appliedtoNARR_nosmoothing_dz10.png',bbox_inches = 'tight')
    
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
    plt.savefig('NARRvsOIB_linearBC.png',bbox_inches = 'tight')

else:
    pass

###############################################################################
# CALCULATE RMSE BASED ON CO-LOCATED GRIDCELLS, NOT BY BINNED ELEVATIONS
    # bin OIB cells into nearest NARR gridcell. Save into an array of shape Zgrid 
        #(in the corresponding cell)
###############################################################################

#for each point in the OIB data, find the closest model gridcell
# gridcell must be within np.sqrt(100**2 + 100**2) = 141.42 m for it to be
# within a NARR gridcell. if the closest cell is > 141.42m then its outside of the model domain.

# mean OIB snowdepth in each gridcell is calculated in the SnowRadar2021/OIB_picked_data.py script
# and imported at the top of this script (Colocated_avg_OIBsnow)

#function that returns lists of the OIB snow, and NARR snow in corresponding gridcells, plus an elevation list
def colocated_OIB_NARR_snow(NARR_accumulation,Zgrid=Zgrid_kw,OIB_accumulation=Colocated_avg_OIBsnow):

    matchingcells = np.where(np.isfinite(OIB_accumulation))
    
    OIBaccumulation = OIB_accumulation[matchingcells]
    NARRaccumulation = NARR_accumulation[matchingcells]
    elevation = Zgrid[matchingcells]

    
    return NARRaccumulation,OIBaccumulation,elevation


## USE ZGRID FOR GLACIER ONLY! this Zgrid is the catchment wide Zgrid. 
binsizes = [100,150,200,250,300,350,400,450,500,550,10,2000]
i = 0
subplot = 1
fig = plt.figure(figsize=(12,11))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_linearBC,snow2021_bin2000BC]:
    fig.add_subplot(3,4,subplot)
    NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)
    NARRsnow_noBC = colocated_OIB_NARR_snow(snow2021_nobc_plussummersnow)[0]
    NARRsnow_ogBC = colocated_OIB_NARR_snow(snow2021_bc_plussummersnow)[0]
    s = 10
    if subplot==1:
        plt.scatter(NARRsnow_noBC,z,color='orange',label='Uncorrected NARR',s=s)
        plt.scatter(NARRsnow_ogBC,z,color='darkmagenta',label='Original Bias Correction',s=s)
        plt.scatter(OIBsnow,z,color='mediumaquamarine',label='OIB Snow Radar (May 2021)',s=s)
        plt.scatter(NARRsnow_newBC,z,color='dodgerblue',label='New Bias Correction',s=s)
        #plt.scatter(OIBsnow,z,color='mediumaquamarine',label='OIB Snow Radar (May 2021)')
    else:
        plt.scatter(NARRsnow_noBC,z,color='orange',s=s)
        plt.scatter(NARRsnow_ogBC,z,color='darkmagenta',s=s)
        plt.scatter(OIBsnow,z,color='mediumaquamarine',s=s)
        plt.scatter(NARRsnow_newBC,z,color='dodgerblue',s=s)
        #plt.scatter(OIBsnow,z,color='mediumaquamarine')
    
    RMSE = np.sqrt(np.nanmean(np.square(np.subtract(np.array(NARRsnow_newBC),np.array(OIBsnow)))))
    plt.text(1.2,1700,'RMSE = ' + str(np.round(RMSE,2)),fontsize=10)
    #plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod',linewidth=3)
    #plt.plot(avgsnow_BC,zbins_BC,color='indigo',linewidth=3)
    #plt.plot(avgsnow_OIB,zbins_OIB,color='darkcyan',linewidth=3)
    #plt.plot(avgsnow_binnedBC,avgz_binnedBC,color='blue',linewidth=3)
    plt.xlabel('Accumulation (m w.e.)',fontsize=12)
    plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
    if binsizes[i] == 10:
        plt.title('linear bias correction',fontsize = 11)
    else:
        plt.title('bin size used in\nbias correction = ' + str(binsizes[i]) + ' m',fontsize=11)
    plt.xlim(0,2.5)
    plt.ylim(1600,2800)

    i+=1
    subplot+=1
fig.legend(fontsize=10,bbox_to_anchor=(0.9, 0.025), ncol=8)
plt.tight_layout()
#plt.savefig('colocated_accumulation_vs_elev.png',bbox_inches = 'tight')

###############################################################################
# scatter plot
###############################################################################

# get cluster of values representing what points belong to which trib:
colocatedcells = np.where(np.isfinite(Colocated_avg_OIBsnow))
tribcluster = tribarray[colocatedcells]
#    codes:
#    nontrib = 0
#    SA = 1
#    SW = 2
#    CA = 3
#    NA = 4
#    Trunk = 5
# the zeros in tribcluster should all actually be part of CA; change to 4
tribcluster[np.where(tribcluster==0)] = 4

binsizes = [100,150,200,250,300,350,400,450,500,550,2000,10]
i = 0
subplot = 1
fig = plt.figure(figsize=(16,12.2))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_bin2000BC,snow2021_linearBC]:
    fig.add_subplot(3,4,subplot)
    NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)
    #NARRsnow_noBC = colocated_OIB_NARR_snow(snow2021_nobc_plussummersnow)[0]
    #NARRsnow_ogBC = colocated_OIB_NARR_snow(snow2021_bc_plussummersnow)[0]
    s = 10
    if subplot==1:
        plt.scatter(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4],s=s,c='royalblue',label='North Arm')
        plt.scatter(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3],s=s,c='orange',label='Central Arm')
        plt.scatter(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1],s=s,c='crimson',label='South Arm')
    else:
        plt.scatter(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4],s=s,c='royalblue')
        plt.scatter(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3],s=s,c='orange')
        plt.scatter(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1],s=s,c='crimson')

    
    RMSEna = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4]))))
    RMSEca = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3]))))
    RMSEsa = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1]))))
    RMSEtot = np.sqrt(np.nanmean(np.square(np.subtract(np.array(NARRsnow_newBC),np.array(OIBsnow)))))
    
    #R2na = r2_score(OIBsnow[tribcluster==4],NARRsnow_newBC[tribcluster==4])
    #R2ca = r2_score(OIBsnow[tribcluster==3],NARRsnow_newBC[tribcluster==3])
    #R2sa = r2_score(OIBsnow[tribcluster==1],NARRsnow_newBC[tribcluster==1])
    #R2tot = r2_score(OIBsnow,NARRsnow_newBC)
    
    R2na = (pearsonr(OIBsnow[tribcluster==4],NARRsnow_newBC[tribcluster==4])[0])**2
    R2ca = (pearsonr(OIBsnow[tribcluster==3],NARRsnow_newBC[tribcluster==3])[0])**2
    R2sa = (pearsonr(OIBsnow[tribcluster==1],NARRsnow_newBC[tribcluster==1])[0])**2
    R2tot = (pearsonr(OIBsnow,NARRsnow_newBC)[0])**2
    
    #plt.text(1.15,0.95,'RMSE$_{NA}$ = ' + str(format(RMSEna,'.2f')),fontsize=14,color='royalblue')
    #plt.text(1.15,0.75,'RMSE$_{CA}$ = ' + str(format(RMSEca,'.2f')),fontsize=14,color='orange')
    #plt.text(1.15,0.55,'RMSE$_{SA}$ = ' + str(format(RMSEsa,'.2f')),fontsize=14,color='crimson')
    #plt.text(1.15,0.35,'RMSE$_{tot}$ = ' + str(format(RMSEtot,'.2f')),fontsize=14)
    
    plt.text(1.27,0.85,'RMSE (m w.e.):',fontsize=13,color='k')
    plt.text(1.45,0.7,'NA = ' + str(format(RMSEna,'.2f')),fontsize=14,color='royalblue')
    plt.text(1.45,0.6,'CA = ' + str(format(RMSEca,'.2f')),fontsize=14,color='orange')
    plt.text(1.45,0.5,'SA = ' + str(format(RMSEsa,'.2f')),fontsize=14,color='crimson')
    plt.text(1.45,0.4,'All = ' + str(format(RMSEtot,'.2f')),fontsize=14,color='k')
    
    #plt.text(1.4,0.95,'R$^{2}_{NA}$ = ' + str(format(R2na,'.3f')),fontsize=10,color='royalblue')
    #plt.text(1.4,0.75,'R$^{2}_{CA}$ = ' + str(format(R2ca,'.3f')),fontsize=10,color='orange')
    #plt.text(1.4,0.55,'R$^{2}_{SA}$ = ' + str(format(R2sa,'.3f')),fontsize=10,color='crimson')
    #plt.text(1.4,0.35,'R$^{2}_{tot}$ = ' + str(format(R2tot,'.3f')),fontsize=10)
    
    plt.xlabel('NARR accumulation (m w.e.)',fontsize=14)
    plt.ylabel('OIB snow depth (m w.e.)',fontsize=14)
    if binsizes[i] == 10:
        plt.title('Linear regression',fontsize = 14)
    else:
        plt.title('Bin size = ' + str(binsizes[i]) + ' m',fontsize=14)
    plt.xlim(0.3,2.1)
    plt.ylim(0.3,2.1)
    plt.plot([0,3],[0,3],color='grey',linestyle='--',zorder=-2)
    plt.xticks(np.arange(0.5,2.1,0.5),np.arange(0.5,2.1,0.5),fontsize=12)
    plt.yticks(np.arange(0.5,2.1,0.5),np.arange(0.5,2.1,0.5),fontsize=12)
    i+=1
    subplot+=1
fig.legend(fontsize=14,bbox_to_anchor=(0.7, 1.03), ncol=3,borderaxespad=0.1)
plt.tight_layout()
plt.savefig('colocated_NARR_vs_OIB.pdf',bbox_inches = 'tight')

NA_rmse = []
CA_rmse = []
SA_rmse = []
total_rmse = []
binsizes = ['100 m Bin','150 m Bin','200 m Bin','250 m Bin','300 m Bin','350 m Bin','400 m Bin','450 m Bin','500 m Bin','550 m Bin','2000 m Bin','Linear\nregression']
for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_bin2000BC,snow2021_linearBC]:
    NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)
    NA_rmse.append(np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4])))))
    CA_rmse.append(np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3])))))
    SA_rmse.append(np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1])))))
    total_rmse.append(np.sqrt(np.nanmean(np.square(np.subtract(np.array(NARRsnow_newBC),np.array(OIBsnow))))))

x = np.linspace(0,10,len(CA_rmse))
plt.figure(figsize=(12,5))
plt.bar(x-0.3,total_rmse,color='k',label='All',width=0.2)
plt.bar(x-0.1,CA_rmse,color='orange',label='CA',width=0.2)
plt.bar(x+0.1,NA_rmse,color='royalblue',label='NA',width=0.2)
plt.bar(x+0.3,SA_rmse,color='crimson',label='SA',width=0.2)
plt.xticks(ticks=np.linspace(0,10,len(CA_rmse)),labels=binsizes,rotation=75)
plt.xlabel('Bias Correction',fontsize=14)
plt.ylabel('RMSE (m w.e.)',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14,bbox_to_anchor=(1, 0.7))
#plt.savefig('RMSEbarchart.pdf',bbox_inches = 'tight')


NA_r2 = []
CA_r2 = []
SA_r2 = []
total_r2 = []
binsizes = ['100m bin','150m bin','200m bin','250m bin','300m bin','350m bin','400m bin','450m bin','500m bin','550m bin','2000m bin','linear']
for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_bin2000BC,snow2021_linearBC]:
#binsizes = ['100m bin','150m bin','450m bin']
#for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin450BC]:
    NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)

    NA_r2.append((pearsonr(OIBsnow[tribcluster==4],NARRsnow_newBC[tribcluster==4])[0])**2)
    CA_r2.append((pearsonr(OIBsnow[tribcluster==3],NARRsnow_newBC[tribcluster==3])[0])**2)
    SA_r2.append((pearsonr(OIBsnow[tribcluster==1],NARRsnow_newBC[tribcluster==1])[0])**2)
    total_r2.append((pearsonr(OIBsnow,NARRsnow_newBC)[0])**2)
    
x = np.linspace(0,12,len(binsizes))
#plt.figure(figsize=(12,5))
plt.figure(figsize=(12,5))
plt.bar(x-0.3,total_r2,color='k',label='R$^{2}_{total}$',width=0.2)
plt.bar(x-0.1,CA_r2,color='orange',label='R$^{2}_{CA}$',width=0.2)
plt.bar(x+0.1,NA_r2,color='royalblue',label='R$^{2}_{NA}$',width=0.2)
plt.bar(x+0.3,SA_r2,color='crimson',label='R$^{2}_{SA}$',width=0.2)
plt.xticks(ticks=x,labels=binsizes,rotation=75)
plt.xlabel('bias correction',fontsize=14)
plt.ylabel('R$^{2}$',fontsize=14)
plt.legend(fontsize=12,bbox_to_anchor=(1, 0.6))
#plt.savefig('R2barchart.png',bbox_inches = 'tight')

###############################################################################
def linfit(x,m,b):
    y = (m*x) + b
    return y

binsizes = [100,150,200,250,300,350,400,450,500,550,10,2000]
i = 0
subplot = 1
fig = plt.figure(figsize=(14,11))
plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
for BC in [snow2021_bin100BC,snow2021_bin150BC,snow2021_bin200BC,snow2021_bin250BC,snow2021_bin300BC,snow2021_bin350BC,snow2021_bin400BC,snow2021_bin450BC,snow2021_bin500BC,snow2021_bin550BC,snow2021_linearBC,snow2021_bin2000BC]:
    fig.add_subplot(3,4,subplot)
    NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)
    s = 10
    if subplot==1:
        plt.scatter(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4],s=s,c='royalblue',label='North Arm')
        plt.scatter(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3],s=s,c='orange',label='Central Arm')
        plt.scatter(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1],s=s,c='crimson',label='South Arm')
    else:
        plt.scatter(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4],s=s,c='royalblue')
        plt.scatter(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3],s=s,c='orange')
        plt.scatter(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1],s=s,c='crimson')

    R2na = (pearsonr(OIBsnow[tribcluster==4],NARRsnow_newBC[tribcluster==4])[0])**2
    R2ca = (pearsonr(OIBsnow[tribcluster==3],NARRsnow_newBC[tribcluster==3])[0])**2
    R2sa = (pearsonr(OIBsnow[tribcluster==1],NARRsnow_newBC[tribcluster==1])[0])**2
    R2tot = (pearsonr(OIBsnow,NARRsnow_newBC)[0])**2
    
    t = linregress(NARRsnow_newBC,OIBsnow)
    x = np.linspace(0,4,5)
    plt.plot(x,linfit(x,t.slope,t.intercept),linestyle='--',color='k')
    na = linregress(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4])
    plt.plot(x,linfit(x,na.slope,na.intercept),linestyle='--',color='royalblue')
    ca = linregress(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3])
    plt.plot(x,linfit(x,ca.slope,ca.intercept),linestyle='--',color='orange')
    sa = linregress(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1])
    plt.plot(x,linfit(x,sa.slope,sa.intercept),linestyle='--',color='crimson')

    plt.text(1.4,0.95,'M$_{NA}$ = ' + str(format(na.slope,'.3f')),fontsize=10,color='royalblue')
    plt.text(1.4,0.75,'M$_{CA}$ = ' + str(format(ca.slope,'.3f')),fontsize=10,color='orange')
    plt.text(1.4,0.55,'M$_{SA}$ = ' + str(format(sa.slope,'.3f')),fontsize=10,color='crimson')
    plt.text(1.4,0.35,'M$_{tot}$ = ' + str(format(t.slope,'.3f')),fontsize=10)
    plt.text(1.4,0.15,'1:1 line',fontsize=10,color='grey')
    
    plt.xlabel('NARR accumulation (m w.e.)',fontsize=10)
    plt.ylabel('OIB snow depth (m w.e.)',fontsize=10)
    if binsizes[i] == 10:
        plt.title('linear bias correction',fontsize = 11)
    else:
        plt.title('bin size used in\nbias correction = ' + str(binsizes[i]) + ' m',fontsize=11)
    plt.xlim(0,2.1)
    plt.ylim(0,2.1)
    plt.plot([0,3],[0,3],color='grey',linestyle='--',zorder=-2)
    
    i+=1
    subplot+=1
fig.legend(fontsize=14,bbox_to_anchor=(0.7, 1.03), ncol=3,borderaxespad=0.1)
plt.tight_layout()
#plt.savefig('colocated_NARR_vs_OIB_bestfits.png',bbox_inches = 'tight')

fig = plt.figure(figsize=(6.5,6.5))
BC = snow2021_bin450BC
NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)
NARRsnow_noBC, OIBsnow, z = colocated_OIB_NARR_snow(snow2021_nobc_plussummersnow)
NARRsnow_oldBC, OIBsnow, z = colocated_OIB_NARR_snow(snow2021_bc_plussummersnow)
plt.scatter(NARRsnow_noBC,OIBsnow,s=s,c='orange',label='Uncorrected NARR')
plt.scatter(NARRsnow_oldBC,OIBsnow,s=s,c='darkmagenta',label='Original Bias Correction')
plt.scatter(NARRsnow_newBC,OIBsnow,s=s,c='deepskyblue',label='New Bias Correction')
RMSE_new = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC,OIBsnow)))) 
RMSE_old = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_oldBC,OIBsnow)))) 
RMSE_none = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_noBC,OIBsnow)))) 
plt.text(1.3,0.72,'RMSE = ' + str(format(RMSE_none,'.2f')),fontsize=15,color='orange')
plt.text(1.3,0.60,'RMSE = ' + str(format(RMSE_old,'.2f')),fontsize=15,color='darkmagenta')
plt.text(1.3,0.48,'RMSE = ' + str(format(RMSE_new,'.2f')),fontsize=15,color='deepskyblue')
#plt.title('New Bias Correction',fontsize=14)
plt.xlim(0,2.1)
plt.ylim(0,2.1)
plt.xlabel('NARR accumulation (m w.e.)',fontsize=16)
plt.ylabel('OIB snow depth (m w.e.)',fontsize=16)
plt.plot([0,3],[0,3],color='grey',linestyle='--',zorder=-2)
plt.legend(fontsize=16, ncol=1,borderaxespad=0.1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
#plt.savefig('Newbiascorrection_colocatedscatter.png',bbox_inches = 'tight')

fig = plt.figure(figsize=(6.5,6.5))
BC = snow2021_bin450BC
NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)
NARRsnow_noBC, OIBsnow, z = colocated_OIB_NARR_snow(snow2021_nobc_plussummersnow)
NARRsnow_oldBC, OIBsnow, z = colocated_OIB_NARR_snow(snow2021_bc_plussummersnow)
plt.scatter(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4],s=s,c='royalblue',label='North Arm')
plt.scatter(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3],s=s,c='orange',label='Central Arm')
plt.scatter(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1],s=s,c='crimson',label='South Arm')
RMSEna = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==4],OIBsnow[tribcluster==4]))))
RMSEca = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==3],OIBsnow[tribcluster==3]))))
RMSEsa = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC[tribcluster==1],OIBsnow[tribcluster==1]))))
RMSEtot = np.sqrt(np.nanmean(np.square(np.subtract(np.array(NARRsnow_newBC),np.array(OIBsnow)))))
plt.text(1.3,0.72,'RMSE = ' + str(format(RMSEca,'.2f')),fontsize=15,color='orange')
plt.text(1.3,0.60,'RMSE = ' + str(format(RMSEna,'.2f')),fontsize=15,color='royalblue')
plt.text(1.3,0.48,'RMSE = ' + str(format(RMSEsa,'.2f')),fontsize=15,color='crimson')
#plt.title('New Bias Correction',fontsize=14)
plt.xlim(0,2.1)
plt.ylim(0,2.1)
plt.xlabel('NARR accumulation (m w.e.)',fontsize=16)
plt.ylabel('OIB snow depth (m w.e.)',fontsize=16)
plt.plot([0,3],[0,3],color='grey',linestyle='--',zorder=-2)
plt.legend(fontsize=16, ncol=1,borderaxespad=0.1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
#plt.savefig('Newbiascorrection_colocatedscatter_tribs.png',bbox_inches = 'tight')



fig = plt.figure(figsize=(5.5,7.5))
dz = 10
BC = snow2021_bin450BC
avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],500,5000,delta_z=dz)
avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],500,5000,delta_z=dz)
avgsnow_binnedBC,avgz_binnedBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid_kw,BC)[1],accumulationvselevation(Zgrid_kw,BC)[0],500,5000,delta_z=dz)
avgsnow_OIB,zbins_OIB = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,500,5000,delta_z=dz)
plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_nobc_plussummersnow)[1],marker='.',color='orange',label='Uncorrected NARR')
plt.scatter(accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid_kw,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta',label='Original Bias Correction')
plt.scatter(accumulationvselevation(Zgrid_kw,BC)[0],accumulationvselevation(Zgrid_kw,BC)[1],marker='.',color='deepskyblue',label='New Bias Correction')
plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='mediumaquamarine',label='OIB Snow Radar (May 2021)')
RMSE = np.sqrt(np.nanmean(np.square(np.subtract(np.array(avgsnow_binnedBC),np.array(avgsnow_OIB)))))    
plt.text(1.75,1500,'RMSE = ' + str(np.round(RMSE,4)),fontsize=14)
plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod',linewidth=2)
plt.plot(avgsnow_BC,zbins_BC,color='indigo',linewidth=2)
plt.plot(avgsnow_OIB,zbins_OIB,color='teal',linewidth=2)
plt.plot(avgsnow_binnedBC,avgz_binnedBC,color='blue',linewidth=2)
#plt.title('New Bias Correction',fontsize=1)
plt.xlabel('Accumulation (m w.e.)',fontsize=16)
plt.ylabel('Elevation (m a.s.l.)',fontsize=16)
plt.title('',fontsize=11)
plt.xlim(0,4.5)
plt.ylim(500,3700)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
fig.legend(fontsize=16,bbox_to_anchor=(1.2,0.3))
plt.tight_layout()
#plt.savefig('Newbiascorrection_accvselev.png',bbox_inches = 'tight')


fig = plt.figure(figsize=(5,6.5))
BC = DP_bin450
plt.plot(OG_corrections,OG_elevations,c='slategray',linestyle='--',linewidth=3, label='Original Bias Correction')
plt.plot(BC(elevs),elevs,linewidth=3,c='k',label='New Bias Correction')
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1],Cz_KW[cluster_KW_NN==1],c='deeppink',label='GL1',s=85)
plt.scatter(NN_multipfact_KW[cluster_KW_NN==4],Cz_KW[cluster_KW_NN==4],c='dodgerblue',label='GL4',s=85)
plt.scatter(NN_multipfact_KW[cluster_KW_NN==18],Cz_KW[cluster_KW_NN==18],c='gold',label='KW 2018',s=85)
plt.scatter(NN_multipfact_KW[cluster_KW_NN==22],Cz_KW[cluster_KW_NN==22],c='blueviolet',label='KW 2022',s=85)
plt.scatter(NN_multipfact_KW[cluster_KW_NN==1000],Cz_KW[cluster_KW_NN==1000],c='cyan',label='Icefield',s=85)
plt.scatter(NN_multipfact_KW[cluster_KW_NN==2018],Cz_KW[cluster_KW_NN==2018],c='springgreen',label='KW Core (Ochwat et al. 2021)',s=85)
plt.xlabel('(C$_{obs}$/C$_{ds}$)',fontsize=16)
plt.xticks(np.arange(0,5,1))
plt.ylabel('Elevation (m a.s.l.)',fontsize=16)
plt.xlim(0,4.5)
plt.ylim(1000,3000)
plt.grid(color='gainsboro',zorder=-10)
fig.legend(fontsize=14,bbox_to_anchor=(1.2,0.49))
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
#plt.savefig('Newbiascorrection_final.pdf',bbox_inches = 'tight')

plt.figure(figsize=(15,5))
plt.subplot(1,2,1)
plt.contourf(Xgrid_kw,Ygrid_kw,np.flipud(snow2021_nobc_plussummersnow), cmap = 'BuPu', levels = np.linspace(0,3.5,36))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
plt.axis('equal')
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=16)
plt.ylabel('Northing',fontsize=16)
plt.subplot(1,2,2)
plt.contourf(Xgrid_kw,Ygrid_kw,np.flipud(snow2021_bc_plussummersnow), cmap = 'BuPu', levels = np.linspace(0,3.5,36))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
plt.axis('equal')
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=16)
plt.ylabel('Northing',fontsize=16)
plt.tight_layout()
#plt.savefig('Corrected_vs_Uncorrected_snow_2021.png',bbox_inches = 'tight')

plt.figure(figsize=(10,5))
plt.contourf(Xgrid_kw,Ygrid_kw,np.flipud(snow2021_nobc_plussummersnow), cmap = 'BuPu', levels = np.linspace(0.3,0.7,41))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
plt.axis('equal')
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=16)
plt.ylabel('Northing',fontsize=16)
plt.tight_layout()
#plt.savefig('2021Snow_NoBC.pdf')


plt.figure(figsize=(12,7))
plt.contourf(np.flipud(snow2021_bc_plussummersnow), cmap = 'BuPu', levels = np.linspace(0,3.3,34))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=16)
plt.ylabel('Northing',fontsize=16)
plt.xticks(ticks=[0,50,100,150,200,250,300],labels=[int(Xgrid[0][0]),int(Xgrid[0][50]),int(Xgrid[0][100]),int(Xgrid[0][150]),int(Xgrid[0][200]),int(Xgrid[0][250]),int(Xgrid[0][300])],fontsize=10,rotation=20)
plt.yticks(ticks=[0,50,100,150,200],labels=[int(Ygrid[:,0][0]),int(Ygrid[:,0][50]),int(Ygrid[:,0][100]),int(Ygrid[:,0][150]),int(Ygrid[:,0][200])],fontsize=10,rotation=20)
#plt.title('Kaskawulsh Catchment Accumulation (2007-2008)')
#plt.savefig('2021Accumulation_OldBC.png')

plt.figure(figsize=(12,7))
plt.contourf(np.flipud(snow2021_bin450BC), cmap = 'BuPu', levels = np.linspace(0,3.3,34))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=16)
plt.ylabel('Northing',fontsize=16)
plt.xticks(ticks=[0,50,100,150,200,250,300],labels=[int(Xgrid[0][0]),int(Xgrid[0][50]),int(Xgrid[0][100]),int(Xgrid[0][150]),int(Xgrid[0][200]),int(Xgrid[0][250]),int(Xgrid[0][300])],fontsize=10,rotation=20)
plt.yticks(ticks=[0,50,100,150,200],labels=[int(Ygrid[:,0][0]),int(Ygrid[:,0][50]),int(Ygrid[:,0][100]),int(Ygrid[:,0][150]),int(Ygrid[:,0][200])],fontsize=10,rotation=20)
#plt.title('Kaskawulsh Catchment Accumulation (2007-2008)')
#plt.savefig('2021Accumulation_NewBC.png')

plt.figure(figsize=(12,7))
difference = snow2021_bin450BC-snow2021_bc_plussummersnow
plt.contourf(np.flipud(difference), cmap = 'coolwarm', levels = np.linspace(-2,2,21))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
legend = plt.colorbar()
legend.ax.set_ylabel('Difference in accumulation (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=16)
plt.ylabel('Northing',fontsize=16)
#nanlocs = np.where(np.isnan(snow2021_bin450BC))
#difference[nanlocs] = -999
#plt.contour(np.flipud(difference), colors = 'black', levels = -999)
plt.xticks(ticks=[0,50,100,150,200,250,300],labels=[int(Xgrid[0][0]),int(Xgrid[0][50]),int(Xgrid[0][100]),int(Xgrid[0][150]),int(Xgrid[0][200]),int(Xgrid[0][250]),int(Xgrid[0][300])],fontsize=10,rotation=20)
plt.yticks(ticks=[0,50,100,150,200],labels=[int(Ygrid[:,0][0]),int(Ygrid[:,0][50]),int(Ygrid[:,0][100]),int(Ygrid[:,0][150]),int(Ygrid[:,0][200])],fontsize=10,rotation=20)
#plt.title('Kaskawulsh Catchment Accumulation (2007-2008)')
plt.title('New Bias Correction - Old Bias Correction')
#plt.savefig('NewBC_OldBC_difference.png')


################################################################################
# FIG FOR EVS
################################################################################
MAE_none = np.nanmean(np.abs(np.subtract(NARRsnow_noBC,OIBsnow)))
MAE_old = np.nanmean(np.abs(np.subtract(NARRsnow_oldBC,OIBsnow)))
MAE_new = np.nanmean(np.abs(np.subtract(NARRsnow_newBC,OIBsnow)))

s=8
arial = {'fontname':'Arial'}
fig = plt.figure(figsize=(5,5))
BC = snow2021_bin450BC
NARRsnow_newBC, OIBsnow, z = colocated_OIB_NARR_snow(BC)
NARRsnow_noBC, OIBsnow, z = colocated_OIB_NARR_snow(snow2021_nobc_plussummersnow)
NARRsnow_oldBC, OIBsnow, z = colocated_OIB_NARR_snow(snow2021_bc_plussummersnow)
plt.scatter(NARRsnow_noBC,OIBsnow,s=s,c='orange',label='Uncorrected NARR')
plt.scatter(NARRsnow_oldBC,OIBsnow,s=s,c='darkmagenta',label='Regionally corrected NARR')
plt.scatter(NARRsnow_newBC,OIBsnow,s=s,c='deepskyblue',label='Locally corrected NARR')
RMSE_new = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_newBC,OIBsnow)))) 
RMSE_old = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_oldBC,OIBsnow)))) 
RMSE_none = np.sqrt(np.nanmean(np.square(np.subtract(NARRsnow_noBC,OIBsnow)))) 
plt.text(1.70,0.60,'MAE (m w.e.):',fontsize=12,color='k',**arial)
plt.text(2.06,0.50,str(format(MAE_none,'.2f')),fontsize=12,color='orange',**arial)
plt.text(2.06,0.40,str(format(MAE_old,'.2f')),fontsize=12,color='darkmagenta',**arial)
plt.text(2.06,0.30,str(format(MAE_new,'.2f')),fontsize=12,color='deepskyblue',**arial)
#plt.title('New Bias Correction',fontsize=14)
plt.xlim(0.25,2.25)
plt.ylim(0.25,2.25)
plt.xlabel('Modelled snow accumulation (m w.e.)',fontsize=12,**arial)
plt.ylabel('Airborne observed snow accumulation (m w.e.)',fontsize=12,**arial)
plt.plot([0,3],[0,3],color='grey',linestyle='--',zorder=-2,label='1:1 line')
plt.legend(fontsize=12, ncol=1,borderaxespad=0.1, prop={'family': 'Arial'},loc='upper right')
plt.xticks(fontsize=11,**arial)
plt.yticks(fontsize=11,**arial)
plt.tight_layout()
#plt.savefig('Snow4Flow_snowbias_300dpi.png',dpi=300,bbox_inches = 'tight')

# Per JH request: calculate distance from KW for the regional obs. 
# get non nan catchment gridcells
# calculate distance from snowpit to each gridcell
# get minimum distance for each snowpit
regional_distance = []
for i in range(0,len(snowpits[sp_regionalobs==1])):
    #print(np.array(snowpits)[sp_regionalobs==1][i])
    sp_x = np.array(sp_easting)[sp_regionalobs==1][i]
    sp_y = np.array(sp_northing)[sp_regionalobs==1][i]
    print(sp_x,sp_y)
    
    x_dist = Xgrid[np.where(~np.isnan(Zgrid))] - sp_x
    y_dist = Ygrid[np.where(~np.isnan(Zgrid))] - sp_y
    distance = np.sqrt((x_dist**2)+(y_dist**2))
    regional_distance.append(np.nanmin(distance))

print(np.min(np.array(regional_distance)/1000))
print(np.max(np.array(regional_distance)/1000))
print(np.mean(np.array(regional_distance)/1000))
print(np.std(np.array(regional_distance)/1000))
###############################################################################
# PLOT ALL THE LOCATIONS OF IN-SITU SNOW DEPTH DATA ON TOP OF ELEVATION MAP
# (PREFERABLY OF THE WHOLE AREA - NOT JUST KASK)
# maybe use arctic DEM if Etienne's DEMs don't cover the full area?

DEM1977_tiffile_v2 = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/1977/1977DEM.tif'
DEM1977_image_v2 = Image.open(DEM1977_tiffile_v2)
DEM1977_v2 = np.array(DEM1977_image_v2)
DEM1977_v2.shape

X = np.arange(549983.289493,646383.289493,40)
Y = np.arange(6772434.55099,6683594.55099,-40)

# add catchment outline
Sfc[np.where(Sfc ==1)] = 0
Sfc[np.where(np.isnan(Sfc))] = 1
# fix DEM so lowell is shown (Arctic DEM?)
plt.figure(figsize=(20,6.5))
plt.subplot(1,2,1)
plt.title('Snowpits Used in Young et al. (2021) Bias Correction',fontsize=14)
plt.contourf(X,Y,DEM1977_v2, cmap = 'Greys_r',levels=np.linspace(600,4800,22))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0)
plt.axis('equal')
plt.margins(x=0)
plt.margins(y=0)
plt.scatter(np.array(sp_easting)[np.where(sp_youngetal == 1)],np.array(sp_northing)[np.where(sp_youngetal == 1)],c=np.float64((np.array(sp_mwe)[np.where(sp_youngetal == 1)])),cmap='BuPu',s=110,edgecolor='black',linewidth=0.5)
legend2 = plt.colorbar()
legend2.ax.set_ylabel('Measured accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.ylabel('Northing',fontsize=14)
plt.xlabel('Easting',fontsize=14)
plt.subplot(1,2,2)
plt.title('Snowpits Used New Bias Correction',fontsize=14)
plt.contourf(X,Y,DEM1977_v2, cmap = 'Greys_r',levels=np.linspace(600,4800,22))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.axis('equal')
plt.margins(x=0)
plt.margins(y=0)
plt.scatter(np.array(sp_easting)[np.where(sp_catchmentdata == 1)],np.array(sp_northing)[np.where(sp_catchmentdata == 1)],c=np.float64((np.array(sp_mwe)[np.where(sp_catchmentdata == 1)])),cmap='BuPu',s=110,edgecolor='black',linewidth=0.5)
legend2 = plt.colorbar()
legend2.ax.set_ylabel('Measured accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0)
plt.ylabel('Northing',fontsize=14)
plt.xlabel('Easting',fontsize=14)
plt.tight_layout()
#plt.savefig('All_snowpit_locations.pdf',bbox_inches='tight')

# Plot ALL snowpits with border representing which BC was used
sp_usealldata[47]=0 # remove Lowell_L for plot to be cleaner
plt.figure(figsize=(10,6))
plt.contourf(X,Y,DEM1977_v2, cmap = 'Greys_r',levels=np.linspace(600,4800,22))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
legend.ax.tick_params(labelsize=14)
plt.axis('equal')
plt.margins(x=0)
plt.margins(y=0)
plt.scatter(np.array(sp_easting)[np.where(sp_usealldata == 1)],np.array(sp_northing)[np.where(sp_usealldata == 1)],c=np.float64((np.array(sp_mwe)[np.where(sp_usealldata == 1)])),cmap='BuPu',s=110,edgecolor='black',linewidth=0.5,zorder=2)
legend2 = plt.colorbar()
legend2.ax.set_ylabel('Measured accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
legend2.ax.tick_params(labelsize=14)
plt.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0,colors='white',zorder=1)
plt.ylabel('Northing',fontsize=14)
plt.xlabel('Easting',fontsize=14)
plt.tight_layout()
#plt.savefig('All_snowpit_locations.pdf',bbox_inches='tight')


plt.figure(figsize=(10,6.5))
#plt.title('Snowpits Used in Young et al. (2021) Bias Correction',fontsize=14)
plt.contourf(X,Y,DEM1977_v2, cmap = 'Greys_r',levels=np.linspace(600,4800,22))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
legend.ax.tick_params(labelsize=14)
plt.axis('equal')
plt.margins(x=0)
plt.margins(y=0)
plt.scatter(np.array(sp_easting)[np.where(sp_youngetal == 1)],np.array(sp_northing)[np.where(sp_youngetal == 1)],c='darkmagenta',s=110,edgecolor='grey',linewidth=0.5,label='Used in Young et al. (2021)\nbias correction only')
plt.scatter(np.array(sp_easting)[np.where(sp_catchmentdata == 1)],np.array(sp_northing)[np.where(sp_catchmentdata == 1)],c='deepskyblue',s=110,edgecolor='grey',linewidth=0.5,label='Used in new bias correction only')
plt.scatter(np.array(sp_easting)[np.where(sp_youngetal+sp_catchmentdata == 2)],np.array(sp_northing)[np.where(sp_youngetal+sp_catchmentdata == 2)],c='gold',s=110,edgecolor='grey',linewidth=0.5,label='Used in both corrections')
#plt.legend(['Data used in Young et al. (2021) bias correction','Data used in new bias correction','Data used in both corrections'])
plt.legend(fontsize=13)
#legend2 = plt.colorbar()
#legend2.ax.set_ylabel('Measured accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.ylabel('Northing',fontsize=14)
plt.xlabel('Easting',fontsize=14)
#plt.savefig('Snowpits_used_in_biascorrections.pdf',bbox_inches='tight')

plt.figure(figsize=(11,6.5))
plt.contourf(Xgrid_kw,np.flipud(Ygrid_kw),Zgrid_kw, cmap = 'Greys_r',levels=np.linspace(600,5000,43))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.contourf(Xgrid_kw,np.flipud(Ygrid_kw),Colocated_avg_OIBsnow)
plt.axis('equal')
plt.margins(x=0)
plt.margins(y=0)
#plt.scatter(np.array(sp_easting)[np.where(sp_usealldata == 1)],np.array(sp_northing)[np.where(sp_usealldata == 1)],c=np.float64((np.array(sp_mwe)[np.where(sp_usealldata == 1)])),cmap='BuPu',s=110,edgecolor='black',linewidth=0.5)
legend2 = plt.colorbar()
legend2.ax.set_ylabel('Measured accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.ylabel('Northing',fontsize=14)
plt.xlabel('Easting',fontsize=14)
plt.scatter(np.array(sp_easting)[np.where(sp_2022 == 1)],np.array(sp_northing)[np.where(sp_2022 == 1)],c='deepskyblue',s=110,edgecolor='grey',linewidth=0.5,label='Used in new bias correction only')
plt.tight_layout()

#Colocated_avg_OIBsnow
#plt.figure(figsize=(9,6))
#plt.title('2021 OIB Accumulation and 2022 Snowpits',fontsize=14)
#plt.contourf(Xgrid_kw,np.flipud(Ygrid_kw),distributed_accumulation, cmap = 'BuPu', levels = np.linspace(0,4,41))
#legend = plt.colorbar()
#legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
#plt.xlabel('Easting',fontsize=14)
#plt.ylabel('Northing',fontsize=14)
#plt.axis('equal')                     # This scales the x/y axes equally
#plt.title('Picks from data/misc/Alaska_seasonal_snow/2021_Alaska_seasonal_snow.mat',fontsize=14)
#plt.scatter(easting_OIB_test,northing_OIB_test,c=thickness, cmap = 'BuPu',s=0.1,marker='o')#linewidth=0.5,edgecolor='black')
#plt.clim(0,4) 
#plt.scatter(np.array(sp_easting)[np.where(sp_2022 == 1)],np.array(sp_northing)[np.where(sp_2022 == 1)],c=np.float64((np.array(sp_mwe)[np.where(sp_2022 == 1)])),cmap='BuPu',s=110,edgecolor='black',linewidth=0.5,vmin=0,vmax=4) 
#plt.tight_layout()
#plt.savefig('2022snowpits_OIBtrackoverlap.png')


################################################################################
# PLOT OIB 2021 FLIGHT PATH OVER KW (WHOLE PATH + DATA THAT WAS ACTUALLY PICKED)
################################################################################

# load .csv file with radar data
radardata = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/Data_20210510_03.csv',delimiter=',',skiprows=1)
# load individual variables from the .csv file, eliminate last line since it is a Null value
lat = radardata[:,0][:-1]
lon = radardata[:,1][:-1]
elevation = radardata[:,4][:-1]
depth = radardata[:,6][:-1]*0.3/100 #corrected with snow density 0.3 g/cm3 --> units = m w.e.
frame = radardata[:,5][:-1]

# convert lat/lon to easting/northing:
easting = []
northing = []
for i in range(0,len(lat)):
    x = utm.from_latlon(lat[i],lon[i])
    easting.append(x[0])
    northing.append(x[1])
    
Zgrid_kw, Xgrid_kw, Ygrid_kw, xbounds, ybounds, Sfc_kw = model_domain(catchment=False)
    

blackdot = mlines.Line2D([], [], color='k', marker='.', linestyle='None',
                          markersize=20, label='OIB Flight Path (May 10 2021)')

purpledot = mlines.Line2D([], [], color='lightblue', marker='.', linestyle='None',
                          markersize=20, label='Picked Data (Li et al., 2023)',alpha=1)

blackdash = mlines.Line2D([], [], color='k', linestyle='dashed',
                          markersize=20, label='Catchment Outline')

plt.figure(figsize=(13,6))
plt.contourf(Xgrid_kw,np.flipud(Ygrid_kw),Zgrid_kw, cmap = 'Greys_r',levels=np.linspace(600,4800,22))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0,colors='k',linestyles='dashed')
plt.axis('equal')
plt.ylabel('Northin (m)',fontsize=14)
plt.xlabel('Easting (m)',fontsize=14)
legend.ax.tick_params(labelsize=14)
plt.scatter(easting,northing,c='k',s=5,label='Flight path')
plt.scatter(Xgrid_kw,np.flipud(Ygrid_kw),c=Colocated_avg_OIBsnow,cmap='BuPu',s=20,vmin=0,vmax=1.8,label='Picked data')
legend2 = plt.colorbar()
legend2.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
legend2.ax.tick_params(labelsize=14)
plt.legend(handles=[blackdot,purpledot,blackdash],loc='lower left',fontsize=14)
#plt.savefig('OIBtracks.pdf',bbox_inches='tight')