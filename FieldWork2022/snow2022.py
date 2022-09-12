# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 18:22:02 2022

script for looking at the may 2022 snow pit data collected by Gwenn, Tryggvi, Giovanni 

@author: katierobinson
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from netCDF4 import Dataset
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something

#load the .csv file with the snowpit data

snowpitdata ='F:/Mass Balance Model/Kaskawulsh-Mass-Balance/FieldWork2022/AccSurvey_KW_2022.csv'

df = pd.read_csv(snowpitdata)
    #df['Tributary'].astype('Int64')
arr = np.array(df)

z_obs = arr[:,4]
swe = arr[:,8] #units = cm
obs_snow = swe/100
easting_obs = arr[:,2]
northing_obs = arr[:,3]
pitnames = arr[:,1]



#compare field data to NARR accumulation
    # 0. get downscaled and bias corrected temp, precip, srad data for 2021-2022
            # download most recent 2022 data ###DONE
            # downscale most recent 2022 data ###in progress as of 2022-08-05: check 'F:/Mass Balance Model/DownscaledNARR_1979-2022' for completed 2022 files ###DONE
            # txt2netcdf the temp and prcp data for 2021-2022 ###Prcp in progress as of 2022-09-06: check 'F:/Mass Balance Model/DownscaledNARR_1979-2022' for completed 2021-2022 netSnow & temp files ###DONE
                #will probably need a new 2022 reffile since we only have 6 months of data and reffile is for 1 year. ###DONE: in Ref_files folder
            # bias correct the 2021-2022 data
    # 1. find nearest neighbour NARR gricell for each observation point (use the NN function?)
    # 2. compare the snow dpeth (in mwe) to both the bias corrected AND non bias-corrected NARR accumulation in that gridcell
    # 3. compare the recorded elevation to the Zgrid elevation of the closest gridcell
    # CHECKS: was there any snow in late summer 2021 that may have contributed at these sites?
             # does the end date of may 10 affect the result? (compare to may 4 vs may 14)

#load NARR accumulation for 2021-2022
Path2KWoutline = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel'
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

# 1. load net snow 2021 and 2022
def seasonalaccumulation(year1,year2,BC):
    if BC == False:
        Fvar = 'Precipitation'
    else:
        Fvar = 'Net snow'
        
    inF1 = Dataset(year1,'r')
    snowyr1 = inF1.variables[Fvar][:] # has shape (time,y,x)
    
    inF2 = Dataset(year2,'r')
    snowyr2 = inF2.variables[Fvar][:] # has shape (time,y,x)
    
    # calculate total accumulation from oct - jan of year 1
    #use np.nansum(axis=0) to get sum through time
    
    #formula = 8*(DOY-1)
    oct1_DOY= 274 
    oct1 = 8*(oct1_DOY-1)
    
    may10_DOY = 132 #what day did snow survey take place in 2022?
    may10 = 8*(may10_DOY-1)
    
    snowoct = np.nansum(snowyr1[oct1:,:,:],axis=0) #units are m.w.e.
    snowmay = np.nansum(snowyr2[:may10,:,:],axis=0)
    
    totalsnow = np.add(snowoct,snowmay)
    nanlocs = np.where(totalsnow == 0)
    totalsnow[nanlocs] = np.nan 
    
    #returns array of shape (Y,X) = (218,328) with total accumulation from Oct1 of year 1 to May10 of year 2
    return totalsnow

snow2022_bc =  seasonalaccumulation('F:/Mass Balance Model/BiasCorrectedNARR_1979-2022/Constant_Z/Prcpkaskonly_BC_2021.nc','F:/Mass Balance Model/BiasCorrectedNARR_1979-2022/Constant_Z/Prcpkaskonly_BC_2022.nc',True)
snow2022_nonbc =  seasonalaccumulation('F:/Mass Balance Model/DownscaledNARR_1979-2022/Constant_Z/netSnowkaskonly2021.nc','F:/Mass Balance Model/DownscaledNARR_1979-2022/Constant_Z/netSnowkaskonly2022.nc',False)

#snow20112012_noBC = seasonalaccumulation('F:/Mass Balance Model/DownscaledNARR/netSnowkaskonly2011.nc','F:/Mass Balance Model/DownscaledNARR/netSnowkaskonly2012.nc',False)
#snow20112012_BC = seasonalaccumulation('F:/Mass Balance Model/BiasCorrectedInputs/Kaskonly_R2S=1/Prcpkaskonly_BC_2011.nc','F:/Mass Balance Model/BiasCorrectedInputs/Kaskonly_R2S=1/Prcpkaskonly_BC_2011.nc',True)

plt.figure(figsize=(12,7))
plt.contourf(np.flipud(snow2022_bc), cmap = 'BuPu', levels = np.linspace(0,3,11))
#plt.contourf((distributed_accumulation_flipped), cmap = 'BuPu', levels = np.linspace(0,10,41))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.xticks(ticks=[0,50,100,150,200,250,300],labels=[int(Xgrid[0][0]),int(Xgrid[0][50]),int(Xgrid[0][100]),int(Xgrid[0][150]),int(Xgrid[0][200]),int(Xgrid[0][250]),int(Xgrid[0][300])],fontsize=10,rotation=20)
plt.yticks(ticks=[0,50,100,150,200],labels=[int(Ygrid[:,0][0]),int(Ygrid[:,0][50]),int(Ygrid[:,0][100]),int(Ygrid[:,0][150]),int(Ygrid[:,0][200])],fontsize=10,rotation=20)
plt.title('Bias corrected NARR Accumulation (Oct 2021 - May 2022)',fontsize=14)
#plt.savefig('NARRsnow_may2022_bc.png',bbox_inches = 'tight')

plt.figure(figsize=(12,7))
plt.contourf(np.flipud(snow2022_nonbc), cmap = 'BuPu', levels = np.linspace(0,3,11))
#plt.contourf((distributed_accumulation_flipped), cmap = 'BuPu', levels = np.linspace(0,10,41))
#plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.xticks(ticks=[0,50,100,150,200,250,300],labels=[int(Xgrid[0][0]),int(Xgrid[0][50]),int(Xgrid[0][100]),int(Xgrid[0][150]),int(Xgrid[0][200]),int(Xgrid[0][250]),int(Xgrid[0][300])],fontsize=10,rotation=20)
plt.yticks(ticks=[0,50,100,150,200],labels=[int(Ygrid[:,0][0]),int(Ygrid[:,0][50]),int(Ygrid[:,0][100]),int(Ygrid[:,0][150]),int(Ygrid[:,0][200])],fontsize=10,rotation=20)
plt.title('Uncorrected NARR Accumulation (Oct 2021 - May 2022)',fontsize=14)
#plt.savefig('NARRsnow_may2022_nonbc.png',bbox_inches = 'tight')

def calculate_XYZAcc(Xgrid,Ygrid,Zgrid,Acc_array):
    '''
    function returns four lists with the x,y,z coordinates of each grid cell 
    that has non-nan amount of snow plus a list with the total snow from 
    the acc array
    '''
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
                northing_EY.append(Ygrid[i,j])
                elevation_EY.append(Zgrid[i,j])
                acc_EY.append(Acc_array[i,j])
    
    return easting_EY,northing_EY,elevation_EY,acc_EY

easting_model, northing_model, z_model, narr_snow_bc = calculate_XYZAcc(Xgrid,np.flipud(Ygrid),Zgrid,snow2022_bc)

plt.figure(figsize=(12,7))
plt.scatter(easting_model,northing_model,c=narr_snow_bc,cmap="BuPu",s=3)
plt.clim(0,2.5)
plt.scatter(easting_obs,northing_obs,c=obs_snow, cmap="BuPu",marker='o',linewidth=0.5,edgecolor='black',s=100) #levels = np.linspace(0,4.5,19))
plt.clim(0,2.5)
legend = plt.colorbar() #colorbar goes to whichever is second
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('Bias Corrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)',fontsize=14)
#plt.savefig('NARR_bcvsSnowpits2022.png',bbox_inches = 'tight')

easting_model, northing_model, z_model, narr_snow_nonbc = calculate_XYZAcc(Xgrid,np.flipud(Ygrid),Zgrid,snow2022_nonbc)

plt.figure(figsize=(12,7))
plt.scatter(easting_model,northing_model,c=narr_snow_nonbc,cmap="BuPu",s=3)
plt.clim(0,2.5)
plt.scatter(easting_obs,northing_obs,c=obs_snow, cmap="BuPu",marker='o',linewidth=0.5,edgecolor='black',s=100) #levels = np.linspace(0,4.5,19))
plt.clim(0,2.5)
legend = plt.colorbar() #colorbar goes to whichever is second
legend.ax.set_ylabel('Accumulation (m w.e.)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('Uncorrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)',fontsize=14)
#plt.savefig('NARR_nonbcvsSnowpits2022.png',bbox_inches = 'tight')

#scatter plot w NARR on Y-Axis
#find nearest neighbour gridcells...
def nearestneighbour(NARRsnowarray,NASAsnowlist=depth,NASAeasting,NASAnorthing=northing,xgrid=Xgrid,ygrid=EY_Ygrid_flipped):    
    '''
    NARRsnowarray = (218,328) accumulation array
    '''
    distance_to_NN = []
    NARRaccumulation = []
    for i in range(0,len(NASAeasting)):
        if np.isnan(NASAsnowlist[i]):
            NARRaccumulation.append(np.nan)
        else:
            x_dist = xgrid - NASAeasting[i]
            y_dist = ygrid - NASAnorthing[i]
            distance = np.sqrt((x_dist**2)+(y_dist**2))
            closest_cell = np.where(distance == np.nanmin(distance))
            distance_to_NN.append(np.nanmin(distance))
            NARRaccumulation.append(NARRsnowarray[closest_cell])
        
    return NARRaccumulation

def nearestneighbour(NARRsnow,snowpitdata=obs_snow,snowpit_e=easting_obs,snowpit_n=northing_obs,snowpit_z=z_obs,xgrid=Xgrid,ygrid=np.flipud(Ygrid),Zgrid=Zgrid):
    distance_to_NN = []
    NARRaccumulation = []
    NARRelevs = []
    for i in range(0,len(snowpitdata)):
        x_dist = xgrid - snowpit_e[i]
        y_dist = ygrid - snowpit_n[i]
        distance = np.sqrt((x_dist**2)+(y_dist**2))
        closest_cell = np.where(distance == np.nanmin(distance))
        while np.isnan(NARRsnow[closest_cell]):
            print('snow is nan in closest cell. searching for next closest cell')
            distance[closest_cell] = 999
            closest_cell = np.where(distance == np.nanmin(distance))
        distance_to_NN.append(np.nanmin(distance))
        NARRaccumulation.append(NARRsnow[closest_cell])
        NARRelevs.append(Zgrid[closest_cell])
        
    return distance_to_NN, NARRaccumulation, NARRelevs

NN_narr_snow_bc = nearestneighbour(snow2022_bc)
NN_narr_snow_nonbc = nearestneighbour(snow2022_nonbc)

plt.figure(figsize=(12,6))
plt.suptitle('Nearest Neighbour Comparison',y=1)
plt.subplot(1,2,1)
plt.scatter(obs_snow,NN_narr_snow_nonbc[1],s=60)
plt.title('Uncorrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)')
plt.xlim(0,2)
plt.ylim(0,2)
plt.plot([0,2],[0,2],linestyle='--',color='k')
plt.xlabel('Snowpits (m w.e.)',fontsize=14)
plt.ylabel('NARR Accumulation (m.w.e.)',fontsize=14)
plt.subplot(1,2,2)
plt.scatter(obs_snow,NN_narr_snow_bc[1],s=60)
plt.title('Bias Corrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)')
plt.xlim(0,2)
plt.ylim(0,2)
plt.plot([0,2],[0,2],linestyle='--',color='k')
plt.xlabel('Snowpits (m w.e.)',fontsize=14)
plt.ylabel('NARR Accumulation (m.w.e.)',fontsize=14)
#plt.savefig('NearestNeighbour_snowpits.png',bbox_inches = 'tight')

plt.figure(figsize=(12,6))
plt.suptitle('Nearest Neighbour Comparison',y=1)
plt.subplot(1,2,1)
for i in range(0,len(z_obs)):
    plt.scatter(obs_snow[i],NN_narr_snow_nonbc[1][i],s=70)
plt.title('Uncorrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)')
plt.xlim(0,2)
plt.ylim(0,2)
plt.legend(pitnames)
plt.plot([0,2],[0,2],linestyle='--',color='k')
plt.xlabel('Snowpits (m w.e.)',fontsize=14)
plt.ylabel('NARR Accumulation (m.w.e.)',fontsize=14)
plt.subplot(1,2,2)
for i in range(0,len(z_obs)):
    plt.scatter(obs_snow[i],NN_narr_snow_bc[1][i],s=70)
plt.title('Bias Corrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)')
plt.xlim(0,2)
plt.ylim(0,2)
plt.plot([0,2],[0,2],linestyle='--',color='k')
plt.xlabel('Snowpits (m w.e.)',fontsize=14)
plt.ylabel('NARR Accumulation (m.w.e.)',fontsize=14)
#plt.savefig('NearestNeighbour_snowpits.png',bbox_inches = 'tight')

plt.figure(figsize=(6,6))
plt.title('Nearest Neighbour Gridcell Comparison',y=1)
#plt.subplot(1,2,1)
plt.scatter(z_obs,NN_narr_snow_nonbc[2],s=70)
#plt.title('Uncorrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)')
plt.xlim(1700,2600)
plt.ylim(1700,2600)
plt.plot([1700,2600],[1700,2600],linestyle='--',color='k')
plt.xlabel('Snowpit elevations (m a.s.l.)',fontsize=14)
plt.ylabel('NARR gridcell elevations (m a.s.l.)',fontsize=14)
#plt.savefig('snowpitelevations.png',bbox_inches = 'tight')

plt.figure(figsize=(6,6))
plt.title('Nearest Neighbour Gridcell Comparison',y=1)
#plt.subplot(1,2,1)
for i in range(0,len(z_obs)):
    plt.scatter(z_obs[i],NN_narr_snow_nonbc[2][i],s=70)
#plt.title('Uncorrected NARR Accumulation (Oct 2021-May 2022) \n vs Snowpit Data (May 2022)')
plt.xlim(1700,2800)
plt.ylim(1700,2800)
plt.legend(pitnames)
plt.plot([1700,2800],[1700,2800],linestyle='--',color='k')
plt.xlabel('Snowpit elevations (m a.s.l.)',fontsize=14)
plt.ylabel('NARR gridcell elevations (m a.s.l.)',fontsize=14)
#plt.savefig('snowpitelevations_coloured.png',bbox_inches = 'tight')