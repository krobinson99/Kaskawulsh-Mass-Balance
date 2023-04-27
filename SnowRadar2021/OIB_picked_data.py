# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 13:01:38 2022

Processing the actual radar picks from the Kaskawulsh. Li et al (2022) https://doi.org/10.5194/egusphere-2022-368
discuss the processing 

@author: katierobinson
"""

import numpy as np
import h5py 
import scipy.io
import mat4py
import matplotlib.pyplot as plt
import utm
import os
import sys
import mat73
from pyproj import Proj
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import KWtributaries
from netCDF4 import Dataset
import simplekml

#############################################################################
# FOUND NEW FILE WITH PICKS FOR KW AND OTHER GLACIERS
# SEE https://data.cresis.ku.edu/data/misc/Alaska_seasonal_snow/ FOR MORE
##############################################################################

picks_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/2021_Alaska_seasonal_snow.mat'
seasonalsnow = scipy.io.loadmat(picks_file)
elevs = seasonalsnow['Elev'][0]
all_lats = seasonalsnow['Lat'][0]
all_lons = seasonalsnow['Lon'][0]
thickness = seasonalsnow['Thickness'][0]
time = seasonalsnow['GPS_time'][0]

plt.figure(figsize=(10,10))
plt.title('2021 Alaska Seasonal Snow Picks')
plt.xlabel('Lon')
plt.ylabel('Lat')
#plt.scatter(all_lons,all_lats,c=elevs,cmap="BuPu", vmin = 1000, vmax=3000)
plt.scatter(all_lons,all_lats,c=thickness,cmap="BuPu", vmin = 0, vmax=5,s=5)
plt.colorbar()

kwlon_l = np.where(all_lons<-139.8)
kwlat_l = np.where(all_lats<60.45)
kwlon_u = np.where(all_lons>-138.7)
kwlat_u = np.where(all_lats>60.85)

#all_lons[kwlon_l] = np.nan
#all_lons[kwlat_l] = np.nan
#all_lons[kwlon_u] = np.nan
#all_lons[kwlat_u] = np.nan


plt.figure(figsize=(10,10))
plt.title('2021 Alaska Seasonal Snow Picks')
plt.xlabel('Lon')
plt.ylabel('Lat')
#plt.scatter(all_lons,all_lats,c=elevs,cmap="BuPu", vmin = 1000, vmax=3000)
plt.scatter(all_lons,all_lats,c=thickness,cmap="BuPu", vmin = 0, vmax=5,s=1)
plt.colorbar()

###############################################################################
### WRITE COORDS TO A KML FILE ####

kml = simplekml.Kml()
for i in range(0,len(all_lons)):
    kml.newpoint(coords=[(all_lons[i],all_lats[i])])  # lon, lat, optional height
#kml.save("2021_Alaska_seasonal_snow.kml")

### 2021_Alaska_seasonal_snow.kml was loaded into ArcMap and coverted to a layer
# then the points just inside KW catchment were selected and coordinates were
# saved in the file OIB_KW_picks.csv 

#LOAD KML POINTS THAT WERE FOUND INSIDE THE CATCHMENT (OIB_KW_PICKS.CSV)
OIB_kw_picks = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/ARC/OIB_KW_picks.csv',delimiter=',',skiprows=1)
#load then match coords to closest coords in original list - make sure indices match
    
kw_lon_arc = OIB_kw_picks[:,0]
kw_lat_arc = OIB_kw_picks[:,1]


kw_lats_final = []
kw_lons_final = []
kw_elevs = []
kw_snowdepth = []
for i in range(0,len(kw_lat_arc)):
    distance_x = np.absolute(all_lons - kw_lon_arc[i])
    distance_y = np.absolute(all_lats - kw_lat_arc[i])
    closest_x = all_lons[np.where(distance_x == np.min(distance_x))[0][0]]
    closest_y = all_lats[np.where(distance_y == np.min(distance_y))[0][0]]
    if (np.where(distance_x == np.min(distance_x))[0][0] - np.where(distance_y == np.min(distance_y))[0][0]) == 0:
        kw_lats_final.append(closest_y)
        kw_lons_final.append(closest_x)
        kw_elevs.append(elevs[np.where(distance_x == np.min(distance_x))[0][0]])
        kw_snowdepth.append(thickness[np.where(distance_x == np.min(distance_x))[0][0]])
    else:
        pass
        #print('mismatch! ' + str(i))
    #print(np.where(distance_x == np.min(distance_x))[0][0],np.where(distance_y == np.min(distance_y))[0][0])
    #print(closest_x,closest_y)
    #problem: indice of closest x and closest y do not match, so hard to know where actual thickness is in the list

plt.figure()
plt.title('2021 KW Seasonal Snow Picks')
plt.xlabel('Lon')
plt.ylabel('Lat')
#plt.scatter(all_lons,all_lats,c=elevs,cmap="BuPu", vmin = 1000, vmax=3000)
plt.scatter(kw_lons_final,kw_lats_final,c=kw_snowdepth,cmap="BuPu", vmin = 0, vmax=5)
plt.colorbar()

#WAIT!!!!! FIRST, CONVERT SNOW DEPTH TO M.W.E.
#LI ET AL. (2021) USE A SNOW SURFACE DENSITY OF 377.36 kg.m3
ps = 377.36
pw = 1000
kw_snow_mwe = np.array(kw_snowdepth) * (ps/pw)

#plt.plot(kw_snowdepth)
#plt.plot(kw_snow_mwe)

#LOAD NARR ACCUMULATION:
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

#Load final model accumulation field (averaged over all years, sims)
KW_accumulation = 'F:/Mass Balance Model/OUTPUTS/Plots/BaselineDebris/allrunAccumulations_BaselineDebris.npy'
uncorrected_accumulation_file = 'F:/Mass Balance Model/Downscaled_files/missing_trib/Kaskonly_Downscaled_NoBC/UncorrectedAccumulation.npy'

# plotting code from the PlotOutputs.py file
overall_Accumulation = np.load(KW_accumulation)
distributed_accumulation = np.nanmean(overall_Accumulation, axis = 0)
#distributed_accumulation_flipped = np.flipud(distributed_accumulation)

uncorrected_acc = np.load(uncorrected_accumulation_file)


#PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
plt.figure(figsize=(12,7))
plt.contourf(np.flipud(distributed_accumulation), cmap = 'BuPu', levels = np.linspace(0,4,41))
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

easting_OIB_test = []
northing_OIB_test = []
for i in range(0,len(all_lons)):
    x = utm.from_latlon(all_lats[i],all_lons[i])
    easting_OIB_test.append(x[0])
    northing_OIB_test.append(x[1])
easting_OIB_test = np.array(easting_OIB_test)
northing_OIB_test = np.array(northing_OIB_test)

    
easting_OIB_test[kwlon_l] = np.nan
easting_OIB_test[kwlat_l] = np.nan
easting_OIB_test[kwlon_u] = np.nan
easting_OIB_test[kwlat_u] = np.nan

plt.figure(figsize=(10,6))
plt.contourf(Xgrid,Ygrid_flipped,distributed_accumulation, cmap = 'BuPu', levels = np.linspace(0,4,41))
legend = plt.colorbar()
legend.ax.set_ylabel('Thickness (m)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.axis('equal')                     # This scales the x/y axes equally
plt.title('Picks from data/misc/Alaska_seasonal_snow/2021_Alaska_seasonal_snow.mat',fontsize=14)
plt.scatter(easting_OIB_test,northing_OIB_test,c=thickness, cmap = 'BuPu',s=0.1,marker='o')#linewidth=0.5,edgecolor='black')
plt.clim(0,4)  
plt.tight_layout()
#plt.savefig('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/TestsforMichaelDaniel/2021_Alaska_seasonal_snow_on_KW.pdf')
################################################################################

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
    
    may10_DOY = 132
    may10 = 8*(may10_DOY-1)
    
    snowoct = np.nansum(snowyr1[oct1:,:,:],axis=0) #units are m.w.e.
    snowmay = np.nansum(snowyr2[:may10,:,:],axis=0)
    
    totalsnow = np.add(snowoct,snowmay)
    nanlocs = np.where(totalsnow == 0)
    totalsnow[nanlocs] = np.nan 
    
    return totalsnow

def calculate_XYZAcc(Xgrid,Ygrid,Zgrid,Acc_array):
    easting = []
    northing = []
    elevation = []
    acc = []
    for i in range(0,len(Xgrid)):
        for j in range(0,len(Xgrid[1])):
            if np.isnan(Zgrid[i,j]):
                pass
            else:
                easting.append(Xgrid[i,j])
                northing.append(Ygrid[i,j])
                elevation.append(Zgrid[i,j])
                acc.append(Acc_array[i,j])
    
    return easting,northing,elevation,acc

def accumulationvselevation(Zgrid,Accumulation_Array):
    elevation = []
    accumulation = []
    for i in range(0,len(Zgrid)):
        for j in range(0,len(Zgrid[1])):
            if np.isnan(Zgrid[i,j]):
                pass
            else:
                elevation.append(Zgrid[i,j])
                accumulation.append(Accumulation_Array[i,j])
    
    return accumulation,elevation

snow2021_nobc = seasonalaccumulation('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/downscaled2021acc/netSnowkaskonly2020.nc','F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/downscaled2021acc/netSnowkaskonly2021.nc',False)
snow2021_bc = seasonalaccumulation('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/downscaled2021acc/Prcpkaskonly_BC_2020.nc','F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/downscaled2021acc/Prcpkaskonly_BC_2021.nc',True)


wintersnow_bc, NARRz = accumulationvselevation(Zgrid,snow2021_bc)
wintersnow_nobc, NARRz = accumulationvselevation(Zgrid,snow2021_nobc)

def howmuchsnowinlatesummer(year1file,BC):
    if BC == False:
        Fvar = 'Precipitation'
    else:
        Fvar = 'Net snow'
        
    inF1 = Dataset(year1file,'r')
    snowyr1 = inF1.variables[Fvar][:] # has shape (time,y,x)
    
    #inF2 = Dataset(year2,'r')
    #snowyr2 = inF2.variables[Fvar][:] # has shape (time,y,x)
    
    # calculate total accumulation from oct - jan of year 1
    #use np.nansum(axis=0) to get sum through time
    
    #formula = 8*(DOY-1)
    aug1_DOY= 214
    aug1 = 8*(aug1_DOY-1)
    
    sept30_DOY = 274
    sept30 = 8*(sept30_DOY-1)
    
    summersnow = np.nansum(snowyr1[aug1:sept30,:,:],axis=0) #units are m.w.e.
    
    nanlocs = np.where(summersnow == 0)
    summersnow[nanlocs] = np.nan 
    
    return summersnow

summersnow_2020_bc = howmuchsnowinlatesummer('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/downscaled2021acc/Prcpkaskonly_BC_2020.nc',True)
snow2021_bc_plussummersnow = np.add(snow2021_bc,summersnow_2020_bc)

summersnow_2020_Nobc = howmuchsnowinlatesummer('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SnowRadar2021/downscaled2021acc/netSnowkaskonly2020.nc',False)
snow2021_nobc_plussummersnow = np.add(snow2021_nobc,summersnow_2020_Nobc)

plt.figure(figsize=(6,8))
plt.scatter(accumulationvselevation(Zgrid,snow2021_nobc)[0],accumulationvselevation(Zgrid,snow2021_nobc)[1],marker='.',color='orange')
plt.scatter(accumulationvselevation(Zgrid,snow2021_bc)[0],accumulationvselevation(Zgrid,snow2021_bc)[1],marker='.',color='darkmagenta')
plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='turquoise')
plt.title('Kaskawulsh Accumulation',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.legend(['NARR Accumulation (Oct2020-May2021) Uncorrected','NARR Accumulation (Oct2020-May2021) Bias Corrected','CReSIS Snow Radar (May 2021)'])
plt.xlim(0,4.5)
plt.ylim(500,4000)
#plt.savefig('NARRCresis2021_AccvsElev.png',bbox_inches = 'tight')


plt.figure(figsize=(6,8))
plt.scatter(accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[1],marker='.',color='orange')
plt.scatter(accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta')
plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='turquoise')
plt.title('Kaskawulsh Accumulation',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.legend(['NARR Accumulation (Aug2020-May2021) Uncorrected','NARR Accumulation (Aug2020-May2021) Bias Corrected','CReSIS Snow Radar (May 2021)'])
plt.xlim(0,4.5)
plt.ylim(500,4000)
#plt.savefig('NARRCresis2021_AccvsElev_summersnowincl.png',bbox_inches = 'tight')

easting_EY,northing_EY,elevation_EY,snow_bc_aug2may_list = calculate_XYZAcc(Xgrid,np.flipud(Ygrid),Zgrid,snow2021_bc_plussummersnow)

easting_OIB = []
northing_OIB = []
for i in range(0,len(kw_lons_final)):
    x = utm.from_latlon(kw_lats_final[i],kw_lons_final[i])
    easting_OIB.append(x[0])
    northing_OIB.append(x[1])

plt.figure(figsize=(12,7))
plt.scatter(easting_EY,northing_EY,c=snow_bc_aug2may_list,cmap="BuPu",vmin = 0, vmax=4)
plt.scatter(easting_OIB,northing_OIB,c=kw_snow_mwe, cmap="BuPu",marker='o',linewidth=3, vmin = 0, vmax=4) #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('CReSIS Snow Radar (May2021) Compared to \n Bias-Corrected NARR Accumulation (Aug2020-May2021)',fontsize=14)
#plt.savefig('NARR_CReSIS_2021_samecolours.png',bbox_inches = 'tight')

plt.figure(figsize=(13,7))
plt.scatter(easting_EY,northing_EY,c=elevation_EY,cmap="bone",vmin = 800, vmax=3000)
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.scatter(easting_OIB,northing_OIB,c=kw_snow_mwe, cmap="BuPu",marker='o',linewidth=3, vmin = 0, vmax=2) #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('Operation IceBridge Snow Radar (May2021)',fontsize=14)
plt.tight_layout()
#plt.savefig('OIBpicks_over_Zgrid.png',bbox_inches = 'tight')

highacc = np.where(np.array(kw_snow_mwe) > 1.5)
plt.figure(figsize=(13,7))
plt.scatter(easting_EY,northing_EY,c=elevation_EY,cmap="bone",vmin = 800, vmax=3000)
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.scatter(np.array(easting_OIB)[highacc],np.array(northing_OIB)[highacc],c=np.array(kw_snow_mwe)[highacc], cmap="BuPu",marker='o',linewidth=3, vmin = 0, vmax=2) #levels = np.linspace(0,4.5,19))
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.title('Operation IceBridge Snow Radar (May2021)',fontsize=14)
plt.tight_layout()

#################################################################################
#PLOT INDIVIDUAL TRIBUTARIES!!!!
#################################################################################

# segment the NASA data into the diff tributaries:
tribarray = KWtributaries()
nanlocs = np.where(np.isnan(Zgrid))
tribarray[nanlocs] = np.nan

def tributary_accvselev(tribarray,Zgrid,snowarray):
    SAs = []
    SAz = []
    SWs = []
    SWz = []
    NAs = []
    NAz = []
    CAs = []
    CAz = []
    Trs = []
    Trz = []
    
    for i in range(0,len(Zgrid)):
        for j in range(0,len(Zgrid[0])):
            if np.isnan(Zgrid[i][j]):
                pass
            else:
                if tribarray[i][j] == 1:
                    SAs.append(snowarray[i][j])
                    SAz.append(Zgrid[i][j])
                elif tribarray[i][j] == 2:
                    SWs.append(snowarray[i][j])
                    SWz.append(Zgrid[i][j])
                elif tribarray[i][j] == 3:
                    CAs.append(snowarray[i][j])
                    CAz.append(Zgrid[i][j])
                elif tribarray[i][j] == 4:
                    NAs.append(snowarray[i][j])
                    NAz.append(Zgrid[i][j])
                elif tribarray[i][j] == 5:
                    Trs.append(snowarray[i][j])
                    Trz.append(Zgrid[i][j])
                else:
                    pass
                
    return [SAs,SAz],[SWs,SWz],[CAs,CAz],[NAs,NAz],[Trs,Trz]

#get tributary data from NASA snow radar:
#s=12260
#e=-1 #12500
#plt.figure()
#plt.title('2021 KW Seasonal Snow Picks')
#plt.xlabel('Lon')
#plt.ylabel('Lat')
#plt.scatter(all_lons,all_lats,c=elevs,cmap="BuPu", vmin = 1000, vmax=3000)
#plt.scatter(kw_lons_final[s:e],kw_lats_final[s:e],c=kw_snow_mwe[s:e],cmap="BuPu", vmin = 0, vmax=5)
#plt.colorbar()

NAs_OIB = kw_snow_mwe[0:7845]
NAz_OIB = kw_elevs[0:7845]

#Trs_OIB = kw_snow_mwe[9000:31000] # no OIB data for trunk :/
#Trz_OIB = kw_elevs[9000:31000]

SAs_OIB = kw_snow_mwe[7850:12255]
SAz_OIB = kw_elevs[7850:12255]

CAs_OIB = kw_snow_mwe[12260:-1]
CAz_OIB = kw_elevs[12260:-1]
                    
plt.figure(figsize=(10,12))
plt.subplot(2,2,1)
plt.title('South Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[0][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[0][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[0][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[0][1],color='orange')
plt.scatter(SAs_OIB,SAz_OIB,color='turquoise')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.subplot(2,2,2)
plt.title('Stairway Glacier',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[1][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[1][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[1][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[1][1],color='orange')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.subplot(2,2,3)
plt.title('Central Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[2][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[2][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[2][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[2][1],color='orange')
plt.scatter(CAs_OIB,CAz_OIB,color='turquoise')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.subplot(2,2,4)
plt.title('North Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[3][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[3][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[3][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[3][1],color='orange')
plt.scatter(NAs_OIB,NAz_OIB,color='turquoise')
plt.xlabel('Accumulation (m w.e.)',fontsize=12)
plt.ylabel('Elevation (m a.s.l.)',fontsize=12)
plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.ylim(700,3700)
plt.xlim(0,4.5)
plt.tight_layout()
#plt.savefig('Tributaries_SnowvsZ.png',bbox_inches = 'tight')

###############################################################################
# PLOT MEAN ACCUMULATION VS ELEVATION LINE OVER EACH TRIBUTARY
##############################################################################

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
        elevation_bins.append(z+(0.5*delta_z))
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


#get average accumulations and elevation
avgsnow_noBC,zbins_noBC = meanaccumulation_vs_z(accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[1],accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[0],500,4000,delta_z=10)
avgsnow_BC,zbins_BC = meanaccumulation_vs_z(accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[1],accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[0],500,4000,delta_z=10)
avgsnow_OIB,zbins_OIB = meanaccumulation_vs_z(kw_elevs,kw_snow_mwe,1500,3000,delta_z=10)

plt.figure(figsize=(6,8))
plt.scatter(accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[1],marker='.',color='orange')
plt.scatter(accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta')
plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='turquoise')
plt.legend(['NARR Accumulation (Aug2020-May2021) Uncorrected','NARR Accumulation (Aug2020-May2021) Bias Corrected','CReSIS Snow Radar (May 2021)'])
plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod')
plt.plot(avgsnow_BC,zbins_BC,color='indigo')
plt.plot(avgsnow_OIB,zbins_OIB,color='darkcyan')
plt.title('Kaskawulsh Accumulation',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,4.5)
plt.ylim(500,4000)
#plt.savefig('NARRCresis2021_AccvsElev_summersnowinclv2.png',bbox_inches = 'tight')

plt.figure(figsize=(6,8))
plt.scatter(accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[0],accumulationvselevation(Zgrid,snow2021_nobc_plussummersnow)[1],marker='.',color='orange')
plt.scatter(accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[0],accumulationvselevation(Zgrid,snow2021_bc_plussummersnow)[1],marker='.',color='darkmagenta')
plt.scatter(kw_snow_mwe,kw_elevs,marker='.',color='turquoise')
plt.legend(['Uncorrected NARR','Original Bias Corrected NARR','OIB Snow Radar (May 2021)'],fontsize=15)
plt.plot(avgsnow_noBC,zbins_noBC,color='darkgoldenrod')
plt.plot(avgsnow_BC,zbins_BC,color='indigo')
plt.plot(avgsnow_OIB,zbins_OIB,color='darkcyan')
#plt.title('Kaskawulsh Accumulation',fontsize=12)
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.xlim(0,4.5)
plt.ylim(500,4000)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.savefig('OIBvsNARRvsElev.png',bbox_inches = 'tight')

#calculate averages for each trib
avgsnow_noBC_SA,zbins_noBC_SA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[0][1],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[0][0],500,4000,delta_z=10)
avgsnow_BC_SA,zbins_BC_SA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[0][1],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[0][0],500,4000,delta_z=10)
avgsnow_OIB_SA,zbins_OIB_SA = meanaccumulation_vs_z(SAz_OIB,SAs_OIB,500,4000,delta_z=10)

avgsnow_noBC_SW,zbins_noBC_SW = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[1][1],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[1][0],500,4000,delta_z=10)
avgsnow_BC_SW,zbins_BC_SW = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[1][1],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[1][0],500,4000,delta_z=10)
#avgsnow_OIB_SA,zbins_OIB_SA = meanaccumulation_vs_z(SAz_OIB,SAs_OIB,500,4000,delta_z=10)

avgsnow_noBC_CA,zbins_noBC_CA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[2][1],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[2][0],500,4000,delta_z=10)
avgsnow_BC_CA,zbins_BC_CA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[2][1],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[2][0],500,4000,delta_z=10)
avgsnow_OIB_CA,zbins_OIB_CA = meanaccumulation_vs_z(CAz_OIB,CAs_OIB,500,4000,delta_z=10)

avgsnow_noBC_NA,zbins_noBC_NA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[3][1],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[3][0],500,4000,delta_z=10)
avgsnow_BC_NA,zbins_BC_NA = meanaccumulation_vs_z(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[3][1],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[3][0],500,4000,delta_z=10)
avgsnow_OIB_NA,zbins_OIB_NA = meanaccumulation_vs_z(NAz_OIB,NAs_OIB,500,4000,delta_z=10)

plt.figure(figsize=(10,12))
plt.subplot(2,2,1)
plt.title('South Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[0][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[0][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[0][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[0][1],color='orange')
plt.scatter(SAs_OIB,SAz_OIB,color='turquoise')
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
#plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.plot(avgsnow_noBC_SA,zbins_noBC_SA,color='darkgoldenrod')
plt.plot(avgsnow_BC_SA,zbins_BC_SA,color='indigo')
plt.plot(avgsnow_OIB_SA,zbins_OIB_SA,color='darkcyan')
plt.ylim(1000,3700)
plt.xlim(0,4.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.subplot(2,2,2)
plt.title('Stairway Glacier',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[1][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[1][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[1][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[1][1],color='orange')
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
#plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.plot(avgsnow_noBC_SW,zbins_noBC_SW,color='darkgoldenrod')
plt.plot(avgsnow_BC_SW,zbins_BC_SW,color='indigo')
plt.ylim(1000,3700)
plt.xlim(0,4.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.subplot(2,2,3)
plt.title('Central Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[2][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[2][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[2][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[2][1],color='orange')
plt.scatter(CAs_OIB,CAz_OIB,color='turquoise')
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
#plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.plot(avgsnow_noBC_CA,zbins_noBC_CA,color='darkgoldenrod')
plt.plot(avgsnow_BC_CA,zbins_BC_CA,color='indigo')
plt.plot(avgsnow_OIB_CA,zbins_OIB_CA,color='darkcyan')
plt.ylim(1000,3700)
plt.xlim(0,4.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.subplot(2,2,4)
plt.title('North Arm',fontsize=14)
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[3][0],tributary_accvselev(tribarray,Zgrid,snow2021_bc_plussummersnow)[3][1],color='darkmagenta')
plt.scatter(tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[3][0],tributary_accvselev(tribarray,Zgrid,snow2021_nobc_plussummersnow)[3][1],color='orange')
plt.scatter(NAs_OIB,NAz_OIB,color='turquoise')
plt.xlabel('Accumulation (m w.e.)',fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
#plt.legend(['NARR (Aug2020-May2021) Bias Corrected','NARR (Aug2020-May2021) Uncorrected','CReSIS Snow Radar (2021)'])
plt.plot(avgsnow_noBC_NA,zbins_noBC_NA,color='darkgoldenrod')
plt.plot(avgsnow_BC_NA,zbins_BC_NA,color='indigo')
plt.plot(avgsnow_OIB_NA,zbins_OIB_NA,color='darkcyan')
plt.ylim(1000,3700)
plt.xlim(0,4.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
#plt.savefig('Tributaries_SnowvsZ_v2.png',bbox_inches = 'tight')
plt.savefig('Tributaries_SnowvsZ_v3.png',bbox_inches = 'tight')

###############################################################################
# CALCULATE RMSE BASED ON CO-LOCATED GRIDCELLS, NOT BY BINNED ELEVATIONS
    # bin OIB cells into nearest NARR gridcell. Save into an array of shape Zgrid 
        #(in the corresponding cell)
###############################################################################

#for each point in the OIB data, find the closest model gridcell
# gridcell must be within np.sqrt(100**2 + 100**2) = 141.42 m for it to be
# within a NARR gridcell. if the closest cell is > 141.42m then its outside of the model domain.

# add OIB snowdpeth to that gridcell it corresponds to,
#also track the number of OIB gridcells that correspond to each NARR gridcell
# use these two arrays to calculate mean snowdepth.
Xgrid[nanlocs] = np.nan
Ygrid_flipped[nanlocs] = np.nan
mindist = []
Colocated_snow_sum = np.zeros(Zgrid.shape)
Colocated_numberofcells = np.zeros(Zgrid.shape) #number of OIB cells corresponding to each NARR gridcell
for i in range(0,len(kw_snow_mwe)):
    x_OIB = easting_OIB[i]
    y_OIB = northing_OIB[i]
    OIB_snow = kw_snow_mwe[i]
    
    x_dists = Xgrid - x_OIB
    y_dists = Ygrid_flipped - y_OIB
    distance = np.sqrt((x_dists**2)+(y_dists**2))
    if np.nanmin(distance) <= np.sqrt(100**2 + 100**2):
        mindist.append(np.nanmin(distance))
        closest_cell = np.where(distance == np.nanmin(distance))
        Colocated_snow_sum[closest_cell] += OIB_snow #add the snow to the corresponding gridcell
        Colocated_numberofcells[closest_cell] += 1
    else:
        pass
Colocated_avg_OIBsnow = Colocated_snow_sum/Colocated_numberofcells
    

Colocated_avg_OIBsnow[nanlocs] = np.nan   
#yields 573 colocated cells! 

plt.figure()
plt.contourf(np.flipud(Colocated_avg_OIBsnow))
plt.colorbar()









