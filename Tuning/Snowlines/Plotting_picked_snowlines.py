# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 17:06:32 2023

A script for plotting picked snowlines on KW

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
import pandas as pd
import sys
from netCDF4 import Dataset
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import KWtributaries

# load in one file at a time (loop through the list of images that are completed)
# change this: load the list of files within each directory in the main dorectory:
# if tehe folder is verified,, open it and get the csv snowlines out.

filepath = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Snowlines/Snowlines_FINAL'
#get list of files in this directory:

def snowlines_from_csv(filepath):
    '''
    filepath is the path the the main directory where all the snowlines are stored.
    '''
    # set up containers to hold outputs from each snowline
    snowline_name = []
    easting = []
    northing = []
    elevation = []
    sfc_type = []
    
    folderlist = os.listdir(filepath)
    for folder in folderlist:
        if folder[-8:] == 'Verified':
            # go inside that folder 
            folderpath = os.path.join(filepath,folder)
            # get all csv files (# of csv's?)
            filelist = os.listdir(folderpath)
            # for each of the csv's, save the name, x, y, z, sfc type. 
            for file in filelist:
                if '.csv' in file:
                    print(file)
                    snowline_name.append(file[:-4])
                    
                    csvfile = os.path.join(folderpath,file)
                    snowlinedata = np.loadtxt(csvfile,delimiter=',',skiprows=1,dtype='str')
                    x = np.array(snowlinedata[:,0],dtype=float)
                    easting.append(x)
                    
                    y = np.array(snowlinedata[:,1],dtype=float)
                    northing.append(y)
                    
                    z = np.array(snowlinedata[:,2],dtype=float)
                    elevation.append(z)
                    
                    sfc = np.array(snowlinedata[:,3])
                    sfc_type.append(sfc)
                    # add line to delete off-glacier cells from all arrays
                
            
        else:
            pass
        
    return np.array(snowline_name), easting, northing, elevation, sfc_type
    
snowline_name, easting, northing, elevation, snowline_sfc_type = snowlines_from_csv(filepath)

# get dates in datetime format for plotting:

# start plotting:
File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kask_catchment.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,4] 
Iy = glacier[:,5] 
Ih = glacier[:,6]  
sfc_type = glacier[:,8]
ELA = glacier[:,10]     

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
Ygrid = np.flipud(Ygrid)

def plot_KW_Zgrid():
    plt.contourf(Xgrid,Ygrid,Zgrid,cmap = 'Greys_r',levels=np.linspace(600,4000,35))
    legend = plt.colorbar()
    plt.axis('equal') 
    legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
    plt.xlabel('Easting',fontsize=14)
    plt.ylabel('Northing',fontsize=14)
    plt.show

plt.figure(figsize=(10,6))
plt.title('All Snowlines',fontsize=14)
plot_KW_Zgrid()
colors = plt.cm.Blues(np.linspace(0,1,len(snowline_name)))
for i in range(0,len(snowline_name)):
    plt.scatter(easting[i],northing[i],s=3,label=snowline_name[i],color=colors[i])
plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))
plt.savefig('allsnowlines.png',bbox_inches='tight')


def get_subset_of_snowlines(years,trib,UL,month='Any'):
    '''
    run snowlines_from_csv function first to get list of snowlines, easting, northing, etc.
    year = desired year(s)
    trib = (string) desired tributary: 'NA', 'CA', 'SW', 'SA'
    U/L = Upper (U) or Lower (L) snowline: 'U', 'L'
    '''
    snowline_subset = []
    x = []
    y = []
    z = []
    sfc = []
    date = []
    
    if month == 'Any':
        for year in years:
            #print(y)
            for n in snowline_name:
                if n[:4] == str(year):
                    #print(n)
                    if n[11:13] == trib:
                        #print(n)
                        if n[14] == UL:
                            #print(n)
                            index = np.where(snowline_name == n)[0][0]
                            snowline_subset.append(n)
                            x.append(easting[index])
                            y.append(northing[index])
                            z.append(elevation[index])
                            sfc.append(snowline_sfc_type[index])
                            date.append(n[0:10])
                            
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
    else:
        for year in years:
            for n in snowline_name:
                if n[:4] == str(year):
                    #print(n)
                    if int(n[5:7]) == month:
                        if n[11:13] == trib:
                            #print(n)
                            if n[14] == UL:
                                #print(n)
                                index = np.where(snowline_name == n)[0][0]
                                snowline_subset.append(n)
                                x.append(easting[index])
                                y.append(northing[index])
                                z.append(elevation[index])
                                sfc.append(snowline_sfc_type[index])
                                date.append(n[0:10])
                                
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
        
    return snowline_subset, x, y, z, date, sfc
    
# Plotting all SA_L snowlines
snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'SA','L')
plt.figure(figsize=(10,6))
plt.title('South Arm Upper Snowlines (2013-2019)',fontsize=14)
plot_KW_Zgrid()
colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
for i in range(0,len(snowlines_subset[0])):
    plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[i])
plt.legend(ncol=3,bbox_to_anchor=(2.0, 1.15))

# Plotting just the snowlines for a specific tributary, month and year (July 2018)
snowlines_subset = get_subset_of_snowlines([2018],'SA','L',7)
plt.figure(figsize=(10,6))
plt.title('South Arm Lower Snowlines (2018)',fontsize=14)
plot_KW_Zgrid()
colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
for i in range(0,len(snowlines_subset[0])):
    plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[i])
plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))

# Plotting all the upper bounds for a specific month in any year:
plt.figure(figsize=(10,6))
plt.title('August Upper Snowlines (2013--2019)',fontsize=14)
plot_KW_Zgrid()
snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'SA','U',8)
colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
for i in range(0,len(snowlines_subset[0])):
    plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[i])
plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))
snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'SW','U',8)
colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
for i in range(0,len(snowlines_subset[0])):
    plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[i])
plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))
snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'NA','U',8)
colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
for i in range(0,len(snowlines_subset[0])):
    plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[i])
plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))
snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'CA','U',8)
colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
for i in range(0,len(snowlines_subset[0])):
    plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[i])
plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))


plt.figure(figsize=(10,6))
plt.title('Upper Snowlines (2013--2019)',fontsize=14)
plot_KW_Zgrid()
colors = ['red','gold','lime','turquoise','blue','blueviolet']
for month in range(5,11):
    snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'SA','U',month)
    #colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
    for i in range(0,len(snowlines_subset[0])):
        plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[month-5])
    plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))
    snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'SW','U',month)
    #colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
    for i in range(0,len(snowlines_subset[0])):
        plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[month-5])
    plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))
    snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'NA','U',month)
    #colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
    for i in range(0,len(snowlines_subset[0])):
        plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[month-5])
    plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))
    snowlines_subset = get_subset_of_snowlines(np.arange(2007,2020),'CA','U',month)
    #colors = plt.cm.Blues(np.linspace(0,0.5,len(snowlines_subset[0])))
    for i in range(0,len(snowlines_subset[0])):
        plt.scatter(snowlines_subset[1][i],snowlines_subset[2][i],s=3,label=snowlines_subset[0][i],color=colors[month-5])
    plt.legend(ncol=3,bbox_to_anchor=(1.25, 1.15))

def mean_snowline(snowlines_subset):
# Get mean snowline elevation:
    snowline_names = snowlines_subset[0]
    elevation = snowlines_subset[3]
    mean_snowline_z = []
    dates = []
    sfc_type = []
    
    # get sfc type bool:
    for sfc in snowlines_subset[5]:
        sfc_list = []
        for i in sfc:
            if len(i) > 1:
                sfc_list.append(int(i[1]))
            else:
                sfc_list.append(int(i))
        sfc_type.append(sfc_list)
    
    i = 0
    for z in elevation:
        zz = np.array(z)
        offglacier = np.where(np.array(sfc_type[i][0])==1)
        zz[offglacier] = np.nan
        mean_z = np.nanmean(zz)
        mean_snowline_z.append(mean_z)
        i += 1
        
    for date in snowline_names:
        dates.append(np.datetime64(date[:10]))
        
    return snowline_names, dates, mean_snowline_z

z_SW_L = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'SW','L'))
z_SW_U = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'SW','U'))

z_SA_L = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'SA','L'))
z_SA_U = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'SA','U'))

z_CA_L = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'CA','L'))
z_CA_U = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'CA','U'))

z_NA_L = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'NA','L'))
z_NA_U = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'NA','U'))

z_TR_L = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'TR','L'))
z_TR_U = mean_snowline(get_subset_of_snowlines(np.arange(2013,2020),'TR','U'))

def combined_dates_barplot(z_U,z_L,AllDates=False):
    if AllDates == False:
        dates = list( dict.fromkeys(sorted(z_U[1]+z_L[1])) )
        width = 0.8
    else:
        dates = list( dict.fromkeys(sorted(z_SW_L[1]+z_SW_U[1]+z_SA_L[1]+z_SA_U[1]+z_CA_L[1]+z_CA_U[1]+z_NA_L[1]+z_NA_U[1]+z_TR_L[1]+z_TR_U[1])) )
        width = 0.7
    
    new_z_L = np.zeros(len(dates))
    new_z_U = np.zeros(len(dates))
    for i in range(0,len(dates)):
        date = dates[i]
        Upper_ind = np.where(z_U[1]==date)[0]
        if len(Upper_ind) == 0:
            pass
        else:
            new_z_U[i] = z_U[2][Upper_ind[0]]
            
        Lower_ind = np.where(z_L[1]==date)[0]
        if len(Lower_ind) == 0:
            pass
        else:
            new_z_L[i] = z_L[2][Lower_ind[0]]
    
    plt.bar(np.arange(0,len(dates)),new_z_U,width=width,label='Upper Boundary',color='lightskyblue')
    plt.bar(np.arange(0,len(dates)),new_z_L,width=width,label='Lower Boundary',color='grey')
    plt.xlim(-1,len(dates))
    plt.ylim(650,2600)
    plt.xticks(ticks=np.arange(0,len(dates),4),labels=dates[::4],rotation=45)
    plt.ylabel('Elevation (m a.s.l.)',fontsize=13)
    #plt.plot([dates[0],dates[-1]],[2000,2000])

fig = plt.figure(figsize=(12,8))
plt.subplot(2,3,1)
plt.title('North Arm',fontsize=14)
combined_dates_barplot(z_NA_U,z_NA_L)
fig.legend(bbox_to_anchor=(0.9,0.4),fontsize=15)
plt.subplot(2,3,2)
plt.title('Central Arm',fontsize=14)
combined_dates_barplot(z_CA_U,z_CA_L)
plt.subplot(2,3,3)
plt.title('Stairway Glacier',fontsize=14)
combined_dates_barplot(z_SW_U,z_SW_L)
plt.subplot(2,3,4)
plt.title('South Arm',fontsize=14)
combined_dates_barplot(z_SA_U,z_SA_L)
plt.subplot(2,3,5)
plt.title('Trunk',fontsize=14)
combined_dates_barplot(z_TR_U,z_TR_L)
plt.tight_layout()
plt.savefig('Snowlines_histogram_combined_dates.png')

fig = plt.figure(figsize=(12,8))
plt.subplot(2,3,1)
plt.title('North Arm',fontsize=14)
combined_dates_barplot(z_NA_U,z_NA_L,AllDates=True)
fig.legend(bbox_to_anchor=(0.9,0.4),fontsize=15)
plt.subplot(2,3,2)
plt.title('Central Arm',fontsize=14)
combined_dates_barplot(z_CA_U,z_CA_L,AllDates=True)
plt.subplot(2,3,3)
plt.title('Stairway Glacier',fontsize=14)
combined_dates_barplot(z_SW_U,z_SW_L,AllDates=True)
plt.subplot(2,3,4)
plt.title('South Arm',fontsize=14)
combined_dates_barplot(z_SA_U,z_SA_L,AllDates=True)
plt.subplot(2,3,5)
plt.title('Trunk',fontsize=14)
combined_dates_barplot(z_TR_U,z_TR_L,AllDates=True)
plt.tight_layout()
plt.savefig('Snowlines_histogram_all_dates.png')


# Calculate mean monthly snowline
# get mean monthly elev:
def mean_monthly_snowline(meansnowlines_list):
    mean_monthly_elev = []
    mean_monthly_std = []
    for m in range(3,12):
        mean_sl_elevs = []
        for d in range(0,len(meansnowlines_list[1])):
            date = meansnowlines_list[1][d]
            if date.astype(object).month == m:
                mean_sl_elevs.append(meansnowlines_list[2][d])
            else:
                pass
        
        mean_monthly_elev.append(np.mean(mean_sl_elevs))
        mean_monthly_std.append(np.std(mean_sl_elevs))

    return mean_monthly_elev, mean_monthly_std

plt.figure(figsize=(8,6))
plt.title('Mean (2013--2019) Monthly Snowline Elevation',fontsize=14)
x = np.arange(1,10)
width = 0.16
plt.bar(x-0.32,mean_monthly_snowline(z_TR_U)[0],width=width,color='mediumpurple')
plt.bar(x-0.32,mean_monthly_snowline(z_TR_L)[0],width=width,color='rebeccapurple',label='Trunk')
plt.bar(x-0.16,mean_monthly_snowline(z_SA_U)[0],width=width,color='lightskyblue')
plt.bar(x-0.16,mean_monthly_snowline(z_SA_L)[0],width=width,color='steelblue',label='South Arm')
plt.bar(x,mean_monthly_snowline(z_SW_U)[0],width=width,color='pink')
plt.bar(x,mean_monthly_snowline(z_SW_L)[0],width=width,color='palevioletred',label='Stairway')
plt.bar(x+0.16,mean_monthly_snowline(z_CA_U)[0],width=width,color='gold')
plt.bar(x+0.16,mean_monthly_snowline(z_CA_L)[0],width=width,color='goldenrod',label='Central Arm')
plt.bar(x+0.32,mean_monthly_snowline(z_NA_U)[0],width=width,color='mediumseagreen')
plt.bar(x+0.32,mean_monthly_snowline(z_NA_L)[0],width=width,color='green',label='North Arm')
plt.xticks(x,labels=['Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov'],fontsize=14,rotation=45)
plt.legend(fontsize=14,bbox_to_anchor=(1.05,1))
plt.ylim(650,2600)
plt.yticks(fontsize=14)
plt.ylabel('Elevation (m a.s.l.)',fontsize=14)
plt.savefig('Mean_monthly_snowline_elevation.png',bbox_inches='tight')




# Plotting model trial runs below:

Path2KWoutline = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel'
File_glacier_in = os.path.join(Path2KWoutline,'kaskonly_deb.txt')
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,2] 
balancefluxes = glacier[:,5]
debris_array = glacier[:,6]

#-------Turn vectors into 3D gridded inputs--------------------
Zgrid_KW, Xgrid_KW, Ygrid_KW, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)

tribarray = KWtributaries()
nanlocs_KW = np.where(np.isnan(Zgrid_KW))
tribarray[nanlocs_KW] = np.nan

plt.figure(figsize=(8,6))
plt.contourf(Xgrid_KW,np.flipud(Ygrid_KW),tribarray,cmap='Accent',levels=np.arange(0,6))
#plt.colorbar()
plt.axis('equal')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
plt.tight_layout()
#plt.savefig('KW_tributaries.pdf')

#    tributary codes:
#    SA = 1
#    SW = 2
#    CA = 3
#    NA = 4
#    Trunk = 5

def rasterize_observed_snowline(date):
    '''
    function to rasterize a satellite image with delineated snowlines.
    Returns raster that can be used to compare with modelled snowlines.
    Key for output raster:
    2 = snow covered
    1 = transition zone
    0 = ice
    nan = off glacier
    -0.5 = NO INFORMATION - DO NOT COUNT IN MOD VS OBS COMPARISON
    '''
    if len(np.where(snowline_name == date + '_NA_U')[0]) == 0:
        NAu = np.nan
    else:
        NAu = np.nanmean(elevation[np.where(snowline_name == date + '_NA_U')[0][0]])
    
    if len(np.where(snowline_name == date + '_NA_L')[0]) == 0:
        NAl = np.nan
    else:
        NAl = np.nanmean(elevation[np.where(snowline_name == date + '_NA_L')[0][0]])
    
    if len(np.where(snowline_name == date + '_CA_U')[0]) == 0:
        CAu = np.nan
    else:
        CAu = np.nanmean(elevation[np.where(snowline_name == date + '_CA_U')[0][0]])
    
    if len(np.where(snowline_name == date + '_CA_L')[0]) == 0:
        CAl = np.nan
    else:
        CAl = np.nanmean(elevation[np.where(snowline_name == date + '_CA_L')[0][0]])
    
    if len(np.where(snowline_name == date + '_SW_U')[0]) == 0:
        SWu = np.nan
    else:
        SWu = np.nanmean(elevation[np.where(snowline_name == date + '_SW_U')[0][0]])
        
    if len(np.where(snowline_name == date + '_SW_L')[0]) == 0:
        SWl = np.nan
    else:
        SWl = np.nanmean(elevation[np.where(snowline_name == date + '_SW_L')[0][0]])
    
    if len(np.where(snowline_name == date + '_SA_U')[0]) == 0:
        SAu = np.nan
    else:
        SAu = np.nanmean(elevation[np.where(snowline_name == date + '_SA_U')[0][0]])
        
    if len(np.where(snowline_name == date + '_SA_L')[0]) == 0:
        SAl = np.nan
    else:
        SAl = np.nanmean(elevation[np.where(snowline_name == date + '_SA_L')[0][0]])
        
    if len(np.where(snowline_name == date + '_TR_U')[0]) == 0:
        TRu = np.nan
        TRu_min = np.nan
    else:
        TRu = np.nanmean(elevation[np.where(snowline_name == date + '_TR_U')[0][0]])
        TRu_min = np.nanmin(elevation[np.where(snowline_name == date + '_TR_U')[0][0]])
        
    if len(np.where(snowline_name == date + '_TR_L')[0]) == 0:
        TRl = np.nan
        TRl_min = np.nan
    else:
        TRl = np.nanmean(elevation[np.where(snowline_name == date + '_TR_L')[0][0]])
        TRl_min = np.nanmin(elevation[np.where(snowline_name == date + '_TR_L')[0][0]])
    
    obs_snowline = np.empty(Zgrid_KW.shape)
    obs_snowline[:] = -0.5
    # 2 = snow covered
    # 1 = transition zone
    # 0 = ice
    # nan = off glacier
    # -0.5 = NO INFORMATION - DO NOT COUNT IN MOD VS OBS COMPARISON

    for i in range(0,len(Zgrid_KW)): # For every cell in the domain
        for j in range(0,len(Zgrid_KW[1])):
            if np.isnan(tribarray[i][j]): # If the cell is off glacier, skip.
                pass
            
            elif tribarray[i][j] == 1: #If the cell is in South Arm:
                upperbound = SAu
                if np.isnan(upperbound): # If there's no upper bound on SA but there is on TR:
                    if ~np.isnan(TRu): 
                        if (TRu_min < np.min(Zgrid_KW[tribarray==1])): # And if the min elev of TR_U is < than SA:
                            upperbound = TRu_min  # Set SA_U = TR_U
                    else:
                        pass
                else:
                    pass
                
                
                lowerbound = SAl
                if np.isnan(lowerbound):
                    if ~np.isnan(TRl):
                        if (TRl < np.min(Zgrid_KW[tribarray==1])):
                            lowerbound = TRl
                    else:
                        pass
                else:
                    pass                
                
            elif tribarray[i][j] == 2: #SW
                lowerbound = SWl
                upperbound = SWu
                # special case: no upper bound on SW, but there is an upper bound on TR.
                # Set SW_U = TR_U
                if np.isnan(upperbound):
                    if ~np.isnan(TRu):
                        if (TRu < np.min(Zgrid_KW[tribarray==2])):
                            upperbound = TRu
                    else:
                        pass
                else:
                    pass
                # special case: no lower bound on SW, but there is an lower bound on TR.
                if np.isnan(lowerbound):
                    if ~np.isnan(TRl):
                        if (TRl < np.min(Zgrid_KW[tribarray==2])):
                            lowerbound = TRl
                    else:
                        pass
                else:
                    pass                     
                
            elif tribarray[i][j] == 3: #CA
                lowerbound = CAl
                upperbound = CAu 
                # special case: no upper bound on CA, but there is an upper bound on TR.
                # Set CA_U = TR_U
                if np.isnan(upperbound):
                    if ~np.isnan(TRu):
                        upperbound = TRu
                    else:
                        pass
                else:
                    pass
            elif tribarray[i][j] == 4: #NA
                lowerbound = NAl
                upperbound = NAu 
                # special case: no upper bound on NA, but there is an upper bound on TR.
                # Set CA_U = TR_U
                if np.isnan(upperbound):
                    if ~np.isnan(TRu):
                        upperbound = TRu
                    else:
                        pass
                else:
                    pass
                # special case: no lower bound on SA, but there is an lower bound on TR.
                if np.isnan(lowerbound):
                    if ~np.isnan(TRl):
                        if (TRl < np.min(Zgrid_KW[tribarray==4])):
                            lowerbound = TRl
                    else:
                        pass
                else:
                    pass                     
                
                
            elif tribarray[i][j] == 5: # TRUNK
                lowerbound = TRl
                upperbound = TRu
                # special case for if there is no upper or lower bound for trunk:
                if np.isnan(lowerbound): 
                    if np.isnan(upperbound):
                        # if there is a lower bound of CA or NA, assume that LB on trunk is the min of the CA or NA lower bounds.
                        if ~np.isnan(CAl):
                            lowerbound = np.min([CAl,NAl])
                        elif ~np.isnan(NAl):
                            lowerbound = np.min([CAl,NAl])
                    else:
                        pass
                else:
                    pass
                
                # another special case:
                # if there is an upper bound on the trunk, any trib without an upper bound can assume that
                # any area above the upper bound of the trunk is also snow covered.
             
            # determine if the cell is above the upper bound (snow = 2), below the lower bound (ice = 0), or in between (transition = 1)
            z = Zgrid_KW[i][j]
            if ~np.isnan(z):
                # if cell is above upper bound, set to snow
                if ~np.isnan(upperbound):
                    if z >= upperbound:
                        obs_snowline[i][j] = 2
                    else:
                        pass
                else:
                    pass
                # if cell is below lower bound, set to ice
                if ~np.isnan(lowerbound):
                    if z <= lowerbound:
                        obs_snowline[i][j] = 0
                    else:
                        pass
                else:
                    pass
                
                # if cell is between upper and lower bounds, set to transition zone
                # there is only a transition zone if BOTH upper and lower bounds exist
                if ~np.isnan(upperbound):
                    if ~np.isnan(lowerbound):
                        if z > lowerbound:
                            if z < upperbound:
                                obs_snowline[i][j] = 1
                            
    obs_snowline[nanlocs_KW] = np.nan
        
    plt.figure(figsize=(8,6))
    plt.pcolor(Xgrid_KW,Ygrid_KW,np.flipud(obs_snowline),cmap='terrain_r',vmin=-0.75,vmax=2.6)
    #plt.colorbar()
    plt.axis('equal')
    plt.xlabel('Easting',fontsize=14)
    plt.ylabel('Northing',fontsize=14)
    plt.title('Observed Snowline: ' + date,fontsize=14)
    plt.tight_layout()
    plt.savefig('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Snowlines/Rasterized_Observed_Snowlines/SnowlineTarget_' + date + '.png')
    plt.show()
    plt.close()

    
    return obs_snowline

dates = list( dict.fromkeys(sorted(z_SW_L[1]+z_SW_U[1]+z_SA_L[1]+z_SA_U[1]+z_CA_L[1]+z_CA_U[1]+z_NA_L[1]+z_NA_U[1]+z_TR_L[1]+z_TR_U[1])) )
for day in dates:
    rasterize_observed_snowline(str(day))
    sys.stdout.flush()





# Get modelled snowline for same date
File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kaskonly.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,5]   
ELA = glacier[:,9]
TOPO = glacier[:,7]    
        
Zgrid1, Xgrid1, Ygrid1, xbounds1, ybounds1 = regridXY_something(Ix, Iy, Ih)
nanlocs_1 = np.where(np.isnan(Zgrid1))

for sim in range(0,12):
    MB_2015_file = 'D:/Model Runs/Debris Tests/VariableThicknessTest_2007-2018/MB2015' + str(sim) + '.nc'
    MB_2016_file = 'D:/Model Runs/Debris Tests/VariableThicknessTest_2007-2018/MB2016' + str(sim) + '.nc'
    MB_in_2015 = Dataset(MB_2015_file, "r")
    MB_in_2016 = Dataset(MB_2016_file, "r")
    
    MB_var = 'MB'
    
    MB_2015 = MB_in_2015.variables[MB_var][:]
    MB_2016 = MB_in_2016.variables[MB_var][:]
    
    dates_2015 = pd.date_range(start="2015-01-01",end="2015-12-31 21:00:00",freq='3H')
    dates_2016 = pd.date_range(start="2016-01-01",end="2016-12-31 21:00:00",freq='3H')
    
    snow_depth = np.zeros(MB_2015[0].shape)
    snow_depth[nanlocs_1] = np.nan
    
    start_date = np.where(dates_2015 == '2015-09-01')[0][0]
    end_date = np.where(dates_2016 == '2016-07-17')[0][0]
    
    for i in range(start_date,len(dates_2015)):
        snow_depth += MB_2015[i]
        snow_depth[np.where(snow_depth < 0)] = 0
    for i in range(0,end_date):
        snow_depth += MB_2016[i]
        snow_depth[np.where(snow_depth < 0)] = 0

    #plt.figure(figsize=(9,6))
    #plt.contourf(Xgrid1,Ygrid1,np.flipud(snow_depth),cmap='BuPu')
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Snow Depth (m w.e.)', rotation=270,fontsize=14,labelpad=20)
    #plt.axis('equal')
    #plt.contour(np.flipud(snow_depth), colors = 'black', levels = 0.1)
    #plt.xlabel('Easting',fontsize=14)
    #plt.ylabel('Northing',fontsize=14)
    #plt.title('Modelled Snow Depth\n2015-09-01 - 2016-07-17',fontsize=14)
    #plt.savefig('ModelledSnowdepth_2016-07-17.png')
    
    modelled_snowline = np.zeros(MB_2015[0].shape)
    
    modelled_snowline[np.where(snow_depth <=0)] = 0
    modelled_snowline[np.where(snow_depth >0)] = 2
    modelled_snowline[nanlocs_1] = np.nan
    
    #plt.figure(figsize=(8,6))
    #plt.pcolor(Xgrid_KW,Ygrid_KW,np.flipud(modelled_snowline),cmap='terrain_r',vmin=-0.75,vmax=2.6)
    #plt.colorbar()
    #plt.axis('equal')
    #plt.xlabel('Easting',fontsize=14)
    #plt.ylabel('Northing',fontsize=14)
    #plt.title('Modelled Snowline: 2016-07-17',fontsize=14)
    #plt.savefig('Modelled_Snowline_2016-07-17.png')
    
    # Compare % of matching cells:
    diff = obs_snowline - modelled_snowline
    passing_cells = (len(np.where(diff == -1)[0]) + len(np.where(diff == 1)[0]) + len(np.where(diff == 0)[0]))
    percent = (passing_cells/26314)*100
    print(percent)
