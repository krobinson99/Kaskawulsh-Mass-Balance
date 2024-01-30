# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 20:37:05 2023

PLotting historical weather data from stations in the Yukon

@author: katierobinson
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
import sys
from scipy.interpolate import interp1d

sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import save_to_netcdf

# Input geometry
Easting_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Xgrid.txt' # Paths to text files defining Easting/Northing coords of every model gridcell
Northing_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/KRH_Ygrid.txt'
Sfc_grid = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_SfcType.txt'      # Path to text file where 1 = off-glacier, 0 = on-glacier, NaN = not in the domain.

Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
Sfc = np.loadtxt(Sfc_grid)

Catchmentoutline = np.array(Sfc)
Catchmentoutline[np.where(np.isfinite(Sfc))] = 0
Catchmentoutline[np.where(np.isnan(Sfc))] = 1

KRH_tributaries = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Tributaries.txt')

All_glacierized_area = np.array(Sfc)
All_glacierized_area[np.where(Sfc==1)] = np.nan


# Load the CSV file into a Pandas DataFrame
Kluanelake = pd.read_csv('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100680_1970-1983_P1M.csv')

# Display the contents of the DataFrame
print(Kluanelake['Total Precip (mm)'])

year_i, year_f = np.array(Kluanelake['Year'])[0], np.array(Kluanelake['Year'])[-1]
KluaneLake_precip = np.array(Kluanelake['Total Precip (mm)']).reshape(((year_f-year_i)+1,12))
KluaneLake_year = np.array(Kluanelake['Year']).reshape(((year_f-year_i)+1,12))
KluaneLake_month = np.array(Kluanelake['Month']).reshape(((year_f-year_i)+1,12))

dates = pd.date_range(start= str(year_i) + '-01-01',end= str(year_f) + '-12-01',freq='MS')

KluaneLake_precip = np.array(Kluanelake['Total Precip (mm)'])
#KluaneLake_year = np.array(Kluanelake['Year'])
#KluaneLake_month = np.array(Kluanelake['Month']).reshape(((year_f-year_i)+1,12))

def Station_data(station_path):
    Kluanelake = pd.read_csv(station_path)
    year_i, year_f = np.array(Kluanelake['Year'])[0], np.array(Kluanelake['Year'])[-1]
    month_i, month_f = np.array(Kluanelake['Month'])[0], np.array(Kluanelake['Month'])[-1]
    dates = pd.date_range(start= str(year_i) + '-' + str(month_i) + '-01',end= str(year_f) + '-' + str(month_f) + '-01',freq='MS')

    print(month_i,month_f)
    KluaneLake_precip = np.array(Kluanelake['Total Precip (mm)'])
    
    plt.title(Kluanelake['Station Name'][0] + ' (' + str(Kluanelake['Latitude (y)'][0]) + ', ' + str(Kluanelake['Longitude (x)'][0]) + ')' , fontsize=14)
    plt.plot(dates,KluaneLake_precip)
    plt.xticks(fontsize=14,rotation=25)
    plt.ylabel('Monthly\nPrecipitation (mm)',fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0,200)
    plt.grid()
    
    return dates

    
plt.figure(figsize=(20,5))
plt.subplot(2,3,1)
Station_data('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100680_1970-1983_P1M.csv')
plt.subplot(2,3,2)
Station_data('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100418_1975-1984_P1M.csv')
plt.subplot(2,3,3)
Station_data('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100182_1966-2007_P1M.csv')
plt.subplot(2,3,4)
Station_data('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100914_1983-1985_P1M.csv')
plt.xticks(ticks=qc_dates[np.arange(0,36,12)],labels=['1983','1984','1985'])
plt.subplot(2,3,5)
Station_data('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100630_1944-2006_P1M.csv')
plt.subplot(2,3,6)
Station_data('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100631_1985-2007_P1M.csv')
plt.tight_layout()

def monthly_mean(station_path):
    Kluanelake = pd.read_csv(station_path)
    year_i, year_f = np.array(Kluanelake['Year'])[0], np.array(Kluanelake['Year'])[-1]
    month_i, month_f = np.array(Kluanelake['Month'])[0], np.array(Kluanelake['Month'])[-1]
    if month_f == 12:
        KluaneLake_year_og = np.array(Kluanelake['Year']).reshape(((year_f-year_i)+1,12))
        KluaneLake_year = np.array(Kluanelake['Year']).reshape(((year_f-year_i)+1,12))[np.where(KluaneLake_year_og>=1900)[0][0]:]
        #KluaneLake_month = np.array(Kluanelake['Month']).reshape(((year_f-year_i)+1,12))[np.where(KluaneLake_year_og>=1979)[0][0]:]
        KluaneLake_precip = np.array(Kluanelake['Total Precip (mm)']).reshape(((year_f-year_i)+1,12))[np.where(KluaneLake_year_og>=1900)[0][0]:]
    else:
        KluaneLake_precip_filled = np.concatenate((np.array(Kluanelake['Total Precip (mm)']),np.ones((12-month_f))*np.nan))
        KluaneLake_precip = KluaneLake_precip_filled.reshape(((year_f-year_i)+1,12))
    
    Precip_monthlymean = np.nanmean(KluaneLake_precip,axis=0)
    
    #plt.figure(figsize=(8,4))
    #plt.title(Kluanelake['Station Name'][0] + ' (' + str(Kluanelake['Latitude (y)'][0]) + ', ' + str(Kluanelake['Longitude (x)'][0]) + ')' , fontsize=14)
    plt.bar(np.arange(1,13),Precip_monthlymean,color='c',zorder=10)
    plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel('Precipitation (mm)',fontsize=14)
    plt.ylim(0,80)
    
plt.figure(figsize=(8,16))
plt.subplot(6,1,1)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100680_1970-1983_P1M.csv')
plt.subplot(6,1,2)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100418_1975-1984_P1M.csv')
plt.subplot(6,1,3)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100182_1966-2007_P1M.csv')
plt.subplot(6,1,4)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100914_1983-1985_P1M.csv')
plt.subplot(6,1,5)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100630_1944-2006_P1M.csv')
plt.subplot(6,1,6)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100631_1985-2007_P1M.csv')
plt.tight_layout()


# =============================================================================
# For Burwash A and Haines Junction YTG, compare monthly precip with downscaled NARR
# Concatenate NARR into monthly precip values
# Get monthly precip values for Burwash and HJ (incl. NaNs)
# Calculate Cobs/Cds
# =============================================================================
def monthly_precip(station_path):
    Kluanelake = pd.read_csv(station_path)
    year_i, year_f = np.array(Kluanelake['Year'])[0], np.array(Kluanelake['Year'])[-1]
    month_i, month_f = np.array(Kluanelake['Month'])[0], np.array(Kluanelake['Month'])[-1]
    if month_f == 12:
        KluaneLake_year_og = np.array(Kluanelake['Year']).reshape(((year_f-year_i)+1,12))
        KluaneLake_year = np.array(Kluanelake['Year']).reshape(((year_f-year_i)+1,12))[np.where(KluaneLake_year_og>=1979)[0][0]:]
        #KluaneLake_month = np.array(Kluanelake['Month']).reshape(((year_f-year_i)+1,12))[np.where(KluaneLake_year_og>=1979)[0][0]:]
        KluaneLake_precip = np.array(Kluanelake['Total Precip (mm)']).reshape(((year_f-year_i)+1,12))[np.where(KluaneLake_year_og>=1979)[0][0]:]
    else:        
        KluaneLake_year_filled = np.concatenate((np.array(Kluanelake['Year']),np.ones((12-month_f))*np.nan))
        KluaneLake_year_allyears = KluaneLake_year_filled.reshape(((year_f-year_i)+1,12))
        KluaneLake_year = KluaneLake_year_allyears[np.where(KluaneLake_year_allyears>=1979)[0][0]:]
        
        KluaneLake_precip_filled = np.concatenate((np.array(Kluanelake['Total Precip (mm)']),np.ones((12-month_f))*np.nan))
        KluaneLake_precip = KluaneLake_precip_filled.reshape(((year_f-year_i)+1,12))[np.where(KluaneLake_year_allyears>=1979)[0][0]:]
    
    return KluaneLake_precip, KluaneLake_year



BurwashA = pd.read_csv('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100182_1966-2007_P1M.csv')
HainesJnct = pd.read_csv('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100631_1985-2007_P1M.csv')

BurwashA_AWS, BurwashA_year = monthly_precip('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100182_1966-2007_P1M.csv')
HainesJnct_AWS, HainesJnct_year = monthly_precip('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100631_1985-2007_P1M.csv')

# Now get the same monthly precip arrays from downscaled NARR precip

inMB = Dataset('D:/Downscaled_files/Burwash_HJ_Stations/Precipitation_Burwash_HJ_1979.nc','r')
precip_array = inMB.variables['Precipitation'][:]
sys.stdout.flush()
inMB.close()

for year in range(1980,2023):
    print(year)
    inMB = Dataset('D:/Downscaled_files/Burwash_HJ_Stations/Precipitation_Burwash_HJ_' + str(year) + '.nc','r')
    precip_array = np.concatenate((precip_array,inMB.variables['Precipitation'][:]))
    precip_array[np.where(precip_array>0.04)] = 0
    sys.stdout.flush()
    inMB.close()
 
date_range = pd.date_range(start= '1979-01-01 00:00:00',end= '2022-12-31 21:00:00',freq='3H')
plt.figure()
plt.plot(date_range,precip_array[:,0,0])

plt.figure()
plt.plot(date_range,precip_array[:,0,1])
    
def monthly_station_precip(Station)    :
    date_range = pd.date_range(start= '1979-01-01 00:00:00',end= '2022-12-31 21:00:00',freq='3H')
    # Create a DataFrame with the dates and precipitation values
    if Station == 'BurwashA':
        data = {'Date': date_range, 'Precipitation': precip_array[:,0,0]}
    elif Station == 'HainesJunction':
        data = {'Date': date_range, 'Precipitation': precip_array[:,0,1]}
    else:
        print('Invalid station name')
        return 'Invalid station name'
        
    df = pd.DataFrame(data)
    # Convert 'Date' column to datetime format (if it's not already)
    df['Date'] = pd.to_datetime(df['Date'])
    # Set 'Date' column as the index
    df.set_index('Date', inplace=True)
    # Resample the data to get monthly sums
    monthly_precipitation = df.resample('M').sum()
    
    return np.array(monthly_precipitation['Precipitation']).reshape(44,12)

Burwash_monthly_precip = monthly_station_precip('BurwashA')
HJ_monthly_precip = monthly_station_precip('HainesJunction')

# Plot monthly downscaled NARR at each station
plt.figure(figsize=(13,6))
plt.subplot(2,2,1)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100182_1966-2007_P1M.csv')
plt.subplot(2,2,2)
plt.title('Downscaled NARR Precip (BURWASH A)', fontsize=14)
plt.bar(np.arange(1,13),np.nanmean(Burwash_monthly_precip,axis=0)*1000,color='darkblue')
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.ylim(0,100)
plt.subplot(2,2,3)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100631_1985-2007_P1M.csv')
plt.subplot(2,2,4)
plt.title('Downscaled NARR Precip (HAINES JUNCTION YTG)', fontsize=14)
plt.bar(np.arange(1,13),np.nanmean(HJ_monthly_precip,axis=0)*1000,color='darkblue')
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.ylim(0,100)
plt.tight_layout()

# Select subset of years from monthly station data
def subset_monthly_narr_precip(year1,year2,station_precip,fullyears):
    return station_precip[(year1-fullyears[0]):len(fullyears)-(fullyears[-1]-year2)]

narryears = np.arange(1979,2022+1)
Burwash_P_1979_2007 = subset_monthly_narr_precip(1979,2007,Burwash_monthly_precip,narryears)
Burwash_P_1979_2007[np.where(np.isnan(BurwashA_AWS))] = np.nan
HJ_P_1985_2007 = subset_monthly_narr_precip(1985,2007,HJ_monthly_precip,narryears)
HJ_P_1985_2007[np.where(np.isnan(HainesJnct_AWS))] = np.nan

plt.figure(figsize=(12,5))
plt.subplot(2,2,1)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100182_1966-2007_P1M.csv')
plt.text(0.3,72,'a) Burwash A (Measured)',fontsize=14,weight='bold')
plt.grid(alpha=0.5)
plt.subplot(2,2,2)
plt.text(0.3,72,'b) Burwash A (Downscaled NARR)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),np.nanmean(Burwash_P_1979_2007,axis=0)*1000,color='midnightblue',zorder=10)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.ylim(0,80)
plt.grid(alpha=0.5)
plt.subplot(2,2,3)
monthly_mean('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/AlternativeAccumulationBiasCorrection/en_climate_monthly_YT_2100631_1985-2007_P1M.csv')
plt.text(0.3,72,'c) Haines Junction YTG (Measured)',fontsize=14,weight='bold')
plt.grid(alpha=0.5)
plt.subplot(2,2,4)
plt.text(0.3,72,'d) Haines Junction YTG (Downscaled NARR)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),np.nanmean(HJ_P_1985_2007,axis=0)*1000,color='midnightblue',zorder=10)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.ylim(0,80)
plt.grid(alpha=0.5,zorder=0)
plt.tight_layout()
#plt.savefig('Monthly_precipitation_HJ_BW.pdf',bbox_inches='tight')


plt.figure(figsize=(10,5))
date_range_BW = pd.date_range(start= '1979-01-01 00:00:00',end= '2007-12-31 21:00:00',freq='M')
plt.subplot(2,1,1)
plt.text(date_range_BW[0],205,'a) Burwash A',fontsize=14,weight='bold')
plt.plot(date_range_BW,Burwash_P_1979_2007.ravel()*1000,label='Downscaled NARR precipitation',color='midnightblue',linewidth=2)
plt.plot(date_range_BW,BurwashA_AWS.ravel(),label='Station precipitation record',color='c',linewidth=2)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(0,200)
plt.ylabel('Precipitation\n(mm/month)',fontsize=14)
plt.legend(fontsize=14)
plt.grid(alpha=0.5)
plt.subplot(2,1,2)
date_range = pd.date_range(start= '1985-01-01 00:00:00',end= '2007-12-31 21:00:00',freq='M')
plt.text(date_range[0],205,'b) Haines Junction YTG',fontsize=14,weight='bold')
plt.plot(date_range,HJ_P_1985_2007.ravel()*1000,label='Downscaled NARR precipitation',color='midnightblue',linewidth=2)
plt.plot(date_range,HainesJnct_AWS.ravel(),label='Station precipitation record',color='c',linewidth=2)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(0,200)
#plt.xticks(date_range_BW[np.arange(12,348,12*4)],date_range_BW[np.arange(12,348,12*4)],rotation=45)
plt.ylabel('Precipitation\n(mm/month)',fontsize=14)
plt.legend(fontsize=14,loc='upper left')
plt.grid(alpha=0.5)
plt.tight_layout()
#plt.savefig('Precipitation_timeseries_HJ_BW.pdf',bbox_inches='tight')



plt.figure(figsize=(12,6))
# =============================================================================
# plt.subplot(3,3,1)
# plt.title('BURWASH A: 1979-2007',fontsize=14)
# plt.bar(np.arange(1,13),np.nanmean(BurwashA_AWS/(Burwash_P_1979_2007*1000),axis=0),color='cadetblue',zorder=10)
# plt.plot(np.arange(1,13),np.nanmean(BurwashA_AWS/(Burwash_P_1979_2007*1000),axis=0),color='k',zorder=20)
# plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
# plt.yticks(fontsize=14)
# plt.hlines(y=1,xmin=0,xmax=13,color='dimgrey',alpha=1,linestyle='--',zorder=15)
# plt.xlim(0.5,12.5)
# plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
# plt.ylim(0,1.8)
# plt.grid(zorder=-50)
# 
# plt.subplot(3,3,2)
# plt.title('HAINES JUNCTION YTG: 1985-2007',fontsize=14)
# plt.bar(np.arange(1,13),np.nanmean(HainesJnct_AWS/(HJ_P_1985_2007*1000),axis=0),color='cadetblue',zorder=10)
# plt.plot(np.arange(1,13),np.nanmean(HainesJnct_AWS/(HJ_P_1985_2007*1000),axis=0),color='k',zorder=20)
# plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
# plt.yticks(fontsize=14)
# plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
# plt.xlim(0.5,12.5)
# plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
# plt.ylim(0,1.8)
# plt.grid(zorder=-50)
# =============================================================================

Burwash_P_1979_2002 = subset_monthly_narr_precip(1979,2002,Burwash_monthly_precip,narryears)
BurwashA_AWS_1979_2002 = subset_monthly_narr_precip(1979,2002,BurwashA_AWS,np.arange(1979,2007+1))

HJ_P_1985_2002 = subset_monthly_narr_precip(1985,2002,HJ_monthly_precip,narryears)
HainesJnct_AWS_1985_2002 = subset_monthly_narr_precip(1985,2002,HainesJnct_AWS,np.arange(1985,2007+1))

plt.subplot(3,2,1)
plt.text(0.6,1.85,'a) Burwash A (Before S.B.)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),np.nanmean(BurwashA_AWS_1979_2002/(Burwash_P_1979_2002*1000),axis=0),color='cadetblue',zorder=10)
plt.plot(np.arange(1,13),np.nanmean(BurwashA_AWS_1979_2002/(Burwash_P_1979_2002*1000),axis=0),color='k',zorder=20)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
plt.xlim(0.5,12.5)
plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
plt.ylim(0,1.8)
plt.grid(zorder=-50)

plt.subplot(3,2,3)
plt.text(0.6,1.85,'b) Haines Junction YTG (Before S.B.)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),np.nanmean(HainesJnct_AWS_1985_2002/(HJ_P_1985_2002*1000),axis=0),color='cadetblue',zorder=10)
plt.plot(np.arange(1,13),np.nanmean(HainesJnct_AWS_1985_2002/(HJ_P_1985_2002*1000),axis=0),color='k',zorder=20)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
plt.xlim(0.5,12.5)
plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
plt.ylim(0,1.8)
plt.grid(zorder=-50)

Burwash_P_2003_2007 = subset_monthly_narr_precip(2003,2007,Burwash_monthly_precip,narryears)
BurwashA_AWS_2003_2007 = subset_monthly_narr_precip(2003,2007,BurwashA_AWS,np.arange(1979,2007+1))

HJ_P_2003_2007 = subset_monthly_narr_precip(2003,2007,HJ_monthly_precip,narryears)
HainesJnct_AWS_2003_2007 = subset_monthly_narr_precip(2003,2007,HainesJnct_AWS,np.arange(1985,2007+1))

plt.subplot(3,2,2)
plt.text(0.6,1.85,'c) Burwash A (After S.B.)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),np.nanmean(BurwashA_AWS_2003_2007/(Burwash_P_2003_2007*1000),axis=0),color='cadetblue',zorder=10)
plt.plot(np.arange(1,13),np.nanmean(BurwashA_AWS_2003_2007/(Burwash_P_2003_2007*1000),axis=0),color='k',zorder=20)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
plt.xlim(0.5,12.5)
plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
plt.ylim(0,1.8)
plt.grid(zorder=-50)

plt.subplot(3,2,4)
plt.text(0.6,1.85,'d) Haines Junction YTG (After S.B.)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),np.nanmean(HainesJnct_AWS_2003_2007/(HJ_P_2003_2007*1000),axis=0),color='cadetblue',zorder=10)
plt.plot(np.arange(1,13),np.nanmean(HainesJnct_AWS_2003_2007/(HJ_P_2003_2007*1000),axis=0),color='k',zorder=20)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
plt.xlim(0.5,12.5)
plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
plt.ylim(0,1.8)
plt.grid(zorder=-50)

def get_distance_from_KW(AWSx,AWSy):
    KWx = np.mean(Xgrid)
    KWy = np.mean(Ygrid)
    
    x_dist = AWSx-KWx
    y_dist = AWSy-KWy
    
    distance = np.sqrt((x_dist**2) + (y_dist**2))
    return distance

Stationx = np.loadtxt('D:\Downscaled_files\Burwash_HJ_Stations\Burwash_HJ_Xgrid.txt')
Stationy = np.loadtxt('D:\Downscaled_files\Burwash_HJ_Stations\Burwash_HJ_Ygrid.txt')

BW_distance = get_distance_from_KW(Stationx[0],Stationy[0])
HJ_distance = get_distance_from_KW(Stationx[1],Stationy[1])

BW_weight = 1/BW_distance
HJ_weight = 1/HJ_distance

Weighted_avg_alldata = ((BW_weight*np.nanmean(BurwashA_AWS/(Burwash_P_1979_2007*1000),axis=0)) + (HJ_weight*np.nanmean(HainesJnct_AWS/(HJ_P_1985_2007*1000),axis=0))) / (BW_weight + HJ_weight)
Weighted_avg_pre2003 =  ((BW_weight*np.nanmean(BurwashA_AWS_1979_2002/(Burwash_P_1979_2002*1000),axis=0)) + (HJ_weight*np.nanmean(HainesJnct_AWS_1985_2002/(HJ_P_1985_2002*1000),axis=0))) / (BW_weight + HJ_weight)
Weighted_avg_post2003 =  ((BW_weight*np.nanmean(BurwashA_AWS_2003_2007/(Burwash_P_2003_2007*1000),axis=0)) + (HJ_weight*np.nanmean(HainesJnct_AWS_2003_2007/(HJ_P_2003_2007*1000),axis=0))) / (BW_weight + HJ_weight)

# =============================================================================
# plt.subplot(3,3,3)
# plt.title('Distance Weighted Average (all data)',fontsize=14)
# plt.bar(np.arange(1,13),Weighted_avg_alldata,color='cadetblue',zorder=10)
# plt.plot(np.arange(1,13),Weighted_avg_alldata,color='k',zorder=20)
# plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
# plt.yticks(fontsize=14)
# plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
# plt.xlim(0.5,12.5)
# plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
# plt.ylim(0,1.8)
# plt.grid(zorder=-50)
# =============================================================================

plt.subplot(3,2,5)
plt.text(0.4,1.85,'e) Distance Weighted Average (Before S.B.)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),Weighted_avg_pre2003,color='cadetblue',zorder=10)
plt.plot(np.arange(1,13),Weighted_avg_pre2003,color='k',zorder=20)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
plt.xlim(0.5,12.5)
plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
plt.ylim(0,1.8)
plt.grid(zorder=-50)

plt.subplot(3,2,6)
plt.text(0.4,1.85,'f) Distance Weighted Average (After S.B.)',fontsize=14,weight='bold')
plt.bar(np.arange(1,13),Weighted_avg_post2003,color='cadetblue',zorder=10)
plt.plot(np.arange(1,13),Weighted_avg_post2003,color='k',zorder=20)
plt.xticks(np.arange(1,13),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=14)
plt.yticks(fontsize=14)
plt.hlines(y=1,xmin=0,xmax=13,color='grey',alpha=1,linestyle='--',zorder=15)
plt.xlim(0.5,12.5)
plt.ylabel('P$_{AWS}$/P$_{ds}$',fontsize=14)
plt.ylim(0,1.8)
plt.grid(zorder=-50)

plt.tight_layout()
#plt.savefig('Paws_Pds_bars.pdf',bbox_inches='tight')


# =============================================================================
# Bias correction:
# =============================================================================
# Create array with shape (2920,239,320) which has the distributed correction
# factors for each year (different pre/post)

# Pre-2003:
BW_Ct = list(np.nanmean(BurwashA_AWS_1979_2002/(Burwash_P_1979_2002*1000),axis=0))
BW_Ct.append(np.nanmean(BurwashA_AWS_1979_2002/(Burwash_P_1979_2002*1000),axis=0)[0])
pre2003_Ct_BW = interp1d(np.asarray([0.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,366.]),np.asarray([BW_Ct]),kind = 'linear') 

HJ_Ct = list(np.nanmean(HainesJnct_AWS_1985_2002/(HJ_P_1985_2002*1000),axis=0))
HJ_Ct.append(np.nanmean(HainesJnct_AWS_1985_2002/(HJ_P_1985_2002*1000),axis=0)[0])
pre2003_Ct_HJ = interp1d(np.asarray([0.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,366.]),np.asarray([HJ_Ct]),kind = 'linear') 

BW_Ct = list(np.nanmean(BurwashA_AWS_2003_2007/(Burwash_P_2003_2007*1000),axis=0))
BW_Ct.append(np.nanmean(BurwashA_AWS_2003_2007/(Burwash_P_2003_2007*1000),axis=0)[0])
post2003_Ct_BW = interp1d(np.asarray([0.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,366.]),np.asarray([BW_Ct]),kind = 'linear') 

HJ_Ct = list(np.nanmean(HainesJnct_AWS_2003_2007/(HJ_P_2003_2007*1000),axis=0))
HJ_Ct.append(np.nanmean(HainesJnct_AWS_2003_2007/(HJ_P_2003_2007*1000),axis=0)[0])
post2003_Ct_HJ = interp1d(np.asarray([0.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,366.]),np.asarray([HJ_Ct]),kind = 'linear') 

BW_x, HJ_x = np.loadtxt('D:\Downscaled_files\Burwash_HJ_Stations\Burwash_HJ_Xgrid.txt')
BW_y, HJ_y = np.loadtxt('D:\Downscaled_files\Burwash_HJ_Stations\Burwash_HJ_Ygrid.txt')
# get distance weighted average correction factor
def distance_weighted_average_Ct(i,j):
    KWx = Xgrid[i][j]
    KWy = Ygrid[i][j]
    
    BW_x_dist = BW_x-KWx
    BW_y_dist = BW_y-KWy
    
    BW_distance = np.sqrt((BW_x_dist**2) + (BW_y_dist**2))
    
    HJ_x_dist = HJ_x-KWx
    HJ_y_dist = HJ_y-KWy
    
    HJ_distance = np.sqrt((HJ_x_dist**2) + (HJ_y_dist**2))

    BW_weight = 1/BW_distance
    HJ_weight = 1/HJ_distance
    
    Weighted_corrfact = ((BW_weight*BW_corrfact) + (HJ_weight*HJ_corrfact)) / (BW_weight + HJ_weight)

    return Weighted_corrfact


pre2003_Ct = np.empty((2928,230,329))
for t in range(0,len(pre2003_Ct)):
    print(t/8)
    pre2003_Ct[t][np.where(np.isnan(Sfc))] = np.nan
    for i in range(0,len(Xgrid)):
        for j in range(0,len(Xgrid[1])):
            if np.isnan(pre2003_Ct[t,i,j]):
                #print('nan')
                pass
            else:
                # Get distance weighted average correction factor
                BW_corrfact = pre2003_Ct_BW(t/8)
                HJ_corrfact = pre2003_Ct_HJ(t/8)
                #print(BW_corrfact,HJ_corrfact)
                
                pre2003_Ct[t,i,j] = distance_weighted_average_Ct(i,j)
    

post2003_Ct = np.empty((2928,230,329))
for t in range(0,len(post2003_Ct)):
    print(t/8)
    post2003_Ct[t][np.where(np.isnan(Sfc))] = np.nan
    for i in range(0,len(Xgrid)):
        for j in range(0,len(Xgrid[1])):
            if np.isnan(post2003_Ct[t,i,j]):
                #print('nan')
                pass
            else:
                # Get distance weighted average correction factor
                BW_corrfact = post2003_Ct_BW(t/8)
                HJ_corrfact = post2003_Ct_HJ(t/8)
                
                post2003_Ct[t,i,j] = distance_weighted_average_Ct(i,j)
    

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,pre2003_Ct[0],cmap='YlGnBu')
plt.axis('equal')
plt.title('Pre-2003 Correction: 01-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,post2003_Ct[0],cmap='YlGnBu')
plt.axis('equal')
plt.title('Post-2003 Correction: 01-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,pre2003_Ct[182*8],cmap='YlGnBu')
plt.axis('equal')
plt.title('Pre-2003 Correction: 07-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,post2003_Ct[182*8],cmap='YlGnBu')
plt.axis('equal')
plt.title('Post-2003 Correction: 07-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,pre2003_Ct[91*8],cmap='YlGnBu')
plt.axis('equal')
plt.title('Pre-2003 Correction: 04-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,post2003_Ct[91*8],cmap='YlGnBu')
plt.axis('equal')
plt.title('Post-2003 Correction: 04-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,pre2003_Ct[274*8],cmap='YlGnBu')
plt.axis('equal')
plt.title('Pre-2003 Correction: 10-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,post2003_Ct[274*8],cmap='YlGnBu')
plt.axis('equal')
plt.title('Post-2003 Correction: 10-01 00:00:00',fontsize=14)
legend = plt.colorbar()
legend.ax.set_ylabel('Distance weighted average C(t)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
plt.tight_layout()

KRH_tributaries[np.isfinite(KRH_tributaries)]=0
KRH_tributaries[np.isnan(KRH_tributaries)]=1

plotletter = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)']
months = ['January','February','March','April','May','June','July','August','September','October','November','December']
timestamps = np.arange(0,365*8,31*8)

fig, axes = plt.subplots(nrows=3, ncols=4,figsize=(15,10))
#fig = plt.figure(figsize=(10,10))
i = 0
for ax in axes.flat:
    #im = ax.imshow(np.random.random((10,10)), vmin=0, vmax=1)
    im = ax.contourf(Xgrid,Ygrid,pre2003_Ct[timestamps[i]],cmap='YlGnBu',levels=np.linspace(0.17,1.38,25))
    #ax.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    #ax.contour(Xgrid,Ygrid,KRH_tributaries,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    ax.axis('equal')
    #ax.set_title(plotletter[i] + ' ' + str(binsizes[i])+ '\n C$_{mean}$ = ' + str(np.round(np.nanmean(biascorrected_snow[i]),2)),fontsize=16)
    ax.set_title(plotletter[i] + ' ' + str(months[i]),fontsize=14)
    #ax.text(Xgrid[-1,55],Ygrid[-1,55]+3000,'C$_{mean}$ = ' + str(np.round(np.nanmean(biascorrected_snow[i]),2)) + ' m w.e. a$^{-1}$',fontsize=14)
    ax.axis('off')
    print(i)
    i+=1

fig.subplots_adjust(right=0.8)
#fig.suptitle('This is a somewhat long figure title', fontsize=16)

cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.02])
legend = fig.colorbar(im, cax=cbar_ax,orientation="horizontal")
legend.ax.set_xlabel('P$_{AWS}$/P$_{ds}$', rotation=0,fontsize=16,labelpad=15)
legend.ax.tick_params(labelsize=16)
legend.set_ticks(np.arange(0.2,1.4,0.1))
legend.set_ticklabels(np.round(np.arange(0.2,1.4,0.1),2))
#plt.savefig('Monthly_Paws_Pds_Pre2003.pdf',bbox_inches='tight')
plt.show()

fig, axes = plt.subplots(nrows=3, ncols=4,figsize=(15,10))
#fig = plt.figure(figsize=(10,10))
i = 0
for ax in axes.flat:
    #im = ax.imshow(np.random.random((10,10)), vmin=0, vmax=1)
    im = ax.contourf(Xgrid,Ygrid,post2003_Ct[timestamps[i]],cmap='YlGnBu',levels=np.linspace(0.17,1.38,25))
    #ax.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    #ax.contour(Xgrid,Ygrid,KRH_tributaries,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    ax.axis('equal')
    #ax.set_title(plotletter[i] + ' ' + str(binsizes[i])+ '\n C$_{mean}$ = ' + str(np.round(np.nanmean(biascorrected_snow[i]),2)),fontsize=16)
    ax.set_title(plotletter[i] + ' ' + str(months[i]),fontsize=14)
    #ax.text(Xgrid[-1,55],Ygrid[-1,55]+3000,'C$_{mean}$ = ' + str(np.round(np.nanmean(biascorrected_snow[i]),2)) + ' m w.e. a$^{-1}$',fontsize=14)
    ax.axis('off')
    print(i)
    i+=1

fig.subplots_adjust(right=0.8)
#fig.suptitle('This is a somewhat long figure title', fontsize=16)

cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.02])
legend = fig.colorbar(im, cax=cbar_ax,orientation="horizontal")
legend.ax.set_xlabel('P$_{AWS}$/P$_{ds}$', rotation=0,fontsize=16,labelpad=15)
legend.ax.tick_params(labelsize=16)
legend.set_ticks(np.arange(0.2,1.4,0.1))
legend.set_ticklabels(np.round(np.arange(0.2,1.4,0.1),2))
#plt.savefig('Monthly_Paws_Pds_Post2003.pdf',bbox_inches='tight')
plt.show()


# =============================================================================
# Bias-correct the downscaled precipitation files
# =============================================================================

Glacier_ID = 'KRH'
OUTPUT_PATH = 'D:\BiasCorrected_files\KRH\AWS_based_correction'
# For year in years
for year in range(1979,2023):
    print(year)
    # Load downscaled (uncorrected) precip file
    filename = 'D:/Downscaled_files/KRH/Precipitation_KRH_' + str(year) + '.nc'
    inP = Dataset(filename,'r')
    P_array = inP.variables['Precipitation'][:]
    sys.stdout.flush()
    inP.close()
    
    # Leap year? Get index for the correction factor array
    if len(P_array) == 2920:
        ind = -8 # array will have shape 2920
    elif len(P_array) == 2928:
        ind = None # array will have shape 2928
    
    if year < 2003:
        Corrected_P = P_array * pre2003_Ct[:ind]
    else:
        Corrected_P = P_array * post2003_Ct[:ind]
        
    # save corrected precip as netcdf file
    Corrected_P_file = os.path.join(OUTPUT_PATH,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc')
    #save_to_netcdf(Corrected_P, 'Precipitation', Corrected_P_file, year, Xgrid, Ygrid, dt=3, precision=32)  

# =============================================================================
# Plot mean annual bias corrected accumulation (before and after 2003)
# =============================================================================
    
BC_INPUT_PATH = 'D:\BiasCorrected_files\KRH\AWS_based_correction'
BCSnow = np.empty((44,230,329))
for year in range(1979,2022+1):
    print(year)
    dataset = Dataset(os.path.join(BC_INPUT_PATH,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    P = dataset.variables['Precipitation'][:]
    sys.stdout.flush()
    
    meanP = np.sum(P,axis=0)
    BCSnow[year-1979,:,:] = meanP

Avg_BCSnow = np.mean(BCSnow,axis=0)

plt.figure(figsize=(9,5))
plt.contourf(Xgrid,Ygrid,Avg_BCSnow,cmap='BuPu',levels=np.linspace(0.28,0.63,16))
#plt.pcolor(Xgrid,np.flipud(Ygrid),catchment_tribarray,cmap='YlGnBu')
plt.axis('equal')
legend = plt.colorbar()
legend.ax.set_ylabel('Accumulation (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
plt.xlabel('Easting (m)',fontsize=14)
plt.ylabel('Northing (m)',fontsize=14)
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
legend.ax.tick_params(labelsize=14)
plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
legend.set_ticks(np.arange(0.3,0.61,0.05))
legend.set_ticklabels(np.round(np.arange(0.3,0.61,0.05),2))
plt.tight_layout()
#plt.savefig('AWSBiasCorrectedAccumulation_1979-2022.pdf',bbox_inches='tight')

# =============================================================================
# Calculate mean annual accumulation in each zone
# =============================================================================
print('Trunk',np.round(np.mean(Avg_BCSnow[np.where(KRH_tributaries==5)]),2))
print('NA',np.round(np.mean(Avg_BCSnow[np.where(KRH_tributaries==4)]),2))
print('CA',np.round(np.mean(Avg_BCSnow[np.where(KRH_tributaries==3)]),2))
print('SW',np.round(np.mean(Avg_BCSnow[np.where(KRH_tributaries==2)]),2))
print('SA',np.round(np.mean(Avg_BCSnow[np.where(KRH_tributaries==1)]),2))

nonKWicemask = np.array(Sfc)
nonKWicemask[np.isfinite(KRH_tributaries)] = 1

print('Non KW ice',np.round(np.mean(Avg_BCSnow[np.where(nonKWicemask==0)]),2))
print('Off-glacier',np.round(np.mean(Avg_BCSnow[np.where(Sfc==1)]),2))
print('All KW',np.round(np.mean(Avg_BCSnow[np.isfinite(KRH_tributaries)]),2))
print('KRH',np.round(np.nanmean(Avg_BCSnow),2))
