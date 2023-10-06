# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:10:01 2023

Functions for calculating and plotting the outputs
from the mass-balance model

@author: katierobinson
"""

import numpy as np
import pandas as pd
import os
import sys
from netCDF4 import Dataset

def average_sim_results():
    '''
    Get the average result from a set of simulations
    Output sim9999
    '''

def load_hydrologic_year(sim,year,varname,var,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,NARRvar=False):
    '''
    Returns a given variable from Oct 1 -- Sept 30 (hydrological year)
    '''
    
    dates_yr1 = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq=str(3)+'H')
    dates_yr2 = pd.date_range(start= str(year+1) + '-01-01 00:00:00',end= str(year+1) + '-12-31 21:00:00',freq=str(3)+'H')
        
    # Concatenate variable from Oct 1 of current year to Sept 30 of following year    
    if NARRvar == True:
        inMB1 = Dataset(os.path.join(NARR_INPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
        inMB2 = Dataset(os.path.join(NARR_INPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year+1) + '.nc'),'r')
    else:
        inMB1 = Dataset(os.path.join(MODEL_OUTPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
        inMB2 = Dataset(os.path.join(MODEL_OUTPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year+1) + '_' + str(sim) + '.nc'),'r')
    
    mb_array_year1 = inMB1.variables[var][:]
    sys.stdout.flush()
    inMB1.close()
    
    mb_array_year2 = inMB2.variables[var][:]
    sys.stdout.flush()
    inMB2.close()
    
    mb_yr1 = mb_array_year1[np.where(dates_yr1 == pd.Timestamp(str(year)+'-10-01T00'))[0][0]:]
    mb_yr2 = mb_array_year2[:np.where(dates_yr2 == pd.Timestamp(str(year+1)+'-10-01T00'))[0][0]]
    mb_hydroyear = np.concatenate((mb_yr1,mb_yr2),axis=0)
    
    return mb_hydroyear

def calculate_mb_components_timeseries(sim,years,R2S,Glacier_grid,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,x,y,PointScale=True,KWaverage=False,Catchmentaverage=False):
    '''
    sim = integer corresponding to the params used to run the model
    years = array of years to be calculated
    R2S = rain to snow threshold
    Glacier_grid = array where KW gridcells = int, otherwise NaN
    x,y = coords for plotting variables at a single gridcell (PointScale = True)
    KWaverage: True if plotting the kaskawulsh wide average, otherwise false
    Catchmentaverage: True if plotting the catchment-wide average, otherwise false
    
    returns:
    daily runoff components:
        net snowmelt (total snow melt minus refreezing)
        superimposed ice melt
        glacier ice melt
        rain runoff (total rain minus refrozen rain)
    daily mass balance components:
        accumulation
        rain that refreezes
    '''
    massbal_l = []
    totalsnowmelt_l = []
    refrozen_melt_l = []
    snowmelt_runoff_l = []
    superimposed_icemelt_l = []
    glacier_icemelt_l = []
    rain_l = []
    rain_runoff_l = []
    rain_refreezing_l = []
    accumulation_l = []
    
    for year in years:
        print('Calculating runoff & mb variables for', year)
        # Get mass balance variable:
        massbal = load_hydrologic_year(sim,year,'Netbalance','Net balance',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        
        # Calculate the snow melt that runs off:
        snowmelt = load_hydrologic_year(sim,year,'Snowmelt','Snow melt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        refreezing = load_hydrologic_year(sim,year,'Refrozenmelt','Refreezing melt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        
        netsnowmelt = np.subtract(snowmelt,refreezing)
        
        # Calculate the ice melt (glacier ice and superimposed ice):
        icemelt = load_hydrologic_year(sim,year,'Icemelt','Ice melt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        total_superimposed_ice = load_hydrologic_year(sim,year,'Superimposedice','Superimposed ice',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        
        glacier_icemelt = np.subtract(icemelt,total_superimposed_ice)
        glacier_icemelt[np.where(glacier_icemelt < 0)] = 0
        
        superimposed_icemelt = np.subtract(icemelt,glacier_icemelt)
        
        # Calculate the runoff from rain and rain that refreezes:
        precip = load_hydrologic_year(sim,year,'Precipitation','Precipitation',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,NARRvar=True)
        temp = load_hydrologic_year(sim,year,'Temperature','Temperature',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,NARRvar=True)
        rain_refreezing = load_hydrologic_year(sim,year,'Refrozenrain','Refreezing rain',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        rain = np.array(precip)
        rain[np.where(temp <= R2S)] = 0
        
        rain_runoff = np.subtract(rain,rain_refreezing)
        
        # Calculate the accumulation:
        accumulation = np.array(precip)
        accumulation[np.where(temp > R2S)] = 0
        
        # Get timeseries from each component at a point location:
        if PointScale == True:
            daily_massbal = np.sum(massbal[:,x,y].reshape(-1,8),axis=1)
            daily_snowmelt = np.sum(snowmelt[:,x,y].reshape(-1,8),axis=1)
            daily_refrozenmelt = np.sum(refreezing[:,x,y].reshape(-1,8),axis=1)
            daily_netsnowmelt = np.sum(netsnowmelt[:,x,y].reshape(-1,8),axis=1)
            daily_glaciericemelt = np.sum(glacier_icemelt[:,x,y].reshape(-1,8),axis=1)
            daily_superimposedicemelt = np.sum(superimposed_icemelt[:,x,y].reshape(-1,8),axis=1)
            daily_rainrunoff = np.sum(rain_runoff[:,x,y].reshape(-1,8),axis=1)
            daily_rain = np.sum(rain[:,x,y].reshape(-1,8),axis=1)
            daily_rainrefreezing = np.sum(rain_refreezing[:,x,y].reshape(-1,8),axis=1)
            daily_accumulation = np.sum(accumulation[:,x,y].reshape(-1,8),axis=1)
        
        # Get timeseries from each component for the catchment wide average
        elif Catchmentaverage == True:
            mean_massbal = np.nanmean(np.nanmean(massbal,axis=1),axis=1)
            mean_snowmelt =  np.nanmean(np.nanmean(snowmelt,axis=1),axis=1)
            mean_refrozenmelt = np.nanmean(np.nanmean(refreezing,axis=1),axis=1)
            mean_netsnowmelt = np.nanmean(np.nanmean(netsnowmelt,axis=1),axis=1)
            mean_glaciericemelt = np.nanmean(np.nanmean(glacier_icemelt,axis=1),axis=1)
            mean_superimposedicemelt = np.nanmean(np.nanmean(superimposed_icemelt,axis=1),axis=1)
            mean_rainrunoff = np.nanmean(np.nanmean(rain_runoff,axis=1),axis=1)
            mean_rain = np.nanmean(np.nanmean(rain,axis=1),axis=1)
            mean_rainrefreezing = np.nanmean(np.nanmean(rain_refreezing,axis=1),axis=1)
            mean_accumulation = np.nanmean(np.nanmean(accumulation,axis=1),axis=1)
            
            daily_massbal = np.sum(mean_massbal.reshape(-1,8),axis=1)
            daily_snowmelt = np.sum(mean_snowmelt[:,x,y].reshape(-1,8),axis=1)
            daily_refrozenmelt = np.sum(mean_refrozenmelt[:,x,y].reshape(-1,8),axis=1)
            daily_netsnowmelt = np.sum(mean_netsnowmelt.reshape(-1,8),axis=1)
            daily_glaciericemelt = np.sum(mean_glaciericemelt.reshape(-1,8),axis=1)
            daily_superimposedicemelt = np.sum(mean_superimposedicemelt.reshape(-1,8),axis=1)
            daily_rainrunoff = np.sum(mean_rainrunoff.reshape(-1,8),axis=1)
            daily_rainrefreezing = np.sum(mean_rainrefreezing.reshape(-1,8),axis=1)
            daily_rain = np.sum(mean_rain.reshape(-1,8),axis=1)
            daily_accumulation = np.sum(mean_accumulation.reshape(-1,8),axis=1)
        
        # Get timeseries from each component for the glacier wide average
        elif KWaverage == True:
            for i in range(0,len(snowmelt)):
                massbal[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                snowmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                refreezing[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                netsnowmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                glacier_icemelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                superimposed_icemelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                rain_runoff[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                rain_refreezing[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                rain[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                accumulation[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
    
            mean_massbal = np.nanmean(np.nanmean(massbal,axis=1),axis=1)
            mean_snowmelt =  np.nanmean(np.nanmean(snowmelt,axis=1),axis=1)
            mean_refrozenmelt = np.nanmean(np.nanmean(refreezing,axis=1),axis=1)
            mean_netsnowmelt = np.nanmean(np.nanmean(netsnowmelt,axis=1),axis=1)
            mean_glaciericemelt = np.nanmean(np.nanmean(glacier_icemelt,axis=1),axis=1)
            mean_superimposedicemelt = np.nanmean(np.nanmean(superimposed_icemelt,axis=1),axis=1)
            mean_rainrunoff = np.nanmean(np.nanmean(rain_runoff,axis=1),axis=1)
            mean_rain = np.nanmean(np.nanmean(rain,axis=1),axis=1)
            mean_rainrefreezing = np.nanmean(np.nanmean(rain_refreezing,axis=1),axis=1)
            mean_accumulation = np.nanmean(np.nanmean(accumulation,axis=1),axis=1)
            
            daily_massbal = np.sum(mean_massbal.reshape(-1,8),axis=1)
            daily_snowmelt = np.sum(mean_snowmelt.reshape(-1,8),axis=1)
            daily_refrozenmelt= np.sum(mean_refrozenmelt.reshape(-1,8),axis=1)
            daily_netsnowmelt = np.sum(mean_netsnowmelt.reshape(-1,8),axis=1)
            daily_glaciericemelt = np.sum(mean_glaciericemelt.reshape(-1,8),axis=1)
            daily_superimposedicemelt = np.sum(mean_superimposedicemelt.reshape(-1,8),axis=1)
            daily_rainrunoff = np.sum(mean_rainrunoff.reshape(-1,8),axis=1)
            daily_rain = np.sum(mean_rain.reshape(-1,8),axis=1)
            daily_rainrefreezing = np.sum(mean_rainrefreezing.reshape(-1,8),axis=1)
            daily_accumulation = np.sum(mean_accumulation.reshape(-1,8),axis=1)
    
        massbal_l.append(daily_massbal)
        totalsnowmelt_l.append(daily_snowmelt)
        refrozen_melt_l.append(daily_refrozenmelt)
        snowmelt_runoff_l.append(daily_netsnowmelt)
        superimposed_icemelt_l.append(daily_superimposedicemelt)
        glacier_icemelt_l.append(daily_glaciericemelt)
        rain_runoff_l.append(daily_rainrunoff)
        rain_refreezing_l.append(daily_rainrefreezing)
        rain_l.append(daily_rain)
        accumulation_l.append(daily_accumulation)
        
    return massbal_l, totalsnowmelt_l, refrozen_melt_l, snowmelt_runoff_l, superimposed_icemelt_l, glacier_icemelt_l, rain_runoff_l, rain_refreezing_l, rain_l, accumulation_l


def calculate_mb_components_distributed(sim,years,R2S,Glacier_grid,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID):
    '''
    sim = integer corresponding to the params used to run the model
    years = array of years to be calculated
    R2S = rain to snow threshold
    Glacier_grid = array where KW gridcells = int, otherwise NaN
    x,y = coords for plotting variables at a single gridcell (PointScale = True)
    KWaverage: True if plotting the kaskawulsh wide average, otherwise false
    Catchmentaverage: True if plotting the catchment-wide average, otherwise false
    
    returns:
    daily runoff components:
        net snowmelt (total snow melt minus refreezing)
        superimposed ice melt
        glacier ice melt
        rain runoff (total rain minus refrozen rain)
    daily mass balance components:
        accumulation
        rain that refreezes
    '''
    
    massbal_l = []
    total_snowmelt_l = []
    refrozen_melt_l = []
    snowmelt_runoff_l = []
    superimposed_icemelt_l = []
    glacier_icemelt_l = []
    rain_l = []
    rain_runoff_l = []
    rain_refreezing_l = []
    accumulation_l = []
    
    for year in years:
        print('Calculating runoff & mb variables for', year)
        
        # Get mass balance variable:
        massbal = load_hydrologic_year(sim,year,'Netbalance','Net balance',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        cumulative_massbal = np.sum(massbal,axis=0)
        
        # Calculate the snow melt that runs off:
        snowmelt = load_hydrologic_year(sim,year,'Snowmelt','Snow melt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        total_snowmelt = np.sum(snowmelt,axis=0)
        
        refreezing = load_hydrologic_year(sim,year,'Refrozenmelt','Refreezing melt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        total_refreezing = np.sum(refreezing,axis=0)
        
        netsnowmelt = np.subtract(snowmelt,refreezing)
        total_netsnowmelt = np.sum(netsnowmelt,axis=0)
        
        # Calculate the ice melt (glacier ice and superimposed ice):
        icemelt = load_hydrologic_year(sim,year,'Icemelt','Ice melt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        total_superimposed_ice = load_hydrologic_year(sim,year,'Superimposedice','Superimposed ice',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        
        glacier_icemelt = np.subtract(icemelt,total_superimposed_ice)
        glacier_icemelt[np.where(glacier_icemelt < 0)] = 0
        total_glacier_icemelt = np.sum(glacier_icemelt,axis=0)
        
        superimposed_icemelt = np.subtract(icemelt,glacier_icemelt)
        total_superimposed_icemelt = np.sum(superimposed_icemelt,axis=0)
        
        # Calculate the runoff from rain and rain that refreezes:
        precip = load_hydrologic_year(sim,year,'Precipitation','Precipitation',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,NARRvar=True)
        temp = load_hydrologic_year(sim,year,'Temperature','Temperature',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,NARRvar=True)
        rain_refreezing = load_hydrologic_year(sim,year,'Refrozenrain','Refreezing rain',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        total_rain_refreezing = np.sum(rain_refreezing,axis=0)
        
        rain = np.array(precip)
        rain[np.where(temp <= R2S)] = 0
        total_rain = np.sum(rain,axis=0)
        
        rain_runoff = np.subtract(rain,rain_refreezing)
        total_rain_runoff = np.sum(rain_runoff,axis=0)
        
        # Calculate the accumulation:
        accumulation = np.array(precip)
        accumulation[np.where(temp > R2S)] = 0
        total_accumulation = np.sum(accumulation,axis=0)
    
        # KR_note: left off here
        massbal_l.append(cumulative_massbal)
        total_snowmelt_l.append(total_snowmelt)
        refrozen_melt_l.append(total_refreezing)
        snowmelt_runoff_l.append(total_netsnowmelt)
        superimposed_icemelt_l.append(total_superimposed_icemelt)
        glacier_icemelt_l.append(total_glacier_icemelt)
        rain_runoff_l.append(total_rain_runoff)
        rain_refreezing_l.append(total_rain_refreezing)
        rain_l.append(total_rain)
        accumulation_l.append(total_accumulation)
        
    return massbal_l, total_snowmelt_l, refrozen_melt_l, snowmelt_runoff_l, superimposed_icemelt_l, glacier_icemelt_l, rain_runoff_l, rain_refreezing_l, rain_l, accumulation_l

