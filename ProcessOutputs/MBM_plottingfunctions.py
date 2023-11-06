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
import matplotlib.pyplot as plt


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
            daily_snowmelt = np.sum(mean_snowmelt.reshape(-1,8),axis=1)
            daily_refrozenmelt = np.sum(mean_refrozenmelt.reshape(-1,8),axis=1)
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

def runoff_timeseries_12years(title,year1, years, daily_runoff_upperlim, cumu_runoff_upperlim, gl_icemelt,snowmelt_runoff, rain_runoff, superimp_icemelt):
    a = []
    fig, axs = plt.subplots(nrows=4, ncols=3,figsize=(12,12))
    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            a.append(ax)
    
    for i, ax in enumerate(a):
        year = i+year1
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
        
        ax.set_title(str(year)+'-'+str(year+1))
        ax.set_ylabel('Runoff (m w.e. d$^{-1}$)')
    
        ax.plot(np.arange(0,len(dates)),gl_icemelt[year-years[0]],c='turquoise',label='Glacier ice melt')    
        ax.plot(np.arange(0,len(dates)),snowmelt_runoff[year-years[0]],c='royalblue',label='Snow melt')
        ax.plot(np.arange(0,len(dates)),rain_runoff[year-years[0]],c='deeppink',label='Rain')
        ax.plot(np.arange(0,len(dates)),superimp_icemelt[year-years[0]],c='darkorange',label='Superimposed ice melt')
        ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
        ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45)
        ax.set_ylim(0,daily_runoff_upperlim)
        ax.grid()
        ax.margins(x=0)
    
        total_runoff = gl_icemelt[year-years[0]] + snowmelt_runoff[year-years[0]] + rain_runoff[year-years[0]] + superimp_icemelt[year-years[0]]
    
        ax0 = ax.twinx()
        ax0.plot(np.arange(0,len(dates)),np.cumsum(total_runoff),c='k',label='cumulative runoff')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(gl_icemelt[year-years[0]]),c='turquoise',label='cumulative runoff')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(snowmelt[year-years[0]]),c='royalblue',label='cumulative runoff')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(superimp_icemelt[year-years[0]]),c='orange')
        ax0.set_ylim(0,cumu_runoff_upperlim)
        ax0.set_ylabel('Cumulative Runoff (m w.e.)')
        ax0.margins(x=0)
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),fontsize=14,bbox_to_anchor=(0.98,1.08), ncol=2, borderaxespad=0.19)
    fig.suptitle(title,fontsize=14,y=1.01) 
    fig.tight_layout()
    #fig.savefig(os.path.join(MODEL_OUTPUTS,'KWaverage_runoff_2007-2018.png'),bbox_inches='tight')

def runoff_piecharts_12years(title,year1,years,gl_icemelt, netsnowmelt, rain_runoff, superimp_icemelt):

    piechart_colours = ['turquoise','royalblue','deeppink','darkorange']
    fig = plt.figure(figsize=(10,14))
    for i in range(1,13):
        year = i+(year1-1)
    
        plt.subplot(4,3,i)
        plt.title(str(year)+'-'+str(year+1),fontsize=14)
        
        total_runoff = np.sum(gl_icemelt[year-years[0]] + netsnowmelt[year-years[0]] + rain_runoff[year-years[0]] + superimp_icemelt[year-years[0]])
        gl_ice_percent = round((np.sum(gl_icemelt[year-years[0]])/total_runoff)*100,1)
        snow_percent = round((np.sum(netsnowmelt[year-years[0]])/total_runoff)*100,1)
        SIice_percent = round((np.sum(superimp_icemelt[year-years[0]])/total_runoff)*100,1)
        rain_percent = round((np.sum(rain_runoff[year-years[0]])/total_runoff)*100,1)
        
        plt.pie([gl_ice_percent,snow_percent,rain_percent,SIice_percent],autopct='%.1f%%',colors=piechart_colours, textprops={'fontsize': 14})
        #plt.pie([gl_ice_percent,snow_percent,rain_percent,SIice_percent],labels=[str(gl_ice_percent) + '%',str(snow_percent) + '%',str(rain_percent) + '%',str(SIice_percent) + '%'],colors=piechart_colours, textprops={'fontsize': 12})
    fig.suptitle(title,fontsize=14,y=1.01) 
    fig.legend(['Glacier ice melt','Snow melt','Rain','Superimposed ice melt'],fontsize=14,bbox_to_anchor=(1,1.03),ncol=2)
    plt.tight_layout()
    
def runoff_timeseries_average(title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt):
    
    ice_sum = np.zeros((365))
    snow_sum = np.zeros((365))
    rain_sum = np.zeros((365))
    SI_sum = np.zeros((365))
    
    for year in avg_years:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt[i][:365])
        snow_sum += np.array(snowmelt_runoff[i][:365])
        rain_sum += np.array(rain_runoff[i][:365])
        SI_sum += np.array(superimp_icemelt[i][:365])
        
    ice_mean = np.array(ice_sum/len(avg_years))
    snow_mean = np.array(snow_sum/len(avg_years))
    rain_mean = np.array(rain_sum/len(avg_years))
    SI_mean = np.array(SI_sum/len(avg_years))
    
    # Plot the average
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')       
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    ax.set_ylabel('Runoff (m w.e. d$^{-1}$)',fontsize=14)
    
    ax.plot(np.arange(0,len(dates)),ice_mean,c='turquoise',label='Glacier ice melt')    
    ax.plot(np.arange(0,len(dates)),snow_mean,c='royalblue',label='Snow melt')
    ax.plot(np.arange(0,len(dates)),rain_mean,c='deeppink',label='Rain')
    ax.plot(np.arange(0,len(dates)),SI_mean,c='darkorange',label='Superimposed ice melt')
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(0,daily_runoff_upperlim)
    ax.grid()
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    
    total_runoff = ice_mean + snow_mean + rain_mean + SI_mean
    
    ax0 = ax.twinx()
    ax0.plot(np.arange(0,len(dates)),np.cumsum(total_runoff),c='k',label='Cumulative runoff',linewidth=3)
    ax0.plot(np.arange(0,len(dates)),np.cumsum(ice_mean),c='turquoise',label='Glacier ice melt',linewidth=3)
    ax0.plot(np.arange(0,len(dates)),np.cumsum(snow_mean),c='royalblue',label='Snow melt',linewidth=3)
    ax0.plot(np.arange(0,len(dates)),np.cumsum(rain_mean),c='deeppink',label='Rain',linewidth=3)
    ax0.plot(np.arange(0,len(dates)),np.cumsum(SI_mean),c='darkorange',label='Superimposed ice melt',linewidth=3)
    ax0.set_ylim(0,cumu_runoff_upperlim)
    ax0.set_ylabel('Cumulative Runoff (m w.e.)',fontsize=14)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.margins(x=0)
    
    handles, labels = ax0.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),loc='upper left',fontsize=14, ncol=1, borderaxespad=0.19)
    fig.tight_layout()
    
    # Pie chart
    piechart_colours = ['turquoise','royalblue','deeppink','darkorange']
    fig = plt.figure(figsize=(5,5))
    plt.title(str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    
    total_runoff = np.sum(ice_mean + snow_mean + SI_mean + rain_mean)
    gl_ice_percent = round((np.sum(ice_mean)/total_runoff)*100,1)
    snow_percent = round((np.sum(snow_mean)/total_runoff)*100,1)
    SIice_percent = round((np.sum(SI_mean)/total_runoff)*100,1)
    rain_percent = round((np.sum(rain_mean)/total_runoff)*100,1)
    
    plt.pie([gl_ice_percent,snow_percent,rain_percent,SIice_percent],autopct='%.1f%%',colors=piechart_colours, textprops={'fontsize': 15})
    #fig.suptitle(title,fontsize=14,y=1.01) 
    #plt.legend(['Glacier ice melt','Snow melt','Rain','Superimposed ice melt'],fontsize=14,ncol=1,bbox_to_anchor=(1.3,0.8))
    plt.tight_layout()
    
def distributed_average_mass_balance(avg_years,all_years,dist_massbal,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        print(year)
        MB_ARRAY[year-avg_years[0],:,:] = dist_massbal[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.mean((np.nanmean(MB_ARRAY,axis=0))[np.isfinite(KRH_tributaries)]),2)
    print(np.nanmin(np.nanmean(MB_ARRAY,axis=0)),np.nanmax(np.nanmean(MB_ARRAY,axis=0)))  
    
    plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),cmap='massbal',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-7,3))
    legend.ax.set_ylabel('Mass Balance (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    #plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]) + '\nmass balance\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    plt.tight_layout()

def massbalance_timeseries_12years(title,year1, years, daily_mb_abs_lim, cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt):
    a = []
    fig, axs = plt.subplots(nrows=4, ncols=3,figsize=(12,12))
    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            a.append(ax)
    
    for i, ax in enumerate(a):
        year = i+year1
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
        
        ax.set_title(str(year)+'-'+str(year+1))
        ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)')
    
        total_accumulation = accumulation[year-years[0]] + refrozen_rain[year-years[0]]
        total_ablation = netsnowmelt[year-years[0]] + superimp_icemelt[year-years[0]] + gl_icemelt[year-years[0]]
    
        ax.plot(np.arange(0,len(dates)),total_accumulation,c='mediumblue',label='Accumulation')    
        ax.plot(np.arange(0,len(dates)),-total_ablation,c='red',label='Ablation')    
        ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
        ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45)
        ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
        ax.grid()
        ax.tick_params(axis='y',labelsize=14)
        ax.margins(x=0)
        
        massbal = (accumulation[year-years[0]] + refrozen_rain[year-years[0]]) - (netsnowmelt[year-years[0]] + superimp_icemelt[year-years[0]] + gl_icemelt[year-years[0]])
    
        d_tr = [np.where(np.cumsum(massbal)[50:] <=0)][0][0]
        if len(d_tr) == 0:
            transition_date = 'N/A'
        else:
            transition_date = str(dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]])[:10]
        ax0 = ax.twinx()
        ax0.plot(np.arange(0,len(dates)),np.cumsum(massbal),c='k',label='$\dot{B}$ = 0 on ' + str(transition_date)[:10])
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(gl_icemelt[year-years[0]]),c='turquoise',label='cumulative runoff')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(snowmelt[year-years[0]]),c='royalblue',label='cumulative runoff')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(superimp_icemelt[year-years[0]]),c='orange')
        ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
        ax0.set_ylabel('Cumulative Mass Balance (m w.e.)')
        ax0.margins(x=0)
        ax0.tick_params(axis='y',labelsize=14)
        ax0.legend(loc='upper left')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),fontsize=14,bbox_to_anchor=(0.98,1.04), ncol=2, borderaxespad=0.19)
    fig.suptitle(title,fontsize=14,y=1.01) 
    fig.tight_layout()
    
    
def massbalance_timeseries_average(title,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt):
    
    snowfall_sum = np.zeros((365))
    refrain_sum = np.zeros((365))
    snowmelt_sum = np.zeros((365))
    SImelt_sum = np.zeros((365))
    gl_melt_sum = np.zeros((365))
    
    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation[i][:365])
        refrain_sum += np.array(refrozen_rain[i][:365])
        snowmelt_sum += np.array(netsnowmelt[i][:365])
        SImelt_sum += np.array(superimp_icemelt[i][:365])
        gl_melt_sum += np.array(gl_icemelt[i][:365]) 
        
    snowfall_mean = np.array(snowfall_sum/len(avg_years))
    refrain_mean = np.array(refrain_sum/len(avg_years))
    snowmelt_mean = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean = np.array(gl_melt_sum/len(avg_years))

    
    # Plot the average
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')       
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=14)
    
    total_accumulation = snowfall_mean + refrain_mean
    total_ablation = snowmelt_mean + SImelt_mean + gl_melt_mean

    ax.plot(np.arange(0,len(dates)),total_accumulation,c='mediumblue',label='Accumulation')    
    ax.plot(np.arange(0,len(dates)),-total_ablation,c='red',label='Ablation')    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax.grid()
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    
    massbal = total_accumulation - total_ablation

    transition_date = dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]]
    ax0 = ax.twinx()
    ax0.plot(np.arange(0,len(dates)),np.cumsum(massbal),c='k',label='Cumulative balance\n$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    #ax0.plot(np.arange(0,len(dates)),np.cumsum(gl_icemelt[year-years[0]]),c='turquoise',label='cumulative runoff')
    #ax0.plot(np.arange(0,len(dates)),np.cumsum(snowmelt[year-years[0]]),c='royalblue',label='cumulative runoff')
    #ax0.plot(np.arange(0,len(dates)),np.cumsum(superimp_icemelt[year-years[0]]),c='orange')
    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass Balance (m w.e.)',fontsize=14)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.legend(fontsize=14,loc='lower left')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
    fig.tight_layout()


def snowrunoff_timeseries_12years(title,year1, years, daily_runoff_lower_lim, daily_runoff_upper_lim, cumu_runoff_lower_lim, cumu_runoff_upper_lim, totalsnowmelt, refreezing):
    a = []
    fig, axs = plt.subplots(nrows=4, ncols=3,figsize=(12,12))
    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            a.append(ax)
    
    for i, ax in enumerate(a):
        year = i+year1
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
        
        ax.set_title(str(year)+'-'+str(year+1))
        ax.set_ylabel('Snow melt/refreezing (m w.e. d$^{-1}$)')
    
        ax.plot(np.arange(0,len(dates)),totalsnowmelt[year-years[0]],c='darkslateblue',label='Total snow melt')    
        ax.plot(np.arange(0,len(dates)),-refreezing[year-years[0]],c='aquamarine',label='Refrozen melt')    
        ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
        ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45)
        ax.set_ylim(daily_runoff_lower_lim,daily_runoff_upper_lim)
        ax.grid()
        ax.margins(x=0)
        
        netsnowmelt = totalsnowmelt[year-years[0]] - refreezing[year-years[0]]
    
        #transition_date = dates[50:][np.where(np.cumsum(netsnowmelt)[50:] <=0)[0][0]]
        ax0 = ax.twinx()
        ax0.plot(np.arange(0,len(dates)),np.cumsum(netsnowmelt),c='k',label='Cumulative runoff from snow melt')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(gl_icemelt[year-years[0]]),c='turquoise',label='cumulative runoff')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(snowmelt[year-years[0]]),c='royalblue',label='cumulative runoff')
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(superimp_icemelt[year-years[0]]),c='orange')
        ax0.set_ylim(cumu_runoff_lower_lim,cumu_runoff_upper_lim)
        ax0.set_ylabel('Cumulative Runoff (m w.e.)')
        ax0.margins(x=0)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),fontsize=14,bbox_to_anchor=(0.98,1.04), ncol=2, borderaxespad=0.19)
    fig.suptitle(title,fontsize=14,y=1.01) 
    fig.tight_layout()
    
def snowrunoff_timeseries_average(title,avg_years,all_years,daily_runoff_lower_lim, daily_runoff_upper_lim, cumu_runoff_lower_lim, cumu_runoff_upper_lim, totalsnowmelt, refreezing):
    
    total_snowmelt_sum = np.zeros((365))
    refreezing_sum = np.zeros((365))
    
    for year in avg_years:
        i = year - all_years[0]
        total_snowmelt_sum += np.array(totalsnowmelt[i][:365])
        refreezing_sum += np.array(refreezing[i][:365])
        
    snowmelt_mean = np.array(total_snowmelt_sum/len(avg_years))
    refreezing_mean = np.array(refreezing_sum/len(avg_years))

    # Plot the average
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')       
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]),fontsize=14)
    ax.set_ylabel('Snow melt/refreezing (m w.e. d$^{-1}$)',fontsize=14)
    
    ax.plot(np.arange(0,len(dates)),snowmelt_mean,c='darkslateblue',label='Total snow melt')    
    ax.plot(np.arange(0,len(dates)),-refreezing_mean,c='aquamarine',label='Refrozen melt')     
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.grid()
    ax.set_ylim(daily_runoff_lower_lim,daily_runoff_upper_lim)
    ax.margins(x=0)

    netsnowmelt = snowmelt_mean - refreezing_mean

    ax0 = ax.twinx()
    ax0.plot(np.arange(0,len(dates)),np.cumsum(netsnowmelt),c='k',label='Cumulative runoff from snow melt')
    ax0.set_ylim(cumu_runoff_lower_lim,cumu_runoff_upper_lim)
    ax0.set_ylabel('Cumulative Runoff (m w.e.)',fontsize=14)
    ax0.margins(x=0)
    #ax0.legend(fontsize=14,loc='upper right')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
    #fig.tight_layout()