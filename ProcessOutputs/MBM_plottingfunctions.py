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
import cmocean
import gc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def load_hydrologic_year(sim,year,varname,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=False,NARRvar=False):
    '''
    Returns a given variable from Oct 1 -- Sept 30 (hydrological year)
    '''
        
    # Concatenate variable from Oct 1 of current year to Sept 30 of following year    
    if NARRvar == True:
        dates_yr1 = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='3H')
        dates_yr2 = pd.date_range(start= str(year+1) + '-01-01 00:00:00',end= str(year+1) + '-12-31 21:00:00',freq='3H')
        inMB1 = Dataset(os.path.join(NARR_INPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
        inMB2 = Dataset(os.path.join(NARR_INPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year+1) + '.nc'),'r')
    else:
        dates_yr1 = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')
        dates_yr2 = pd.date_range(start= str(year+1) + '-01-01 00:00:00',end= str(year+1) + '-12-31 21:00:00',freq='D')
        inMB1 = Dataset(os.path.join(MODEL_OUTPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
        inMB2 = Dataset(os.path.join(MODEL_OUTPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year+1) + '_' + str(sim) + '.nc'),'r')
    
    if stddev == True:
        varname = str(varname[:-4])
        
    mb_array_year1 = inMB1.variables[varname][:]
    sys.stdout.flush()
    inMB1.close()
    
    mb_array_year2 = inMB2.variables[varname][:]
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
    snowmelt_l = []
    refrozenmelt_l = []
    netsnowmelt_l = []
    glaciermelt_l = []
    SImelt_l = []
    rain_l = []
    refrozenrain_l = []
    rainrunoff_l = []
    accumulation_l = []
    snowdepth_l = []
    Ptau_l = []
    SI_l = []
    temp_l = []
    
    for year in years:
        print('Calculating runoff & mb variables for', year)
        # Get mass balance variable:
        massbal = load_hydrologic_year(sim,year,'Netbalance',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        snowmelt = load_hydrologic_year(sim,year,'Snowmelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        refrozenmelt = load_hydrologic_year(sim,year,'Refrozenmelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        netsnowmelt = load_hydrologic_year(sim,year,'Netsnowmelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        glaciermelt = load_hydrologic_year(sim,year,'Glaciericemelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        SImelt = load_hydrologic_year(sim,year,'Superimposedicemelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        rain = load_hydrologic_year(sim,year,'Rain',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        refrozenrain = load_hydrologic_year(sim,year,'Refrozenrain',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        rainrunoff = load_hydrologic_year(sim,year,'Rainrunoff',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        accumulation = load_hydrologic_year(sim,year,'Accumulation',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        snowdepth = load_hydrologic_year(sim,year,'Snowdepth',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        Ptau = load_hydrologic_year(sim,year,'Ptau',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        SI = load_hydrologic_year(sim,year,'Superimposedice',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        temp = load_hydrologic_year(sim,year,'Temperature',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,NARRvar=True)
        
        # Get timeseries from each component at a point location:
        if PointScale == True:
            mean_massbal = np.array(massbal[:,x,y])
            mean_snowmelt = np.array(snowmelt[:,x,y])
            mean_refrozenmelt = np.array(refrozenmelt[:,x,y])
            mean_netsnowmelt = np.array(netsnowmelt[:,x,y])
            mean_glaciericemelt = np.array(glaciermelt[:,x,y])
            mean_superimposedicemelt = np.array(SImelt[:,x,y])
            mean_rain = np.array(rain[:,x,y])
            mean_refrozenrain = np.array(refrozenrain[:,x,y])
            mean_rainrunoff = np.array(rainrunoff[:,x,y])
            mean_accumulation = np.array(accumulation[:,x,y])
            mean_snowdepth = np.array(snowdepth[:,x,y])
            mean_Ptau = np.array(Ptau[:,x,y])
            mean_SI = np.array(SI[:,x,y])
            mean_temp = np.array(temp[:,x,y])
            
        
        # Get timeseries from each component for the catchment wide average
        elif Catchmentaverage == True:            
            mean_massbal = np.nanmean(massbal,axis=(1,2))
            mean_snowmelt = np.nanmean(snowmelt,axis=(1,2))
            mean_refrozenmelt = np.nanmean(refrozenmelt,axis=(1,2))
            mean_netsnowmelt = np.nanmean(netsnowmelt,axis=(1,2))
            mean_glaciericemelt = np.nanmean(glaciermelt,axis=(1,2))
            mean_superimposedicemelt = np.nanmean(SImelt,axis=(1,2))
            mean_rain = np.nanmean(rain,axis=(1,2))
            mean_refrozenrain = np.nanmean(refrozenrain,axis=(1,2))
            mean_rainrunoff = np.nanmean(rainrunoff,axis=(1,2))
            mean_accumulation = np.nanmean(accumulation,axis=(1,2))
            mean_snowdepth = np.nanmean(snowdepth,axis=(1,2))
            mean_Ptau = np.nanmean(Ptau,axis=(1,2))
            mean_SI = np.nanmean(SI,axis=(1,2))
            mean_temp = np.nanmean(temp,axis=(1,2))
            
        # Get timeseries from each component for the glacier wide average
        elif KWaverage == True:
            for i in range(0,len(snowmelt)):
                massbal[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                snowmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                refrozenmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                netsnowmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                glaciermelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                SImelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                rain[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                refrozenrain[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                rainrunoff[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                accumulation[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                snowdepth[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                Ptau[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                SI[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                temp[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
    
            mean_massbal = np.nanmean(massbal,axis=(1,2))
            mean_snowmelt = np.nanmean(snowmelt,axis=(1,2))
            mean_refrozenmelt = np.nanmean(refrozenmelt,axis=(1,2))
            mean_netsnowmelt = np.nanmean(netsnowmelt,axis=(1,2))
            mean_glaciericemelt = np.nanmean(glaciermelt,axis=(1,2))
            mean_superimposedicemelt = np.nanmean(SImelt,axis=(1,2))
            mean_rain = np.nanmean(rain,axis=(1,2))
            mean_refrozenrain = np.nanmean(refrozenrain,axis=(1,2))
            mean_rainrunoff = np.nanmean(rainrunoff,axis=(1,2))
            mean_accumulation = np.nanmean(accumulation,axis=(1,2))
            mean_snowdepth = np.nanmean(snowdepth,axis=(1,2))
            mean_Ptau = np.nanmean(Ptau,axis=(1,2))
            mean_SI = np.nanmean(SI,axis=(1,2))
            mean_temp = np.nanmean(temp,axis=(1,2))
        else:
            pass
        
        massbal_l.append(mean_massbal)
        snowmelt_l.append(mean_snowmelt)
        refrozenmelt_l.append(mean_refrozenmelt)
        netsnowmelt_l.append(mean_netsnowmelt)
        glaciermelt_l.append(mean_glaciericemelt)
        SImelt_l.append(mean_superimposedicemelt)
        rain_l.append(mean_rain)
        refrozenrain_l.append(mean_refrozenrain)
        rainrunoff_l.append(mean_rainrunoff)
        accumulation_l.append(mean_accumulation)
        snowdepth_l.append(mean_snowdepth)
        Ptau_l.append(mean_Ptau)
        SI_l.append(mean_SI)
        temp_l.append(mean_temp)
        gc.collect()
        
    return massbal_l, snowmelt_l, refrozenmelt_l, netsnowmelt_l, glaciermelt_l, SImelt_l, rain_l, refrozenrain_l, rainrunoff_l, accumulation_l, snowdepth_l, Ptau_l, SI_l, temp_l

def save_timeserieslist_as_txt(years,outputlist,OUTPUTPATH,modelname,variable,domain):
    ''' 
    outputlist = the list of daily timeseries for each year (e.g. massbal_kw)
    OUTPUTPATH = where the text files should be saved
    modelname (str): the name of the model (e.g. REF_MODEL, ROUNCE_DEBRIS)
    variable (str): the name of the var (e.g. massbal, totalsnowmelt)
    domain (str): what area these outputs cover (e.g. krh, allgl, kw, nonKWice)
    '''
    mb_array = np.zeros((len(years),366))
    for year in years:
        #print(year,year-years[0])
        i = year-years[0]
        if len(outputlist[i]) == 366:
            mb_array[i,:] = outputlist[i]
        else:
            mb_array[i,:] = np.concatenate((outputlist[i],np.array([np.nan])))
            
    np.savetxt(os.path.join(OUTPUTPATH,modelname + '_' + variable + '_' + domain + '_' + str(years[0]) + '-' + str(years[-1]) + '.txt'),np.array(mb_array))

def full_save(all_outputs,years,OUTPUTPATH,modelname,varnames,domain):
    i  = 0
    for output in all_outputs:
        print(varnames[i])
        save_timeserieslist_as_txt(years,output,OUTPUTPATH,modelname,varnames[i],domain)
        i += 1
        


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
    snowmelt_l = []
    refrozenmelt_l = []
    netsnowmelt_l = []
    glaciermelt_l = []
    SImelt_l = []
    rain_l = []
    refrozenrain_l = []
    rainrunoff_l = []
    accumulation_l = []
    snowdepth_l = []
    Ptau_l = []
    SI_l = []
    temp_l = []
    
    for year in years:
        print('Calculating runoff & mb variables for', year)
        # Get mass balance variable:
        massbal = load_hydrologic_year(sim,year,'Netbalance',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        snowmelt = load_hydrologic_year(sim,year,'Snowmelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        refrozenmelt = load_hydrologic_year(sim,year,'Refrozenmelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        netsnowmelt = load_hydrologic_year(sim,year,'Netsnowmelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        glaciermelt = load_hydrologic_year(sim,year,'Glaciericemelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        SImelt = load_hydrologic_year(sim,year,'Superimposedicemelt',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        rain = load_hydrologic_year(sim,year,'Rain',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        refrozenrain = load_hydrologic_year(sim,year,'Refrozenrain',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        rainrunoff = load_hydrologic_year(sim,year,'Rainrunoff',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        accumulation = load_hydrologic_year(sim,year,'Accumulation',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        snowdepth = load_hydrologic_year(sim,year,'Snowdepth',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        Ptau = load_hydrologic_year(sim,year,'Ptau',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        SI = load_hydrologic_year(sim,year,'Superimposedice',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID)
        temp = load_hydrologic_year(sim,year,'Temperature',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,NARRvar=True)
        
        massbal_l.append(np.sum(massbal,axis=0))
        snowmelt_l.append(np.sum(snowmelt,axis=0))
        refrozenmelt_l.append(np.sum(refrozenmelt,axis=0))
        netsnowmelt_l.append(np.sum(netsnowmelt,axis=0))
        glaciermelt_l.append(np.sum(glaciermelt,axis=0))
        SImelt_l.append(np.sum(SImelt,axis=0))
        rain_l.append(np.sum(rain,axis=0))
        refrozenrain_l.append(np.sum(refrozenrain,axis=0))
        rainrunoff_l.append(np.sum(rainrunoff,axis=0))
        accumulation_l.append(np.sum(accumulation,axis=0))
        snowdepth_l.append(snowdepth[-1])
        Ptau_l.append(Ptau[-1])
        SI_l.append(SI[-1])
        temp_l.append(np.mean(temp,axis=0))
        
    return massbal_l, snowmelt_l, refrozenmelt_l, netsnowmelt_l, glaciermelt_l, SImelt_l, rain_l, refrozenrain_l, rainrunoff_l, accumulation_l, snowdepth_l, Ptau_l, SI_l, temp_l

def save_distributed_outputs_as_txt(outputs,varnames,years,OUTPUT_PATH,modelname):
    i  = 0
    for output in outputs:
        print(varnames[i])
        # Create a directory
        new_directory = os.path.join(OUTPUT_PATH,varnames[i] + '_distributed')
        os.makedirs(new_directory, exist_ok=True)
        for year in years:
            print(year)
            np.savetxt(os.path.join(new_directory,modelname + '_' + varnames[i] + '_distributed_' + str(year) + '.txt'),output[year-years[0]])
        i += 1

def runoff_timeseries_12years(units,title,year1, years, daily_runoff_upperlim, cumu_runoff_upperlim, gl_icemelt,snowmelt_runoff, rain_runoff, superimp_icemelt,area_map):
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year    
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2
    
    
    a = []
    fig, axs = plt.subplots(nrows=4, ncols=3,figsize=(12,12))
    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            a.append(ax)
    
    for i, ax in enumerate(a):
        year = i+year1
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')[:365]
        total_runoff = gl_icemelt[year-years[0]][:365] + snowmelt_runoff[year-years[0]][:365] + rain_runoff[year-years[0]][:365] + superimp_icemelt[year-years[0]][:365]
    
        ax.set_title(str(year)+'-'+str(year+1))
    
        ax.plot(np.arange(0,len(dates)),gl_icemelt[year-years[0]][:365]*dc,c='turquoise',label='Glacier ice melt')    
        ax.plot(np.arange(0,len(dates)),snowmelt_runoff[year-years[0]][:365]*dc,c='royalblue',label='Snow melt')
        ax.plot(np.arange(0,len(dates)),rain_runoff[year-years[0]][:365]*dc,c='deeppink',label='Rain')
        ax.plot(np.arange(0,len(dates)),superimp_icemelt[year-years[0]][:365]*dc,c='darkorange',label='Superimposed ice melt')
        ax.plot(np.arange(0,len(dates)),total_runoff*dc,c='k',label='Total runoff')
        ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
        ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45)
        ax.set_ylim(0,daily_runoff_upperlim)
        ax.grid()
        ax.set_xlim(182,365)
        ax.margins(x=0)
    
        #total_runoff = gl_icemelt[year-years[0]][:365] + snowmelt_runoff[year-years[0]][:365] + rain_runoff[year-years[0]][:365] + superimp_icemelt[year-years[0]][:365]
    
        #ax0 = ax.twinx()
        #ax0.plot(np.arange(0,len(dates)),np.cumsum(total_runoff)*yc,c='k',label='cumulative runoff')
        #ax0.set_ylim(0,cumu_runoff_upperlim)
        #ax0.margins(x=0)
        
        if units=='mwe':
            ax.set_ylabel('Runoff (m w.e. day$^{-1}$)',fontsize=11)
            #ax0.set_ylabel('Cumulative Runoff (m w.e. a$^{-1}$)',fontsize=11)
        elif units=='m3':
            ax.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=11)
            #ax0.set_ylabel('Cumulative Runoff (km$^3$ a$^{-1}$)',fontsize=11)        
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),fontsize=14,bbox_to_anchor=(0.98,1.08), ncol=3, borderaxespad=0.19)
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
        
        total_runoff = np.nansum(gl_icemelt[year-years[0]] + netsnowmelt[year-years[0]] + rain_runoff[year-years[0]] + superimp_icemelt[year-years[0]])
        gl_ice_percent = round((np.nansum(gl_icemelt[year-years[0]])/total_runoff)*100,1)
        snow_percent = round((np.nansum(netsnowmelt[year-years[0]])/total_runoff)*100,1)
        SIice_percent = round((np.nansum(superimp_icemelt[year-years[0]])/total_runoff)*100,1)
        rain_percent = round((np.nansum(rain_runoff[year-years[0]])/total_runoff)*100,1)
        
        plt.pie([gl_ice_percent,snow_percent,rain_percent,SIice_percent],autopct='%.1f%%',colors=piechart_colours, textprops={'fontsize': 14})
        #plt.pie([gl_ice_percent,snow_percent,rain_percent,SIice_percent],labels=[str(gl_ice_percent) + '%',str(snow_percent) + '%',str(rain_percent) + '%',str(SIice_percent) + '%'],colors=piechart_colours, textprops={'fontsize': 12})
    fig.suptitle(title,fontsize=14,y=1.01) 
    fig.legend(['Glacier ice melt','Snow melt','Rain','Superimposed ice melt'],fontsize=14,bbox_to_anchor=(1,1.03),ncol=2)
    plt.tight_layout()
     
    
def runoff_contribution_timeseries(all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):
    ''' 
    Same function as above, just calculated using the distributed totals instead,
    want to make sure I get the same result both ways
    '''
    
    Total_runoff_annual = []
    Glacier_ice_melt = []
    Snow_melt = []
    Rain = []
    SI_melt = []
    for year in all_years:
        #print(year)
        Total_runoff = netsnowmelt_dist[year-all_years[0]] + gl_icemelt_dist[year-all_years[0]] + superimp_icemelt_dist[year-all_years[0]] + rain_runoff_dist[year-all_years[0]]
        Total_runoff_annual.append(np.nansum(Total_runoff))
        
        Glacier_ice_melt.append(np.nansum(gl_icemelt_dist[year-all_years[0]]))
        Snow_melt.append(np.nansum(netsnowmelt_dist[year-all_years[0]]))
        Rain.append(np.nansum(rain_runoff_dist[year-all_years[0]]))
        SI_melt.append(np.nansum(superimp_icemelt_dist[year-all_years[0]]))
        
    #return np.array(Total_runoff_annual), np.array(Glacier_ice_melt), np.array(Snow_melt), np.array(Rain), np.array(SI_melt)
    gl_ice_percent = (np.array(Glacier_ice_melt)/np.array(Total_runoff_annual))*100
    snow_percent = (np.array(Snow_melt)/np.array(Total_runoff_annual))*100
    SIice_percent = (np.array(SI_melt)/np.array(Total_runoff_annual))*100
    rain_percent = (np.array(Rain)/np.array(Total_runoff_annual))*100
    
    fig = plt.figure(figsize=(9,3.5))
    plt.plot(all_years,gl_ice_percent,color='turquoise',label = 'Glacier ice runoff')
    plt.plot(all_years,gl_ice_percent+snow_percent,color='royalblue',label = 'Snow runoff')
    plt.plot(all_years,gl_ice_percent+snow_percent+rain_percent,color='deeppink',label = 'Rain')
    plt.plot(all_years,gl_ice_percent+snow_percent+rain_percent+SIice_percent,color='darkorange',label = 'Refrozen ice runoff')
    
    plt.fill_between(all_years,np.zeros(all_years.shape),gl_ice_percent,color='turquoise',alpha=0.35)
    plt.fill_between(all_years,gl_ice_percent,gl_ice_percent+snow_percent,color='royalblue',alpha=0.35)
    plt.fill_between(all_years,gl_ice_percent+snow_percent,gl_ice_percent+snow_percent+rain_percent,color='deeppink',alpha=0.35)
    plt.fill_between(all_years,gl_ice_percent+snow_percent+rain_percent,np.ones(all_years.shape)*100,color='darkorange',alpha=0.35)
    plt.xlim(1980,2021)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0,101)
    fig.legend(ncol=4,fontsize=14,bbox_to_anchor=(1.03,1.1))
    plt.ylabel('Contribution to\ntotal runoff (%)',fontsize=14)
    plt.grid()  
    plt.tight_layout()
    
    fig = plt.figure(figsize=(9,3.5))
    plt.hist(gl_ice_percent,all_years, histtype='step', stacked=True, fill=False,color='turquoise',label = 'Glacier ice runoff')
    plt.plot(all_years,gl_ice_percent+snow_percent,color='royalblue',label = 'Snow runoff')
    plt.plot(all_years,gl_ice_percent+snow_percent+rain_percent,color='deeppink',label = 'Rain')
    plt.plot(all_years,gl_ice_percent+snow_percent+rain_percent+SIice_percent,color='darkorange',label = 'Refrozen ice runoff')
    
    plt.fill_between(all_years,np.zeros(all_years.shape),gl_ice_percent,color='turquoise',alpha=0.35)
    plt.fill_between(all_years,gl_ice_percent,gl_ice_percent+snow_percent,color='royalblue',alpha=0.35)
    plt.fill_between(all_years,gl_ice_percent+snow_percent,gl_ice_percent+snow_percent+rain_percent,color='deeppink',alpha=0.35)
    plt.fill_between(all_years,gl_ice_percent+snow_percent+rain_percent,np.ones(all_years.shape)*100,color='darkorange',alpha=0.35)
    plt.xlim(1980,2021)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0,101)
    fig.legend(ncol=4,fontsize=14,bbox_to_anchor=(1.03,1.1))
    plt.ylabel('Contribution to\ntotal runoff (%)',fontsize=14)
    plt.grid()  
    plt.tight_layout()
    
    fig = plt.figure(figsize=(9.5,3.5))
    m, b = np.polyfit(all_years[1:],gl_ice_percent[1:],deg=1)
    plt.plot(all_years,m*all_years + b,linestyle='--',c='teal',linewidth=2)
    plt.plot(all_years,gl_ice_percent,c='turquoise',linewidth=2,label='Glacier ice melt\n(Trend = ' + f'{m:.2f}' + ' % a$^{-1}$)' )
    print('m_ice =',m)
    m, b = np.polyfit(all_years[1:],snow_percent[1:],deg=1)
    plt.plot(all_years,m*all_years + b,linestyle='--',c='navy',linewidth=2)
    plt.plot(all_years,snow_percent,c='royalblue',linewidth=2,label='Snow melt\n(Trend = ' + f'{m:.2f}' + ' % a$^{-1}$)' )
    print('m_snow =',m)
    m, b = np.polyfit(all_years[1:],rain_percent[1:],deg=1)
    plt.plot(all_years,m*all_years + b,linestyle='--',c='mediumvioletred',linewidth=2)
    plt.plot(all_years,rain_percent,c='deeppink',linewidth=2,label='Rain\n(Trend = ' + f'{m:.2f}' + ' % a$^{-1}$)')
    print('m_ice =',m)
    m, b = np.polyfit(all_years[1:],SIice_percent[1:],deg=1)
    plt.plot(all_years,m*all_years + b,linestyle='--',c='chocolate',linewidth=2)
    plt.plot(all_years,SIice_percent,c='darkorange',linewidth=2,label='Refrozen ice melt\n(Trend = ' + f'{m:.2f}' + ' % a$^{-1}$)')
    print('m_SI =',m)
    #plt.plot(years,SIice_percent,c='darkorange')
    #plt.plot(years,rain_percent,c='deeppink')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim(1980,2021)
    plt.ylim(0,70)
    #plt.ylim(20,70)
    plt.ylabel('Contribution to\ntotal runoff (%)',fontsize=14)
    fig.legend(fontsize=14,ncol=2,bbox_to_anchor=(0.8,1.25))
    plt.grid()
    #plt.savefig('D:\Model Runs\REF_MODEL\Plots\RefModel_Runoff_Contribution_Timeseries.pdf',bbox_inches='tight')
    
def glacier_contribution_to_annualrunoff(all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):
    
    Total_runoff_annual = []
    Glacierized_area_runoff_annual = []
    for year in all_years:
        #print(year)
        Total_runoff = netsnowmelt_dist[year-all_years[0]] + gl_icemelt_dist[year-all_years[0]] + superimp_icemelt_dist[year-all_years[0]] + rain_runoff_dist[year-all_years[0]]
        Total_runoff_annual.append(np.nansum(Total_runoff))
        
        #Glacierized_area_runoff = np.nansum(Total_runoff[Sfc==0])
        #Glacierized_area_runoff_annual.append(Glacierized_area_runoff)
        
        Glacierized_area_runoff_annual.append(np.nansum(gl_icemelt_dist[year-all_years[0]]))
        
    return Total_runoff_annual, Glacierized_area_runoff_annual

    
def distributed_average_mass_balance(model_name,avg_years,all_years,dist_massbal,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        print(year)
        MB_ARRAY[year-avg_years[0],:,:] = dist_massbal[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.mean((np.nanmean(MB_ARRAY,axis=0))[np.isfinite(KRH_tributaries)]),2)
    print(np.nanmin(np.nanmean(MB_ARRAY,axis=0)),np.nanmax(np.nanmean(MB_ARRAY,axis=0)))  
    
    plt.figure(figsize=(9,5))
    #plt.title(model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),cmap='massbal',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,5,2))
    legend.ax.set_ylabel('Mass Balance (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=16)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=1)
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    #plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nmass balance\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)    
    # fig.patch.set_facecolor('#f8f5f0')
    plt.axis('off')
    plt.tight_layout()
    #plt.gca().set_facecolor('#f8f5f0')
    
    return np.nanmean(MB_ARRAY,axis=0)
    
def AAR_over_time(years,dist_massbal_ref,dist_massbal_uncorrectedacc,massbal_dist_rounce):

    AARs_ref = []
    for year in years:
        dist_massbal_ref[year-years[0]][Sfc==1]=np.nan
        AAR = (np.where(dist_massbal_ref[year-years[0]] > 0)[0].shape[0])/(np.where(np.isfinite(dist_massbal_ref[year-years[0]]))[0].shape[0])
        AARs_ref.append(AAR)
        
    AARs_uncorracc = []
    for year in years:
        dist_massbal_uncorrectedacc[year-years[0]][Sfc==1]=np.nan
        AAR = (np.where(dist_massbal_uncorrectedacc[year-years[0]] > 0)[0].shape[0])/(np.where(np.isfinite(dist_massbal_uncorrectedacc[year-years[0]]))[0].shape[0])
        AARs_uncorracc.append(AAR)

    AARs_rouncedebris = []
    for year in years:
        massbal_dist_rounce[year-years[0]][Sfc==1]=np.nan
        AAR = (np.where(massbal_dist_rounce[year-years[0]] > 0)[0].shape[0])/(np.where(np.isfinite(massbal_dist_rounce[year-years[0]]))[0].shape[0])
        AARs_rouncedebris.append(AAR)
        
    AARs_uncorracc[34] = AARs_ref[34]
        
    fig = plt.figure(figsize=(7,2.5))
    plt.plot(years,AARs_ref,c='k',linewidth=2,label='Site-specific\naccumulation\ncorrection')
    #m, b = np.polyfit(years,np.array(gl_ice_percent)[:,0],deg=1)
    #plt.plot(years,m*years + b,linestyle='--',c='mediumturquoise',linewidth=2,label='Trend = ' + str(np.round(m,2)) + ' days a$^{-1}$')
    #print('m_ice =',m)
    plt.plot(years,AARs_uncorracc,c='skyblue',linewidth=2,label='Uncorrected\naccumulation')
    #m, b = np.polyfit(years,np.array(snow_percent)[:,0],deg=1)
    #plt.plot(years,m*years + b,linestyle='--',c='navy',linewidth=2,label='Trend = ' + str(np.round(m,2)) + ' days a$^{-1}$')
    #print('m_snow =',m)
    #plt.plot(years,AARs_rouncedebris,c='saddlebrown',linewidth=2,label='Standard debris treatment (global dataset)')
    #plt.plot(years,SIice_percent,c='darkorange')
    #plt.plot(years,rain_percent,c='deeppink')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0.20,0.9)
    plt.ylabel('Accumulation Area Ratio',fontsize=14)
    plt.legend(fontsize=16,bbox_to_anchor=(1,1))
    plt.grid()
    
    fig.patch.set_facecolor('#f8f5f0')
    
    # fig.patch.set_facecolor('#f8f5f0')
    #plt.axis('off')
    #plt.gca().set_facecolor('#f8f5f0')
    return AARs_ref, AARs_uncorracc, AARs_rouncedebris
    

def distributed_runoff(avg_years,all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] = netsnowmelt_dist[year-all_years[0]] + gl_icemelt_dist[year-all_years[0]] + superimp_icemelt_dist[year-all_years[0]] + rain_runoff_dist[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.nanmean(np.nanmean(MB_ARRAY,axis=0)),2)
    
    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))
    
    fig = plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,Runoff,cmap=cmocean.cm.ice_r,levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.linspace(0,20,11))
    legend.ax.set_ylabel('Total runoff (m w.e. a$^{-1}$)', rotation=270,fontsize=16,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=16)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    #plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nMean runoff\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    fig.patch.set_facecolor('#f8f5f0')
    plt.axis('off')
    plt.gca().set_facecolor('#f8f5f0')
    plt.tight_layout()
    
def distributed_glaciericemelt(avg_years,all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] = gl_icemelt_dist[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.nanmean(np.nanmean(MB_ARRAY,axis=0)),2)

    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))
    
    plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,Runoff,cmap=cmocean.cm.ice_r,levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.linspace(0,20,11))
    legend.ax.set_ylabel('Runoff from glacier ice melt (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nMean glacier ice melt\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    plt.tight_layout()
    
def distributed_snowrunoff(avg_years,all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] =  netsnowmelt_dist[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.nanmean(np.nanmean(MB_ARRAY,axis=0)),2)
    
    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))
    
    plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,Runoff,cmap=cmocean.cm.ice_r,levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.linspace(0,1,11))
    legend.ax.set_ylabel('Runoff from snow melt (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nMean snow melt\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    plt.tight_layout()
    
def distributed_rainrunoff(avg_years,all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] =  rain_runoff_dist[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.nanmean(np.nanmean(MB_ARRAY,axis=0)),2)
    
    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))
    
    plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,Runoff,cmap=cmocean.cm.ice_r,levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.linspace(0,1,11))
    legend.ax.set_ylabel('Runoff from rain (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nMean rain runoff\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    plt.tight_layout()
    
def distributed_SImelt(avg_years,all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] =  superimp_icemelt_dist[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.nanmean(np.nanmean(MB_ARRAY,axis=0)),2)
    
    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))
    
    plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,Runoff,cmap=cmocean.cm.ice_r,levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.linspace(0,0.1,11))
    legend.ax.set_ylabel('Runoff from superimposed ice (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nMean superimposed ice melt\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    plt.tight_layout()
    
def distributed_runoff_components(avg_years,all_years,Sfc,runoff_component):
    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] =  runoff_component[year-all_years[0]]

    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))

    return Runoff   

def distributed_runoff_allcomponents(years,Sfc,Xgrid,Ygrid,gl_icemelt_dist_ref,netsnowmelt_dist_ref,rain_runoff_dist_ref,superimp_icemelt_dist_ref):
    fig, axs = plt.subplots(2, 2, figsize=(9,8))
    
    # Plot 1 with inset
    icemelt = distributed_runoff_components(np.arange(1980,2021+1),years,Sfc,gl_icemelt_dist_ref)
    contour1 = axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18))
    axs[0, 0].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[0, 0].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[0, 0].axis('off')
    axs[0, 0].axis('equal')
    axs[0, 0].set_title('a)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider1 = make_axes_locatable(axs[0, 0])
    cax1 = divider1.append_axes("bottom", size="5%", pad=0.1)
    cb1 = plt.colorbar(contour1, cax=cax1, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb1.set_label('Runoff from glacier ice (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb1.set_ticks(np.linspace(0,20,11))
    cb1.ax.tick_params(labelsize=14)
    
    snowmelt = distributed_runoff_components(np.arange(1980,2021+1),years,Sfc,netsnowmelt_dist_ref)
    contour2 = axs[0, 1].contourf(Xgrid,Ygrid,snowmelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,0.7,18))
    axs[0, 1].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[0, 1].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[0, 1].axis('off')
    axs[0, 1].axis('equal')
    axs[0, 1].set_title('b)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider2 = make_axes_locatable(axs[0, 1])
    cax2 = divider2.append_axes("bottom", size="5%", pad=0.1)
    cb2 = plt.colorbar(contour2, cax=cax2, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb2.set_label('Runoff from snow (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb2.set_ticks(np.linspace(0,1,11))
    cb2.ax.tick_params(labelsize=14)
    
    rain = distributed_runoff_components(np.arange(1980,2021+1),years,Sfc,rain_runoff_dist_ref)
    contour3 = axs[1, 0].contourf(Xgrid,Ygrid,rain,cmap=cmocean.cm.ice_r,levels=np.linspace(0,0.3,18))
    axs[1, 0].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[1, 0].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[1, 0].axis('off')
    axs[1, 0].axis('equal')
    axs[1, 0].set_title('c)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider3 = make_axes_locatable(axs[1, 0])
    cax3 = divider3.append_axes("bottom", size="5%", pad=0.1)
    cb3 = plt.colorbar(contour3, cax=cax3, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb3.set_label('Runoff from rain (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb3.set_ticks(np.arange(0,0.4,0.05))
    cb3.ax.tick_params(labelsize=14)
    
    SI = distributed_runoff_components(np.arange(1980,2021+1),years,Sfc,superimp_icemelt_dist_ref)
    contour4 = axs[1, 1].contourf(Xgrid,Ygrid,SI,cmap=cmocean.cm.ice_r,levels=np.linspace(0,0.1,18))
    axs[1, 1].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[1, 1].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[1, 1].axis('off')
    axs[1, 1].axis('equal')
    axs[1, 1].set_title('d)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider4 = make_axes_locatable(axs[1, 1])
    cax4 = divider4.append_axes("bottom", size="5%", pad=0.1)
    cb4 = plt.colorbar(contour4, cax=cax4, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb4.set_label('Runoff from refrozen ice (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb4.set_ticks(np.arange(0,0.2,0.02))
    cb4.ax.tick_params(labelsize=14)
    
    #plt.savefig('D:\Model Runs\REF_MODEL\Plots\Distributed_runoff_components.pdf',bbox_inches='tight')
 

def massbalance_timeseries_12years(title,year1, years, daily_mb_abs_lim, cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt):
    a = []
    fig, axs = plt.subplots(nrows=4, ncols=3,figsize=(12,12))
    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            a.append(ax)
    
    for i, ax in enumerate(a):
        year = i+year1
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')[:365]
        
        ax.set_title(str(year)+'-'+str(year+1))
        ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)')
    
        total_accumulation = accumulation[year-years[0]][:365] + refrozen_rain[year-years[0]][:365]
        total_ablation = netsnowmelt[year-years[0]][:365] + superimp_icemelt[year-years[0]][:365] + gl_icemelt[year-years[0]][:365]
    
        ax.plot(np.arange(0,len(dates)),total_accumulation,c='mediumblue',label='Accumulation')    
        ax.plot(np.arange(0,len(dates)),-total_ablation,c='red',label='Ablation')    
        ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
        ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45)
        ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
        ax.grid()
        ax.tick_params(axis='y',labelsize=14)
        ax.margins(x=0)
        
        massbal = (accumulation[year-years[0]][:365] + refrozen_rain[year-years[0]][:365]) - (netsnowmelt[year-years[0]][:365] + superimp_icemelt[year-years[0]][:365] + gl_icemelt[year-years[0]][:365])
    
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
        
        print(year,np.cumsum(massbal)[-1])
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),fontsize=14,bbox_to_anchor=(0.98,1.04), ncol=2, borderaxespad=0.19)
    fig.suptitle(title,fontsize=14,y=1.01) 
    fig.tight_layout()
    
    
def massbalance_timeseries_average(title,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt, accumulation_std, refrozen_rain_std, netsnowmelt_std, superimp_icemelt_std, gl_icemelt_std,area_map):
    
    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation[i][:365])
        refrain_sum += np.array(refrozen_rain[i][:365])
        snowmelt_sum += np.array(netsnowmelt[i][:365])
        SImelt_sum += np.array(superimp_icemelt[i][:365])
        gl_melt_sum += np.array(gl_icemelt[i][:365]) 
        
        snowfallstd_sum += np.array(accumulation_std[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std[i][:365]) 
        
    snowfall_mean = np.array(snowfall_sum/len(avg_years))
    refrain_mean = np.array(refrain_sum/len(avg_years))
    snowmelt_mean = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean = np.array(gl_melt_sum/len(avg_years))
      
    snowfallstd_mean = np.array(snowfallstd_sum/len(avg_years))
    refrainstd_mean = np.array(refrainstd_sum/len(avg_years))
    snowmeltstd_mean = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean = np.array(gl_meltstd_sum/len(avg_years))

    total_accumulation = snowfall_mean + refrain_mean
    total_accumulation_std = snowfallstd_mean + refrainstd_mean
    
    total_ablation = snowmelt_mean + SImelt_mean + gl_melt_mean
    total_ablation_std = snowmeltstd_mean + SImeltstd_mean + gl_meltstd_mean

    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')   
    time = np.arange(0,len(dates))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9

    # Plot the average
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(9,3.5))
    #ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=15)
    ax.plot(time,total_accumulation,c='mediumblue',label='Accumulation')    
    ax.plot(time,-total_ablation,c='red',label='Ablation')    
    plt.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    plt.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax.grid()
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    
    massbal = total_accumulation - total_ablation
    max_mb = (total_accumulation + total_accumulation_std) - (total_ablation - total_ablation_std)
    min_mb = (total_accumulation - total_accumulation_std) - (total_ablation + total_ablation_std)

    transition_date = dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]]
    ax0 = ax.twinx()
    ax0.plot(time,np.cumsum(massbal)*area,c='k',label='Cumulative balance\n$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    print(np.cumsum(massbal)*area)
    #plt.fill_between(time,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='std. dev.')

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass Change (Gt yr$^{-1}$)',fontsize=15)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.legend(fontsize=15,loc='lower left')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=15,loc='upper left', ncol=1, borderaxespad=0.19)
    fig.patch.set_facecolor('#f8f5f0')
    fig.tight_layout()

def cumulative_massbalance(title,years,plotting_years,daily_mb_abs_lim,cumu_mb_abs_lim,accumulation,refrozen_rain,gl_icemelt,netsnowmelt,SImelt,area_map):
    
    year_i_index = plotting_years[0] - years[0]
    year_f_index = len(years) - (years[-1] - plotting_years[-1])
    #year_i_index = 0
    #year_f_index = 43
    
    all_accumulation = np.concatenate(accumulation[year_i_index:year_f_index], axis=0) + np.concatenate(refrozen_rain[year_i_index:year_f_index], axis=0)
    all_ablation = np.concatenate(gl_icemelt[year_i_index:year_f_index], axis=0) + np.concatenate(netsnowmelt[year_i_index:year_f_index], axis=0) + np.concatenate(SImelt[year_i_index:year_f_index], axis=0)
    
    acc_final = all_accumulation[~np.isnan(all_accumulation)]
    abl_final = all_ablation[~np.isnan(all_ablation)]
    
    massbal = acc_final - abl_final    
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9
    
    dates = pd.date_range(start= str(plotting_years[0]) + '-10-01 00:00:00',end= str(plotting_years[-1]+1) + '-09-30 21:00:00',freq='1D')   
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(9,4.5))
    ax.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=14)
    ax.plot(dates,acc_final,c='mediumblue',label='Accumulation')    
    ax.plot(dates,-abl_final,c='red',label='Ablation')  
    #ax.plot(np.arange(-10,-8),np.arange(-10,-8),c='k',label='Cumulative Mass Balance',linewidth=3)  
    #plt.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    #plt.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    #ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    #ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax.grid()
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.tick_params(axis='x',labelsize=14)
    ax.margins(x=0)
    
    ##massbal = total_accumulation - total_ablation
    #max_mb = (total_accumulation + total_accumulation_std) - (total_ablation - total_ablation_std)
    #min_mb = (total_accumulation - total_accumulation_std) - (total_ablation + total_ablation_std)

    #transition_date = dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]]
    ax0 = ax.twinx()
    ax0.plot(dates,np.cumsum(massbal)*area,c='k',linewidth=3,label='Cumulative Mass Balance')
    print((np.cumsum(massbal)*area)[-1])
    #plt.fill_between(dates,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(dates,-12,-10,color='grey',alpha=0.35,label='std. dev.')

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass Balance (Gt)',fontsize=14)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.legend(fontsize=14,loc='lower left')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
    fig.tight_layout()
    
    cumu_mb_Gt = (np.cumsum(massbal)*area)[-1]
    cumu_mb_mwe =  np.cumsum(massbal)[-1]
    
    return [cumu_mb_Gt,cumu_mb_mwe]
    
def annual_massbalance_barchart(title,years,plotting_years,mb_abs_lim,accumulation,refrozen_rain,gl_icemelt,netsnowmelt,SImelt):
    
    year_i_index = plotting_years[0] - years[0]
    year_f_index = len(years) - (years[-1] - plotting_years[-1])
    #year_i_index = 0
    #year_f_index = 43
    
    all_accumulation = np.nansum(accumulation[year_i_index:year_f_index], axis=1) + np.nansum(refrozen_rain[year_i_index:year_f_index], axis=1)
    all_ablation = np.nansum(gl_icemelt[year_i_index:year_f_index], axis=1) + np.nansum(netsnowmelt[year_i_index:year_f_index], axis=1) + np.nansum(SImelt[year_i_index:year_f_index], axis=1)
    
    #acc_final = all_accumulation[~np.isnan(all_accumulation)]
    #abl_final = all_ablation[~np.isnan(all_ablation)]
    
    massbal = all_accumulation - all_ablation    
    
    #area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9
    #dates = pd.date_range(start= str(plotting_years[0]) + '-10-01 00:00:00',end= str(plotting_years[-1]+1) + '-09-30 21:00:00',freq='1D')   
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(9,4.5))
    ax.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Balance (m w.e. a$^{-1}$)',fontsize=14)
    ax.bar(plotting_years,all_accumulation,color='mediumblue',label='Accumulation',zorder=10)    
    ax.bar(plotting_years,-all_ablation,color='red',label='Ablation',zorder=10)  
    #ax.plot(np.arange(-10,-8),np.arange(-10,-8),c='k',label='Cumulative Mass Balance',linewidth=3)  
    #plt.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    #plt.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    #ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    #ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(-mb_abs_lim,mb_abs_lim)
    ax.grid(zorder=-10)
    ax.axhline(y=0,xmin=0,xmax=len(plotting_years),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.tick_params(axis='x',labelsize=14)
    ax.margins(x=0.01)
    
    ax.bar(plotting_years,massbal,color='k',linewidth=3,label='Mass Balance',zorder=20)
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
    fig.tight_layout()
    
    return plotting_years,massbal, all_accumulation, -all_ablation
    
def cumulative_and_annual_massbal(title,years,plotting_years,daily_mb_abs_lim,cumu_mb_abs_lim,mb_abs_lim,accumulation,refrozen_rain,gl_icemelt,netsnowmelt,SImelt,area_map):
    year_i_index = plotting_years[0] - years[0]
    year_f_index = len(years) - (years[-1] - plotting_years[-1])
    
    all_accumulation = np.concatenate(accumulation[year_i_index:year_f_index], axis=0) + np.concatenate(refrozen_rain[year_i_index:year_f_index], axis=0)
    all_ablation = np.concatenate(gl_icemelt[year_i_index:year_f_index], axis=0) + np.concatenate(netsnowmelt[year_i_index:year_f_index], axis=0) + np.concatenate(SImelt[year_i_index:year_f_index], axis=0)
    
    acc_final = all_accumulation[~np.isnan(all_accumulation)]
    abl_final = all_ablation[~np.isnan(all_ablation)]
    
    massbal = acc_final - abl_final    
    
    dates = pd.date_range(start= str(plotting_years[0]) + '-10-01 00:00:00',end= str(plotting_years[-1]+1) + '-09-30 21:00:00',freq='1D')   
    filler_dates = len(dates) - (((plotting_years.repeat(365).reshape(len(plotting_years),365)) + np.linspace(0,0.995,365)).flatten()).shape[0]
    datess = np.concatenate((((plotting_years.repeat(365).reshape(len(plotting_years),365)) + np.linspace(0,0.995,365)).flatten(),np.arange(plotting_years[-1]+1,plotting_years[-1]+1+(filler_dates*0.002733),0.002733)))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,figsize=(10,7))
    #ax1.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax1.set_ylabel('Mass Change\n(m w.e. d$^{-1}$)',fontsize=14)
    ax1.plot(datess,acc_final,c='dodgerblue',label='Accumulation')    
    ax1.plot(datess,-abl_final,c='crimson',label='Ablation')  
    ax1.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax1.grid()
    ax1.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax1.tick_params(axis='y',labelsize=14)
    ax1.tick_params(axis='x',labelsize=14)
    #ax1.margins(x=0)
    ax1.text(1980,0.06,'a)',fontsize=14,weight='bold')
    ax1.set_xlim(1979,2022)
    
    ax0 = ax1.twinx()
    ax0.plot(datess,np.cumsum(massbal)*area,c='k',linewidth=3,label='Cumulative Mass Balance')
    print((np.cumsum(massbal)*area)[-1])

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass\nBalance (Gt)',fontsize=14)
    #ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.legend(fontsize=14,loc='lower left')
        
    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    
    ax1.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper right', ncol=2, borderaxespad=0.19)
    
    all_accumulation = np.nansum(accumulation[year_i_index:year_f_index], axis=1) + np.nansum(refrozen_rain[year_i_index:year_f_index], axis=1)
    all_ablation = np.nansum(gl_icemelt[year_i_index:year_f_index], axis=1) + np.nansum(netsnowmelt[year_i_index:year_f_index], axis=1) + np.nansum(SImelt[year_i_index:year_f_index], axis=1)
    massbal = all_accumulation - all_ablation 
    
    #ax2.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax2.set_ylabel('Mass Balance\n(m w.e. a$^{-1}$)',fontsize=14)
    ax2.bar(plotting_years,all_accumulation,color='dodgerblue',label='Accumulation',zorder=10)    
    ax2.bar(plotting_years,-all_ablation,color='crimson',label='Ablation',zorder=10)  
    ax2.set_ylim(-mb_abs_lim,mb_abs_lim)
    ax2.grid(zorder=-10)
    ax2.axhline(y=0,xmin=0,xmax=len(plotting_years),linestyle='--',c='k')
    ax2.tick_params(axis='y',labelsize=14)
    ax2.tick_params(axis='x',labelsize=14)
    #ax2.margins(x=0.01)
    ax2.text(1980,2,'b)',fontsize=14,weight='bold')
    ax2.set_xlim(1979,2022)
    
    ax3 = ax2.twinx()
    #ax3.bar(plotting_years,massbal*area*0,color='white',linewidth=3,label='Mass Balance',zorder=-20)
    ax3.set_ylim(-mb_abs_lim*area,mb_abs_lim*area)
    ax3.set_yticks(np.round([-1*area,-2*area,0,1*area,2*area],1))
    ax3.tick_params(axis='y',labelsize=14)
    ax3.set_ylabel('Mass Balance\n(Gt a$^{-1}$)',fontsize=14)
    
    ax2.bar(plotting_years,massbal,color='k',linewidth=3,label='Mass Balance',zorder=30)
    #ax3.axhline(y=0,xmin=0,xmax=len(plotting_years),linestyle='--',c='k')    
    
    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax2.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper right', ncol=3, borderaxespad=0.19)
    fig.tight_layout()
    
    return all_accumulation, -all_ablation, massbal
    
def annual_massbalance_barchart_5yraverages(title,years,plotting_years,mb_abs_lim,accumulation,refrozen_rain,gl_icemelt,netsnowmelt,SImelt):
    
    year_i_index = plotting_years[0] - years[0]
    year_f_index = len(years) - (years[-1] - plotting_years[-1])
    #year_i_index = 0
    #year_f_index = 43
    
    all_accumulation = np.nansum(accumulation[year_i_index:year_f_index], axis=1) + np.nansum(refrozen_rain[year_i_index:year_f_index], axis=1)
    all_ablation = np.nansum(gl_icemelt[year_i_index:year_f_index], axis=1) + np.nansum(netsnowmelt[year_i_index:year_f_index], axis=1) + np.nansum(SImelt[year_i_index:year_f_index], axis=1)
    
    avg_accumulation = np.mean(all_accumulation.reshape(-1, 5), axis=1)
    avg_ablation = np.mean(all_ablation.reshape(-1, 5), axis=1)
    #acc_final = all_accumulation[~np.isnan(all_accumulation)]
    #abl_final = all_ablation[~np.isnan(all_ablation)]
    
    massbal = all_accumulation - all_ablation
    avg_massbal =  np.mean(massbal.reshape(-1, 5), axis=1)
    
    #dates = pd.date_range(start= str(plotting_years[0]) + '-10-01 00:00:00',end= str(plotting_years[-1]+1) + '-09-30 21:00:00',freq='1D')   
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    ax.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Balance (m w.e. a$^{-1}$)',fontsize=14)
    ax.bar(np.arange(0,len(avg_accumulation)),avg_accumulation,color='mediumblue',label='Accumulation',zorder=10)    
    ax.bar(np.arange(0,len(avg_accumulation)),-avg_ablation,color='red',label='Ablation',zorder=10)  
    #ax.plot(np.arange(-10,-8),np.arange(-10,-8),c='k',label='Cumulative Mass Balance',linewidth=3)  
    #plt.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    #plt.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    ax.set_xticks(ticks=np.arange(0,len(avg_accumulation)))
    ax.set_xticklabels(['1980-\n1985','1985-\n1990','1990-\n1995','1995-\n2000','2000-\n2005','2005-\n2010','2010-\n2015','2015-\n2020'],rotation=0,fontsize=14)
    ax.set_ylim(-mb_abs_lim,mb_abs_lim)
    ax.grid(zorder=-10)
    ax.axhline(y=0,xmin=0,xmax=len(plotting_years),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.tick_params(axis='x',labelsize=14)
    ax.margins(x=0.01)
    
    ax.bar(np.arange(0,len(avg_accumulation)),avg_massbal,color='k',linewidth=3,label='Mass Balance',zorder=20)
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
    fig.tight_layout()
    
def annual_massbalance_barchart_10yraverages(title,years,plotting_years,mb_abs_lim,accumulation,refrozen_rain,gl_icemelt,netsnowmelt,SImelt):
    
    year_i_index = plotting_years[0] - years[0]
    year_f_index = len(years) - (years[-1] - plotting_years[-1])
    #year_i_index = 0
    #year_f_index = 43
    
    all_accumulation = np.nansum(accumulation[year_i_index:year_f_index], axis=1) + np.nansum(refrozen_rain[year_i_index:year_f_index], axis=1)
    all_ablation = np.nansum(gl_icemelt[year_i_index:year_f_index], axis=1) + np.nansum(netsnowmelt[year_i_index:year_f_index], axis=1) + np.nansum(SImelt[year_i_index:year_f_index], axis=1)
    
    avg_accumulation = np.mean(all_accumulation.reshape(-1, 10), axis=1)
    avg_ablation = np.mean(all_ablation.reshape(-1, 10), axis=1)
    #acc_final = all_accumulation[~np.isnan(all_accumulation)]
    #abl_final = all_ablation[~np.isnan(all_ablation)]
    
    massbal = all_accumulation - all_ablation
    avg_massbal =  np.mean(massbal.reshape(-1, 10), axis=1)
    
    #dates = pd.date_range(start= str(plotting_years[0]) + '-10-01 00:00:00',end= str(plotting_years[-1]+1) + '-09-30 21:00:00',freq='1D')   
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,5))
    ax.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Balance (m w.e. a$^{-1}$)',fontsize=14)
    ax.bar(np.arange(0,len(avg_accumulation)),avg_accumulation,color='mediumblue',label='Accumulation',zorder=10)    
    ax.bar(np.arange(0,len(avg_accumulation)),-avg_ablation,color='red',label='Ablation',zorder=10)  
    #ax.plot(np.arange(-10,-8),np.arange(-10,-8),c='k',label='Cumulative Mass Balance',linewidth=3)  
    #plt.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    #plt.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    ax.set_xticks(ticks=np.arange(0,len(avg_accumulation)))
    ax.set_xticklabels(['1980-\n1990','1990-\n2000','2000-\n2010','2010-\n2020'],rotation=0,fontsize=14)
    ax.set_ylim(-mb_abs_lim,mb_abs_lim)
    ax.grid(zorder=-10)
    ax.axhline(y=0,xmin=0,xmax=len(plotting_years),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.tick_params(axis='x',labelsize=14)
    ax.margins(x=0.01)
    
    ax.bar(np.arange(0,len(avg_accumulation)),avg_massbal,color='k',linewidth=3,label='Mass Balance',zorder=20)
        
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
    
def annualrunoff_stackedbar(title,plotting_years,years,netsnowmelt,gl_icemelt,superimp_icemelt,rain_runoff,area_map):
    # get total runoff from each component for all years:
    yrly_snowmelt = []
    yrly_icemelt = []
    yrly_SImelt = []
    yrly_rain = []
    for year in plotting_years:
        i = year-years[0]
        # get total of each component
        yrly_snowmelt.append(np.nansum(netsnowmelt[i]))
        yrly_icemelt.append(np.nansum(gl_icemelt[i]))
        yrly_SImelt.append(np.nansum(superimp_icemelt[i]))
        yrly_rain.append(np.nansum(rain_runoff[i]))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9
    
    total_runoff = (np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt) + np.array(yrly_icemelt))*area
    m, b = np.polyfit(plotting_years,total_runoff,deg=1)
    print('Trend in total runoff = ' + str(np.round(m,2))+ ' km$^{3}$ a$^{-1}$')
    
    width = 0.75       # the width of the bars: can also be len(x) sequence
    
    fig, ax = plt.subplots(figsize=(10,4))
    ax.bar(plotting_years, np.array(yrly_rain)*area, width, label='Rain',color='deeppink',zorder=5)
    ax.bar(plotting_years, np.array(yrly_SImelt)*area, width, bottom=np.array(yrly_rain)*area,label='Refrozen ice melt',color='darkorange',zorder=5)
    ax.bar(plotting_years, np.array(yrly_snowmelt)*area, width, bottom=(np.array(yrly_rain) + np.array(yrly_SImelt))*area,label='Snow melt',color='royalblue',zorder=5)
    ax.bar(plotting_years, np.array(yrly_icemelt)*area, width, bottom=(np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt))*area,label='Glacier ice melt',color='turquoise',zorder=5)
    #plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2,zorder=20,label='Trend = +' + str(np.round(m,2))+ ' km$^{3}$ a$^{-1}$')
    plt.plot(plotting_years,m*plotting_years + b,linestyle='--',c='k',linewidth=2,zorder=20)
    ax.set_xlim(1979,2022)
    
    ax.set_ylabel('Runoff (Gt a$^{-1}$)',fontsize=14)
    plt.xticks(np.arange(plotting_years[0],plotting_years[-1]+1,5),fontsize=14,rotation=0)
    plt.yticks(fontsize=14)
    plt.ylim(0,3.75)
    #ax.set_title('Monthly Runoff (2007-2018) \n Debris Case',fontsize = 15)
    ax.legend(fontsize=14,ncol=4)
    plt.grid(zorder=-30)

    plt.tight_layout() 
    
    return total_runoff

def mb_vs_runoff(years,netsnowmelt,gl_icemelt,superimp_icemelt,rain_runoff,area_map,accumulation,refrozen_rain,):
    year_i_index = 0
    year_f_index = len(years)
    
    # get total runoff from each component for all years:
    yrly_snowmelt = []
    yrly_icemelt = []
    yrly_SImelt = []
    yrly_rain = []
    for year in years:
        i = year-years[0]
        # get total of each component
        yrly_snowmelt.append(np.nansum(netsnowmelt[i]))
        yrly_icemelt.append(np.nansum(gl_icemelt[i]))
        yrly_SImelt.append(np.nansum(superimp_icemelt[i]))
        yrly_rain.append(np.nansum(rain_runoff[i]))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9
    
    total_runoff = (np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt) + np.array(yrly_icemelt))*area
    
    all_accumulation = (np.nansum(accumulation[year_i_index:year_f_index], axis=1) + np.nansum(refrozen_rain[year_i_index:year_f_index], axis=1))*area
    all_ablation = (np.nansum(gl_icemelt[year_i_index:year_f_index], axis=1) + np.nansum(netsnowmelt[year_i_index:year_f_index], axis=1) + np.nansum(superimp_icemelt[year_i_index:year_f_index], axis=1))*area
    massbal = all_accumulation - all_ablation 
    
    fig, axs = plt.subplots(3, 1, figsize=(6, 10))
    scatter1 = axs[0].scatter(total_runoff,all_accumulation, c=years, cmap='gist_heat_r')
    axs[0].set_xlabel('Runoff (Gt a$^{-1}$)',fontsize=14)
    axs[0].set_ylabel('Accumulation (Gt a$^{-1}$)',fontsize=14)
    m, b = np.polyfit(total_runoff,all_accumulation,deg=1)
    axs[0].plot(total_runoff, (m*total_runoff)+b, c='k', linestyle='-')
    print(m)
    
    axs[1].scatter(total_runoff,all_ablation, c=years, cmap='gist_heat_r')
    axs[1].set_xlabel('Runoff (Gt a$^{-1}$)',fontsize=14)
    axs[1].set_ylabel('Ablation (Gt a$^{-1}$)',fontsize=14)
    m, b = np.polyfit(total_runoff,all_ablation,deg=1)
    axs[1].plot(total_runoff, (m*total_runoff)+b, c='k', linestyle='-')
    print(m)
    
    axs[2].scatter(total_runoff,massbal, c=years, cmap='gist_heat_r')
    axs[2].set_xlabel('Runoff (Gt a$^{-1}$)',fontsize=14)
    axs[2].set_ylabel('Mass Balance (Gt a$^{-1}$)',fontsize=14)
    m, b = np.polyfit(total_runoff,massbal,deg=1)
    axs[2].plot(total_runoff, (m*total_runoff)+b, c='k', linestyle='-')
    print(m)
    axs[2].grid()
    
    cax = fig.add_axes([0, 0.1, 0.2, 0.02])  # [x_position, y_position, width, height]

    # Add a single colorbar on the right side using the separate axes
    cbar = fig.colorbar(scatter1, cax=cax, orientation='horizontal')
    
    fig.tight_layout()
    
    plt.figure()
    
    # left off here
    return total_runoff, all_accumulation, all_ablation

    
def date_of_zero_balance(years,mb):
    # transition date
    transition_dates = []
    transition_DOYs = []
    for year in years:
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
        massbal = mb[year-years[0]]
        
        d_tr = [np.where(np.cumsum(massbal)[50:] <=0)][0][0]
        if len(d_tr) == 0:
            transition_date = 'N/A'
            transition_DOY = np.nan
        else:
            transition_date = (dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]])
            transition_DOY = np.where(np.cumsum(massbal)[50:] <=0)[0][0]
        
        transition_dates.append(transition_date)
        transition_DOYs.append(transition_DOY+50)
        
    monthday = []
    for i in np.arange(int(np.nanmin(transition_DOYs)),int(np.nanmax(transition_DOYs))+11,10):
        monthday.append( str(dates[i])[5:10])
          
    m, b = np.polyfit(years[1:][np.isfinite(transition_DOYs[1:])],np.array(transition_DOYs[1:])[np.isfinite(transition_DOYs[1:])],deg=1)
    
    fig = plt.figure(figsize=(9,3.5))
    #plt.title(title + '\nDate of zero balance',fontsize=14)
    plt.scatter(years[1:],transition_DOYs[1:],c='darkblue',s=60)
    plt.yticks(np.arange(int(np.nanmin(transition_DOYs)),int(np.nanmax(transition_DOYs))+11,10),monthday,fontsize=14)
    plt.ylim(int(np.nanmin(transition_DOYs))-10,int(np.nanmax(transition_DOYs))+10)
    plt.xticks(fontsize=14)
    plt.ylabel('date of $\dot{B}$ = 0',fontsize=14)
    plt.grid()
    #fig.patch.set_facecolor('#f8f5f0')
    plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2,label='Trend = ' + str(np.round(m,2)) + ' days a$^{-1}$')
    plt.legend(fontsize=14,loc='upper right')
    
    width=8.5
    for year in years:
        #print(year,transition_DOYs[year-years[0]])
        if np.isnan(transition_DOYs[year-years[0]]):
            plt.vlines(year, ymin=0, ymax=365, linewidth=width, color='k',alpha=0.3)

def date_of_zero_balance_with_mbcurve(b0_years,mb,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries,years,avg_years,daily_mb_abs_lim,cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt, accumulation_std, refrozen_rain_std, netsnowmelt_std, superimp_icemelt_std, gl_icemelt_std,area_map):
    # transition date
    transition_dates = []
    transition_DOYs = []
    for year in b0_years:
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
        massbal = mb[year-b0_years[0]]
        
        d_tr = [np.where(np.cumsum(massbal)[50:] <=0)][0][0]
        if len(d_tr) == 0:
            transition_date = 'N/A'
            transition_DOY = np.nan
        else:
            transition_date = (dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]])
            transition_DOY = np.where(np.cumsum(massbal)[50:] <=0)[0][0]
        
        transition_dates.append(transition_date)
        transition_DOYs.append(transition_DOY+50)
        
    monthday = []
    for i in np.arange(int(np.nanmin(transition_DOYs)),int(np.nanmax(transition_DOYs))+11,10):
        monthday.append( str(dates[i])[5:10])
          
    m, b = np.polyfit(b0_years[1:][np.isfinite(transition_DOYs[1:])],np.array(transition_DOYs[1:])[np.isfinite(transition_DOYs[1:])],deg=1)
    
    #fig = plt.figure(figsize=(9,3.5))
    fig, (a2, a1) = plt.subplots(2,1,gridspec_kw={'height_ratios':[1,0.8]},figsize=(9,6.9))
    #plt.title(title + '\nDate of zero balance',fontsize=14)
    a1.scatter(b0_years[1:],transition_DOYs[1:],c='darkblue',s=60)
    a1.set_yticks(np.arange(int(np.nanmin(transition_DOYs)),int(np.nanmax(transition_DOYs))+11,10))
    a1.set_yticklabels(monthday)
    a1.tick_params(axis='both',labelsize=14)
    a1.set_ylim(int(np.nanmin(transition_DOYs))-10,int(np.nanmax(transition_DOYs))+10)
    a1.set_ylabel('date of $\dot{B}$ = 0',fontsize=14)
    a1.grid()
    a1.text(1978,347,'b)',fontsize=14,weight='bold')
    #fig.patch.set_facecolor('#f8f5f0')
    a1.plot(b0_years,m*b0_years + b,linestyle='--',c='k',linewidth=2,label='Trend = ' + str(np.round(m,2)) + ' days a$^{-1}$')
    a1.legend(fontsize=14,loc='upper right')
    
    width=8.5
    for year in b0_years:
        #print(year,transition_DOYs[year-years[0]])
        if np.isnan(transition_DOYs[year-b0_years[0]]):
            a1.vlines(year, ymin=0, ymax=365, linewidth=width, color='k',alpha=0.3)

    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - years[0]
        snowfall_sum += np.array(accumulation[i][:365])
        refrain_sum += np.array(refrozen_rain[i][:365])
        snowmelt_sum += np.array(netsnowmelt[i][:365])
        SImelt_sum += np.array(superimp_icemelt[i][:365])
        gl_melt_sum += np.array(gl_icemelt[i][:365]) 
        
        snowfallstd_sum += np.array(accumulation_std[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std[i][:365]) 
        
    snowfall_mean = np.array(snowfall_sum/len(avg_years))
    refrain_mean = np.array(refrain_sum/len(avg_years))
    snowmelt_mean = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean = np.array(gl_melt_sum/len(avg_years))
      
    snowfallstd_mean = np.array(snowfallstd_sum/len(avg_years))
    refrainstd_mean = np.array(refrainstd_sum/len(avg_years))
    snowmeltstd_mean = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean = np.array(gl_meltstd_sum/len(avg_years))

    total_accumulation = snowfall_mean + refrain_mean
    total_accumulation_std = snowfallstd_mean + refrainstd_mean
    
    total_ablation = snowmelt_mean + SImelt_mean + gl_melt_mean
    total_ablation_std = snowmeltstd_mean + SImeltstd_mean + gl_meltstd_mean

    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')   
    time = np.arange(0,len(dates))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9

    #ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    a2.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=15)
    a2.plot(time,total_accumulation,c='mediumblue',label='Accumulation')    
    a2.plot(time,-total_ablation,c='red',label='Ablation')    
    a2.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    a2.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    a2.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    a2.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=0,fontsize=14)
    a2.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    a2.grid()
    a2.text(10,0.028,'a)',fontsize=14,weight='bold')
    a2.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    a2.tick_params(axis='y',labelsize=14)
    a2.margins(x=0)
    
    massbal = total_accumulation - total_ablation
    max_mb = (total_accumulation + total_accumulation_std) - (total_ablation - total_ablation_std)
    min_mb = (total_accumulation - total_accumulation_std) - (total_ablation + total_ablation_std)

    transition_date = dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]]
    ax0 = a2.twinx()
    ax0.plot(time,np.cumsum(massbal)*area,c='k',label='Cumulative balance\n$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    print(np.cumsum(massbal)*area)
    #plt.fill_between(time,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='std. dev.')

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass Change (Gt yr$^{-1}$)',fontsize=15)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.legend(fontsize=15,loc='lower left')
        
    handles, labels = a2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    a2.legend(by_label.values(), by_label.keys(),fontsize=15,loc='upper left', ncol=1, borderaxespad=0.19)
    fig.tight_layout()
    #fig.savefig('D:/Model Runs/REF_MODEL/Plots/Refmodel_B0date_MBcurve.pdf',bbox_inches='tight')
    return b0_years[1:],transition_DOYs[1:]
    
def calculate_stddev_timeseries(sim,years,R2S,Glacier_grid,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,x,y,PointScale=True,KWaverage=False,Catchmentaverage=False):
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
    snowmelt_l = []
    refrozenmelt_l = []
    netsnowmelt_l = []
    glaciermelt_l = []
    SImelt_l = []
    rain_l = []
    refrozenrain_l = []
    rainrunoff_l = []
    accumulation_l = []
    snowdepth_l = []
    Ptau_l = []
    SI_l = []
    
    for year in years:
        print('Calculating runoff & mb variables for', year)
        # Get mass balance variable:
        massbal = load_hydrologic_year(sim,year,'Netbalance_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        snowmelt = load_hydrologic_year(sim,year,'Snowmelt_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        refrozenmelt = load_hydrologic_year(sim,year,'Refrozenmelt_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        netsnowmelt = load_hydrologic_year(sim,year,'Netsnowmelt_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        glaciermelt = load_hydrologic_year(sim,year,'Glaciericemelt_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        SImelt = load_hydrologic_year(sim,year,'Superimposedicemelt_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        rain = load_hydrologic_year(sim,year,'Rain_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        refrozenrain = load_hydrologic_year(sim,year,'Refrozenrain_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        rainrunoff = load_hydrologic_year(sim,year,'Rainrunoff_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        accumulation = load_hydrologic_year(sim,year,'Accumulation_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        snowdepth = load_hydrologic_year(sim,year,'Snowdepth_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        Ptau = load_hydrologic_year(sim,year,'Ptau_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        SI = load_hydrologic_year(sim,year,'Superimposedice_std',NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=True)
        
        # Get timeseries from each component at a point location:
        if PointScale == True:
            mean_massbal = massbal[:,x,y]
            mean_snowmelt = snowmelt[:,x,y]
            mean_refrozenmelt = refrozenmelt[:,x,y]
            mean_netsnowmelt = netsnowmelt[:,x,y]
            mean_glaciericemelt = glaciermelt[:,x,y]
            mean_superimposedicemelt = SImelt[:,x,y]
            mean_rain = rain[:,x,y]
            mean_refrozenrain = refrozenrain[:,x,y]
            mean_rainrunoff = rainrunoff[:,x,y]
            mean_accumulation = accumulation[:,x,y]
            mean_snowdepth = snowdepth[:,x,y]
            mean_Ptau = Ptau[:,x,y]
            mean_SI = SI[:,x,y]
            
        
        # Get timeseries from each component for the catchment wide average
        elif Catchmentaverage == True:            
            mean_massbal = np.nanmean(np.nanmean(massbal,axis=1),axis=1)
            mean_snowmelt = np.nanmean(np.nanmean(snowmelt,axis=1),axis=1)
            mean_refrozenmelt = np.nanmean(np.nanmean(refrozenmelt,axis=1),axis=1)
            mean_netsnowmelt = np.nanmean(np.nanmean(netsnowmelt,axis=1),axis=1)
            mean_glaciericemelt = np.nanmean(np.nanmean(glaciermelt,axis=1),axis=1)
            mean_superimposedicemelt = np.nanmean(np.nanmean(SImelt,axis=1),axis=1)
            mean_rain = np.nanmean(np.nanmean(rain,axis=1),axis=1)
            mean_refrozenrain = np.nanmean(np.nanmean(refrozenrain,axis=1),axis=1)
            mean_rainrunoff = np.nanmean(np.nanmean(rainrunoff,axis=1),axis=1)
            mean_accumulation = np.nanmean(np.nanmean(accumulation,axis=1),axis=1)
            mean_snowdepth = np.nanmean(np.nanmean(snowdepth,axis=1),axis=1)
            mean_Ptau = np.nanmean(np.nanmean(Ptau,axis=1),axis=1)
            mean_SI = np.nanmean(np.nanmean(SI,axis=1),axis=1)

        # Get timeseries from each component for the glacier wide average
        elif KWaverage == True:
            for i in range(0,len(snowmelt)):
                massbal[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                snowmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                refrozenmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                netsnowmelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                glaciermelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                SImelt[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                rain[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                refrozenrain[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                rainrunoff[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                accumulation[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                snowdepth[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                Ptau[i][np.where(~np.isfinite(Glacier_grid))] = np.nan
                SI[i][np.where(~np.isfinite(Glacier_grid))] = np.nan

            mean_massbal = np.nanmean(np.nanmean(massbal,axis=1),axis=1)
            mean_snowmelt = np.nanmean(np.nanmean(snowmelt,axis=1),axis=1)
            mean_refrozenmelt = np.nanmean(np.nanmean(refrozenmelt,axis=1),axis=1)
            mean_netsnowmelt = np.nanmean(np.nanmean(netsnowmelt,axis=1),axis=1)
            mean_glaciericemelt = np.nanmean(np.nanmean(glaciermelt,axis=1),axis=1)
            mean_superimposedicemelt = np.nanmean(np.nanmean(SImelt,axis=1),axis=1)
            mean_rain = np.nanmean(np.nanmean(rain,axis=1),axis=1)
            mean_refrozenrain = np.nanmean(np.nanmean(refrozenrain,axis=1),axis=1)
            mean_rainrunoff = np.nanmean(np.nanmean(rainrunoff,axis=1),axis=1)
            mean_accumulation = np.nanmean(np.nanmean(accumulation,axis=1),axis=1)
            mean_snowdepth = np.nanmean(np.nanmean(snowdepth,axis=1),axis=1)
            mean_Ptau = np.nanmean(np.nanmean(Ptau,axis=1),axis=1)
            mean_SI = np.nanmean(np.nanmean(SI,axis=1),axis=1)

        massbal_l.append(mean_massbal)
        snowmelt_l.append(mean_snowmelt)
        refrozenmelt_l.append(mean_refrozenmelt)
        netsnowmelt_l.append(mean_netsnowmelt)
        glaciermelt_l.append(mean_glaciericemelt)
        SImelt_l.append(mean_superimposedicemelt)
        rain_l.append(mean_rain)
        refrozenrain_l.append(mean_refrozenrain)
        rainrunoff_l.append(mean_rainrunoff)
        accumulation_l.append(mean_accumulation)
        snowdepth_l.append(mean_snowdepth)
        Ptau_l.append(mean_Ptau)
        SI_l.append(mean_SI)

    return massbal_l, snowmelt_l, refrozenmelt_l, netsnowmelt_l, glaciermelt_l, SImelt_l, rain_l, refrozenrain_l, rainrunoff_l, accumulation_l, snowdepth_l, Ptau_l, SI_l

def runoff_timeseries_average_discharge_withstddev(units,title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt,gl_icemelt_std, snowmelt_runoff_std, rain_runoff_std, superimp_icemelt_std,area_map):
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt[i][:365])
        snow_sum += np.array(snowmelt_runoff[i][:365])
        rain_sum += np.array(rain_runoff[i][:365])
        SI_sum += np.array(superimp_icemelt[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std[i][:365])
        rainstd_sum += np.array(rain_runoff_std[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean = np.array(ice_sum/len(avg_years))
    snow_mean = np.array(snow_sum/len(avg_years))
    rain_mean = np.array(rain_sum/len(avg_years))
    SI_mean = np.array(SI_sum/len(avg_years))
    total_runoff = ice_mean + snow_mean + rain_mean + SI_mean
    
    icestd_mean = np.array(icestd_sum/len(avg_years))
    snowstd_mean = np.array(snowstd_sum/len(avg_years))
    rainstd_mean = np.array(rainstd_sum/len(avg_years))
    SIstd_mean = np.array(SIstd_sum/len(avg_years))
    total_std = icestd_mean + snowstd_mean + rainstd_mean + SIstd_mean
    
    time = np.arange(0,len(total_runoff))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    ax.plot(time,ice_mean*dc,c='turquoise',label='Glacier ice melt')    
    ax.plot(time,snow_mean*dc,c='royalblue',label='Snow melt')
    ax.plot(time,rain_mean*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean*dc,c='darkorange',label='Superimposed ice melt')
    
    plt.fill_between(time,(ice_mean-icestd_mean)*dc,(ice_mean+icestd_mean)*dc,color='turquoise',alpha=0.35)    
    plt.fill_between(time,(snow_mean-snowstd_mean)*dc,(snow_mean+snowstd_mean)*dc,color='royalblue',alpha=0.35)
    plt.fill_between(time,(rain_mean-rainstd_mean)*dc,(rain_mean+rainstd_mean)*dc,color='deeppink',alpha=0.35)
    plt.fill_between(time,(SI_mean-SIstd_mean)*dc,(SI_mean+SIstd_mean)*dc,color='darkorange',alpha=0.35)
    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(0,daily_runoff_upperlim)
    ax.grid()
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)

    ax0 = ax.twinx()
    ax0.plot(time,np.cumsum(total_runoff)*yc,c='k',label='Cumulative runoff',linewidth=3)
    ax0.plot(time,np.cumsum(ice_mean)*yc,c='turquoise',label='Glacier ice melt',linewidth=3)
    ax0.plot(time,np.cumsum(snow_mean)*yc,c='royalblue',label='Snow melt',linewidth=3)
    ax0.plot(time,np.cumsum(rain_mean)*yc,c='deeppink',label='Rain',linewidth=3)
    ax0.plot(time,np.cumsum(SI_mean)*yc,c='darkorange',label='Superimposed ice melt',linewidth=3)
    ax0.set_ylim(0,cumu_runoff_upperlim)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.margins(x=0)
    
    if units=='mwe':
        ax.set_ylabel('Runoff (m w.e. day$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (m w.e. a$^{-1}$)',fontsize=14)
    elif units=='m3':
        ax.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (km$^3$ a$^{-1}$)',fontsize=14)
    
    plt.fill_between(time,np.cumsum(total_runoff-total_std)*yc,np.cumsum(total_runoff+total_std)*yc,color='k',alpha=0.35)
    plt.fill_between(time,np.cumsum(ice_mean-icestd_mean)*yc,np.cumsum(ice_mean+icestd_mean)*yc,color='turquoise',alpha=0.35)    
    plt.fill_between(time,np.cumsum(snow_mean-snowstd_mean)*yc,np.cumsum(snow_mean+snowstd_mean)*yc,color='royalblue',alpha=0.35)
    plt.fill_between(time,np.cumsum(rain_mean-rainstd_mean)*yc,np.cumsum(rain_mean+rainstd_mean)*yc,color='deeppink',alpha=0.35)
    plt.fill_between(time,np.cumsum(SI_mean-SIstd_mean)*yc,np.cumsum(SI_mean+SIstd_mean)*yc,color='darkorange',alpha=0.35)
    ax0.fill_between(time,-12,-10,color='grey',alpha=0.35,label='Standard deviation')
    
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
    

def runoff_timeseries_average_SLRformat(units,title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt,gl_icemelt_std, snowmelt_runoff_std, rain_runoff_std, superimp_icemelt_std,area_map):
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt[i][:365])
        snow_sum += np.array(snowmelt_runoff[i][:365])
        rain_sum += np.array(rain_runoff[i][:365])
        SI_sum += np.array(superimp_icemelt[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std[i][:365])
        rainstd_sum += np.array(rain_runoff_std[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean = np.array(ice_sum/len(avg_years))
    snow_mean = np.array(snow_sum/len(avg_years))
    rain_mean = np.array(rain_sum/len(avg_years))
    SI_mean = np.array(SI_sum/len(avg_years))
    total_runoff = ice_mean + snow_mean + rain_mean + SI_mean
    
    icestd_mean = np.array(icestd_sum/len(avg_years))
    snowstd_mean = np.array(snowstd_sum/len(avg_years))
    rainstd_mean = np.array(rainstd_sum/len(avg_years))
    SIstd_mean = np.array(SIstd_sum/len(avg_years))
    total_std = icestd_mean + snowstd_mean + rainstd_mean + SIstd_mean
    
    time = np.arange(0,len(total_runoff))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2

    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8.5,4.5))
    #fig.patch.set_facecolor('#f8f5f0')
    #ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    #ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    ax.plot(time,ice_mean*dc,c='turquoise',label='Glacier ice melt')    
    ax.plot(time,snow_mean*dc,c='royalblue',label='Snow melt')
    ax.plot(time,rain_mean*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean*dc,c='darkorange',label='Refrozen ice melt')
    
    #plt.fill_between(time,(total_runoff-total_std)*dc,(total_runoff+total_std)*dc,color='k',alpha=0.35)
    plt.fill_between(time,(ice_mean-icestd_mean)*dc,(ice_mean+icestd_mean)*dc,color='turquoise',alpha=0.35)    
    plt.fill_between(time,(snow_mean-snowstd_mean)*dc,(snow_mean+snowstd_mean)*dc,color='royalblue',alpha=0.35)
    plt.fill_between(time,(rain_mean-rainstd_mean)*dc,(rain_mean+rainstd_mean)*dc,color='deeppink',alpha=0.35)
    plt.fill_between(time,(SI_mean-SIstd_mean)*dc,(SI_mean+SIstd_mean)*dc,color='darkorange',alpha=0.35)
    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335,365])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',''],rotation=45,fontsize=14)
    plt.xlim(0,420)

    ax.set_ylim(0,daily_runoff_upperlim)
    ax.grid()
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)

    ax0 = ax.twinx()
    #ax0.plot(time,np.cumsum(total_runoff)*yc,c='k',label='Cumulative runoff',linewidth=3)
    #ax0.plot(time,np.cumsum(ice_mean)*yc,c='turquoise',label='Glacier ice melt',linewidth=3)
    #ax0.plot(time,np.cumsum(snow_mean)*yc,c='royalblue',label='Snow melt',linewidth=3)
    #ax0.plot(time,np.cumsum(rain_mean)*yc,c='deeppink',label='Rain',linewidth=3)
    #ax0.plot(time,np.cumsum(SI_mean)*yc,c='darkorange',label='Superimposed ice melt',linewidth=3)
    ax0.set_ylim(0,cumu_runoff_upperlim)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.margins(x=0)
    
    if units=='mwe':
        ax.set_ylabel('Runoff (m w.e. day$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (m w.e. a$^{-1}$)',fontsize=14)
    elif units=='m3':
        ax.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (Gt a$^{-1}$)',fontsize=14)
    
    width=10
    plt.vlines(395, ymin=(np.cumsum(total_runoff-total_std)*yc)[-1], ymax=(np.cumsum(total_runoff+total_std)*yc)[-1], linewidth=width, color='k',alpha=0.3)
    plt.hlines(y=(np.cumsum(total_runoff)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='k',linewidth=4,label='Cumulative total runoff')
    
    plt.vlines(395, ymin=(np.cumsum(ice_mean-icestd_mean)*yc)[-1], ymax=(np.cumsum(ice_mean+icestd_mean)*yc)[-1], linewidth=width, color='turquoise',alpha=0.3)
    plt.hlines(y=(np.cumsum(ice_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='turquoise',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(snow_mean-snowstd_mean)*yc)[-1], ymax=(np.cumsum(snow_mean+snowstd_mean)*yc)[-1], linewidth=width, color='royalblue',alpha=0.3)
    plt.hlines(y=(np.cumsum(snow_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='royalblue',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(rain_mean-rainstd_mean)*yc)[-1], ymax=(np.cumsum(rain_mean+rainstd_mean)*yc)[-1], linewidth=width, color='deeppink',alpha=0.3)
    plt.hlines(y=(np.cumsum(rain_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='deeppink',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(SI_mean-SIstd_mean)*yc)[-1], ymax=(np.cumsum(SI_mean+SIstd_mean)*yc)[-1], linewidth=width, color='darkorange',alpha=0.3)
    plt.hlines(y=(np.cumsum(SI_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='darkorange',linewidth=4)
    #plt.fill_between(time,np.cumsum(total_runoff-total_std)*yc,np.cumsum(total_runoff+total_std)*yc,color='k',alpha=0.35)
    #plt.fill_between(time,np.cumsum(ice_mean-icestd_mean)*yc,np.cumsum(ice_mean+icestd_mean)*yc,color='turquoise',alpha=0.35)    
    #plt.fill_between(time,np.cumsum(snow_mean-snowstd_mean)*yc,np.cumsum(snow_mean+snowstd_mean)*yc,color='royalblue',alpha=0.35)
    #plt.fill_between(time,np.cumsum(rain_mean-rainstd_mean)*yc,np.cumsum(rain_mean+rainstd_mean)*yc,color='deeppink',alpha=0.35)
    #plt.fill_between(time,np.cumsum(SI_mean-SIstd_mean)*yc,np.cumsum(SI_mean+SIstd_mean)*yc,color='darkorange',alpha=0.35)
    ax.hlines(y=-10, xmin=375 - width / 2, xmax=375 + width / 2,colors='k',linewidth=4,label='Cumulative total runoff')
    ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='Standard deviation')
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),loc='upper left',fontsize=14, ncol=1, borderaxespad=0.19)
    
    # Create box plots for cumulative values in the second subplot
    ax1 = plt.axes([0.05, 0.12, 0.45, 0.45])
    #box_plot_data = [cum_ice_melt, cum_snow_melt, cum_superimposed_ice_melt, cum_rain]
    #box_colors = ['blue', 'orange', 'green', 'red']
    #box = ax2.boxplot(box_plot_data, positions=[1, 2, 3, 4], widths=0.15, patch_artist=True, showfliers=False)

    piechart_colours = ['turquoise','royalblue','deeppink','darkorange']
    #fig = plt.figure(figsize=(5,5))
    #
    total_runoff_sum = np.sum(ice_mean + snow_mean + SI_mean + rain_mean)
    gl_ice_percent = round((np.sum(ice_mean)/total_runoff_sum)*100,1)
    snow_percent = round((np.sum(snow_mean)/total_runoff_sum)*100,1)
    SIice_percent = round((np.sum(SI_mean)/total_runoff_sum)*100,1)
    rain_percent = round((np.sum(rain_mean)/total_runoff_sum)*100,1)
    
    #ax1.pie([gl_ice_percent,snow_percent,rain_percent,SIice_percent],autopct='%.1f%%',colors=piechart_colours, textprops={'fontsize': 14})
   
    non_zero_percents = [p for p in [gl_ice_percent, snow_percent, rain_percent, SIice_percent] if p != 0]

    if all(p == 0 for p in [gl_ice_percent, snow_percent, rain_percent, SIice_percent]):
        ax1.pie([1], colors=['k'], labels=['No Data'], textprops={'fontsize': 14})
    else:
        #ax1.pie(non_zero_percents, colors=piechart_colours, textprops={'fontsize': 14}, autopct=lambda pct: '' if pct == 0 else f'{pct:.1f}%')
        ax1.pie(non_zero_percents, colors=piechart_colours)

    ax.text(151,102,str(gl_ice_percent) + '%',fontsize=14,weight='bold',color='turquoise')
    ax.text(151,87,str(snow_percent) + '%',fontsize=14,weight='bold',color='royalblue')
    ax.text(151,72,str(rain_percent) + '%',fontsize=14,weight='bold',color='deeppink')
    ax.text(151,57,str(SIice_percent) + '%',fontsize=14,weight='bold',color='darkorange')

    #fig.suptitle(title,fontsize=14,y=1.01) 
    #plt.legend(['Glacier ice melt','Snow melt','Rain','Superimposed ice melt'],fontsize=14,ncol=1,bbox_to_anchor=(1.3,0.8))
    #plt.tight_layout()
    
    fig.tight_layout()
    #plt.savefig('D:/Model Runs/REF_MODEL/Plots/Hydrograph_Catchmentwide_1980-2022_REFMODEL.pdf',bbox_inches='tight')

    print('glacier ice:',str(np.round((np.cumsum(ice_mean)*yc)[-1],2)) + ' Gt a$^{-1}$')
    print('snow:',str(np.round((np.cumsum(snow_mean)*yc)[-1],2)) + ' Gt a$^{-1}$')
    print('rain:',str(np.round((np.cumsum(rain_mean)*yc)[-1],2)) + ' Gt a$^{-1}$')
    print('refrozen ice:',str(np.round((np.cumsum(SI_mean)*yc)[-1],2)) + ' Gt a$^{-1}$')
    print('total runoff:',str(np.round(total_runoff_sum*yc,2)) + ' Gt a$^{-1}$')
    
    print('standard deviations')
    print('total runoff: ',((np.cumsum(total_runoff+total_std)*yc)[-1] - (np.cumsum(total_runoff-total_std)*yc)[-1])/2)
    print('glacier ice runoff: ',((np.cumsum(ice_mean+icestd_mean)*yc)[-1] - (np.cumsum(ice_mean-icestd_mean)*yc)[-1])/2)
    print('snow runoff: ',((np.cumsum(snow_mean+snowstd_mean)*yc)[-1] - (np.cumsum(snow_mean-snowstd_mean)*yc)[-1])/2)
    print('rain runoff: ',((np.cumsum(rain_mean+rainstd_mean)*yc)[-1] -(np.cumsum(rain_mean-rainstd_mean)*yc)[-1])/2)
    print('refrozen ice runoff: ',((np.cumsum(SI_mean+SIstd_mean)*yc)[-1] - (np.cumsum(SI_mean-SIstd_mean)*yc)[-1])/2)

    return ice_mean*dc,snow_mean*dc

def runoff_percentile_curves(units,title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt,gl_icemelt_std, snowmelt_runoff_std, rain_runoff_std, superimp_icemelt_std,area_map):
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    decade_starts = [1980,1990,2000,2010]

    plt.figure()
    for y_i in decade_starts:
        for year in np.arange(y_i,y_i+10):
            i = year - all_years[0]
            ice_sum += np.array(gl_icemelt[i][:365])
            snow_sum += np.array(snowmelt_runoff[i][:365])
            rain_sum += np.array(rain_runoff[i][:365])
            SI_sum += np.array(superimp_icemelt[i][:365])
            
            icestd_sum += np.array(gl_icemelt_std[i][:365])
            snowstd_sum += np.array(snowmelt_runoff_std[i][:365])
            rainstd_sum += np.array(rain_runoff_std[i][:365])
            SIstd_sum += np.array(superimp_icemelt_std[i][:365])
            
        # Get mean daily values and mean daily standard deviations across all years
        ice_mean = np.array(ice_sum/len(avg_years))
        snow_mean = np.array(snow_sum/len(avg_years))
        rain_mean = np.array(rain_sum/len(avg_years))
        SI_mean = np.array(SI_sum/len(avg_years))
        total_runoff = ice_mean + snow_mean + rain_mean + SI_mean
        
        plt.plot(np.cumsum(total_runoff)/np.sum(total_runoff),label=str(y_i) + 's')
    
    plt.xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335],labels=['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],fontsize=14)
    plt.xlim(210,365)
    plt.yticks(fontsize=14)
    plt.ylabel('Fraction of total runoff',fontsize=14)
    plt.legend(fontsize=14)
    plt.grid()
        
    
def gaussian_curve(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(x - mean)**2 / (2 * stddev**2))


def runoff_timeseries_with_gaussian_fit(units,title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt,gl_icemelt_std, snowmelt_runoff_std, rain_runoff_std, superimp_icemelt_std,area_map):
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        #print(year)
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt[i][:365])
        snow_sum += np.array(snowmelt_runoff[i][:365])
        rain_sum += np.array(rain_runoff[i][:365])
        SI_sum += np.array(superimp_icemelt[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std[i][:365])
        rainstd_sum += np.array(rain_runoff_std[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean = np.array(ice_sum/len(avg_years))
    snow_mean = np.array(snow_sum/len(avg_years))
    rain_mean = np.array(rain_sum/len(avg_years))
    SI_mean = np.array(SI_sum/len(avg_years))
    total_runoff = ice_mean + snow_mean + rain_mean + SI_mean
    
    icestd_mean = np.array(icestd_sum/len(avg_years))
    snowstd_mean = np.array(snowstd_sum/len(avg_years))
    rainstd_mean = np.array(rainstd_sum/len(avg_years))
    SIstd_mean = np.array(SIstd_sum/len(avg_years))
    total_std = icestd_mean + snowstd_mean + rainstd_mean + SIstd_mean
    
    time = np.arange(0,len(total_runoff))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8.5,4.5))
    #fig.patch.set_facecolor('#f8f5f0')
    #ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    #ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    ax.plot(time,ice_mean*dc,c='turquoise',label='Glacier ice melt')    
    ax.plot(time,snow_mean*dc,c='royalblue',label='Snow melt')
    ax.plot(time,rain_mean*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean*dc,c='darkorange',label='Refrozen ice melt')
    ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    
    #plt.fill_between(time,(total_runoff-total_std)*dc,(total_runoff+total_std)*dc,color='k',alpha=0.35)
    #plt.fill_between(time,(ice_mean-icestd_mean)*dc,(ice_mean+icestd_mean)*dc,color='turquoise',alpha=0.35)    
    #plt.fill_between(time,(snow_mean-snowstd_mean)*dc,(snow_mean+snowstd_mean)*dc,color='royalblue',alpha=0.35)
    #plt.fill_between(time,(rain_mean-rainstd_mean)*dc,(rain_mean+rainstd_mean)*dc,color='deeppink',alpha=0.35)
    #plt.fill_between(time,(SI_mean-SIstd_mean)*dc,(SI_mean+SIstd_mean)*dc,color='darkorange',alpha=0.35)
    #plt.fill_between(time,(total_runoff-total_std)*dc,(total_runoff+total_std)*dc,color='k',alpha=0.35)
     
    initial_guess = [np.max(total_runoff*dc), np.argmax(total_runoff*dc), 20]
    params, covariance = curve_fit(gaussian_curve, time, total_runoff*dc, p0=initial_guess)
    ax.plot(time, gaussian_curve(time, *params), label='Fitted Gaussian Curve', color='grey', linewidth=3,linestyle='--')

    initial_guess = [np.max(ice_mean*dc), np.argmax(ice_mean*dc), 20]
    params, covariance = curve_fit(gaussian_curve, time, ice_mean*dc, p0=initial_guess)
    ax.plot(time, gaussian_curve(time, *params), color='teal', linewidth=3,linestyle='--')   

    initial_guess = [np.max(snow_mean*dc), np.argmax(snow_mean*dc), 20]
    params, covariance = curve_fit(gaussian_curve, time, snow_mean*dc, p0=initial_guess)
    ax.plot(time, gaussian_curve(time, *params), color='navy', linewidth=3,linestyle='--')   
    
    initial_guess = [np.max(rain_mean*dc), np.argmax(rain_mean*dc), 20]
    params, covariance = curve_fit(gaussian_curve, time, rain_mean*dc, p0=initial_guess)
    ax.plot(time, gaussian_curve(time, *params), color='mediumvioletred', linewidth=3,linestyle='--')   
    
    initial_guess = [np.max(SI_mean*dc), np.argmax(SI_mean*dc), 20]
    params, covariance = curve_fit(gaussian_curve, time, SI_mean*dc, p0=initial_guess)
    ax.plot(time, gaussian_curve(time, *params), color='darkred', linewidth=3,linestyle='--')   
    
    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335,365])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',''],rotation=45,fontsize=14)
    plt.xlim(0,370)
    ax.set_ylim(0,daily_runoff_upperlim)
    ax.grid()
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    
    if units=='mwe':
        ax.set_ylabel('Runoff (m w.e. day$^{-1}$)',fontsize=14)
        #ax0.set_ylabel('Cumulative Runoff (m w.e. a$^{-1}$)',fontsize=14)
    elif units=='m3':
        ax.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=14)
        #ax0.set_ylabel('Cumulative Runoff (Gt$^3$ a$^{-1}$)',fontsize=14)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),loc='upper left',fontsize=14, ncol=1, borderaxespad=0.19)
    
    initial_guess = [np.max(SI_mean*dc), np.argmax(SI_mean*dc), 20]
    params, covariance = curve_fit(gaussian_curve, time, SI_mean*dc, p0=initial_guess)
    ax.plot(time, gaussian_curve(time, *params), label='Fitted Gaussian Curve', color='grey', linewidth=3,linestyle='--')

    return time, (gaussian_curve(time, *params))
  

def runoff_timeseries_with_zerophaseshift(ax,window_size, poly_order,units,title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt,gl_icemelt_std, snowmelt_runoff_std, rain_runoff_std, superimp_icemelt_std,area_map):
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        #print(year)
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt[i][:365])
        snow_sum += np.array(snowmelt_runoff[i][:365])
        rain_sum += np.array(rain_runoff[i][:365])
        SI_sum += np.array(superimp_icemelt[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std[i][:365])
        rainstd_sum += np.array(rain_runoff_std[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean = np.array(ice_sum/len(avg_years))
    snow_mean = np.array(snow_sum/len(avg_years))
    rain_mean = np.array(rain_sum/len(avg_years))
    SI_mean = np.array(SI_sum/len(avg_years))
    total_runoff = ice_mean + snow_mean + rain_mean + SI_mean
    
    icestd_mean = np.array(icestd_sum/len(avg_years))
    snowstd_mean = np.array(snowstd_sum/len(avg_years))
    rainstd_mean = np.array(rainstd_sum/len(avg_years))
    SIstd_mean = np.array(SIstd_sum/len(avg_years))
    total_std = icestd_mean + snowstd_mean + rainstd_mean + SIstd_mean
    
    time = np.arange(0,len(total_runoff))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2
    #fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8.5,4.5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14,weight='bold',loc='left')

    ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    ax.plot(time,ice_mean*dc,c='turquoise',label='Glacier ice runoff')    
    ax.plot(time,snow_mean*dc,c='royalblue',label='Snow runoff')
    ax.plot(time,rain_mean*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean*dc,c='darkorange',label='Refrozen ice runoff')
    #ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    
    ax.fill_between(time,(total_runoff-total_std)*dc,(total_runoff+total_std)*dc,color='grey',alpha=0.35,label='Std. dev.')
    ax.fill_between(time,(ice_mean-icestd_mean)*dc,(ice_mean+icestd_mean)*dc,color='turquoise',alpha=0.35)    
    ax.fill_between(time,(snow_mean-snowstd_mean)*dc,(snow_mean+snowstd_mean)*dc,color='royalblue',alpha=0.35)
    ax.fill_between(time,(rain_mean-rainstd_mean)*dc,(rain_mean+rainstd_mean)*dc,color='deeppink',alpha=0.35)
    ax.fill_between(time,(SI_mean-SIstd_mean)*dc,(SI_mean+SIstd_mean)*dc,color='darkorange',alpha=0.35)
    
    # Implement a Savitzky-Golay filter for smoothing
    runoff_smoothed = savgol_filter(total_runoff*dc, window_size, poly_order)
    ax.plot(time, runoff_smoothed, color='dimgrey', linewidth=3,linestyle='--')

    icemelt_smoothed = savgol_filter(ice_mean*dc, window_size, poly_order)
    ax.plot(time, icemelt_smoothed, color='teal', linewidth=3,linestyle='--')

    snowmelt_smoothed = savgol_filter(snow_mean*dc, window_size, poly_order)
    ax.plot(time, snowmelt_smoothed, color='navy', linewidth=3,linestyle='--')

    rain_smoothed = savgol_filter(rain_mean*dc, window_size, poly_order)
    ax.plot(time, rain_smoothed, color='mediumvioletred', linewidth=3,linestyle='--')

    SI_smoothed = savgol_filter(SI_mean*dc, window_size, poly_order)
    ax.plot(time, SI_smoothed, color='darkred', linewidth=3,linestyle='--')

    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335,365])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',''],rotation=0,fontsize=14)
    plt.xlim(182,370)
    ax.set_ylim(0,daily_runoff_upperlim)
    ax.grid()
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)

    total_runoff_sum = np.sum(ice_mean + snow_mean + SI_mean + rain_mean)
    gl_ice_percent = round((np.sum(ice_mean)/total_runoff_sum)*100,1)
    snow_percent = round((np.sum(snow_mean)/total_runoff_sum)*100,1)
    SIice_percent = round((np.sum(SI_mean)/total_runoff_sum)*100,1)
    rain_percent = round((np.sum(rain_mean)/total_runoff_sum)*100,1)
    
    ax.text(185,140,str(gl_ice_percent) + '%',fontsize=14,weight='bold',color='turquoise')
    ax.text(185,105,str(snow_percent) + '%',fontsize=14,weight='bold',color='royalblue')
    ax.text(191,70,str(rain_percent) + '%',fontsize=14,weight='bold',color='deeppink')
    ax.text(191,35,str(SIice_percent) + '%',fontsize=14,weight='bold',color='darkorange')


    dates = pd.date_range(start= '2000-10-01 00:00:00',end= '2001-09-30 21:00:00',freq='1D')
    Mice_dominates = time[212:][np.where(icemelt_smoothed[212:]>snowmelt_smoothed[212:])][0]
    print(Mice_dominates,str(dates[Mice_dominates])[5:10])
    ax.text(363,381,'M$_{gl. ice}$ > M$_{snow}$\non ' + str(dates[Mice_dominates])[5:10],fontsize=14,color='darkslategrey',horizontalalignment='right')
    
    runoff_smoothed = savgol_filter(total_runoff*dc, window_size, poly_order)
    return time, runoff_smoothed,total_runoff*dc

def plot_pie_chart(ax, sizes, labels):
    piechart_colours = ['turquoise','royalblue','deeppink','darkorange']
    ax.pie(sizes, colors=piechart_colours)
    
def decadal_hydrographs_smoothed_curves():
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2,figsize=(10,6))
    runoff_timeseries_with_zerophaseshift(ax1,51,3,'m3','a) ',np.arange(1980,1990),years,470,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
    ax1.set_xlim(182,365)
    ax1.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=14)
    axins1 = inset_axes(ax1,width="60%", height="60%", bbox_to_anchor=(-0.52, 0.038, 1, 1), bbox_transform=ax1.transAxes)
    
    runoff_timeseries_with_zerophaseshift(ax2,51,3,'m3','b) ',np.arange(1990,2000),years,470,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
    ax2.set_xlim(182,365)
    axins2 = inset_axes(ax2,width="60%", height="60%", bbox_to_anchor=(-0.52, 0.038, 1, 1), bbox_transform=ax2.transAxes)
    
    runoff_timeseries_with_zerophaseshift(ax3,51,3,'m3','c) ',np.arange(2000,2010),years,470,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
    ax3.set_xlim(182,365)
    ax3.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=14)
    axins3 = inset_axes(ax3,width="60%", height="60%", bbox_to_anchor=(-0.52, 0.038, 1, 1), bbox_transform=ax3.transAxes)
    
    runoff_timeseries_with_zerophaseshift(ax4,51,3,'m3','d) ',np.arange(2010,2020),years,470,2.75,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc)
    ax4.set_xlim(182,365)
    axins4 = inset_axes(ax4,width="60%", height="60%", bbox_to_anchor=(-0.52, 0.038, 1, 1), bbox_transform=ax4.transAxes)
    
    handles, labels = ax4.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(0.98,1.1),fontsize=14, ncol=3, borderaxespad=0.19)
     
    plot_pie_chart(axins1, sizes=[59.3,33.0,5.5,2.3], labels=['A', 'B', 'C', 'D'])
    plot_pie_chart(axins2, sizes=[60.5,31.7,5.5,2.3], labels=['A', 'B', 'C', 'D'])
    plot_pie_chart(axins3, sizes=[62.0,30.2,5.8,2.0], labels=['A', 'B', 'C', 'D'])
    plot_pie_chart(axins4, sizes=[62.7,28.5,6.8,2.0], labels=['A', 'B', 'C', 'D'])
    
    fig.tight_layout()
    #fig.savefig('D:/Model Runs/REF_MODEL/Plots/Refmodel_decadal_hydrographs_smoothed.pdf',bbox_inches='tight')


# =============================================================================
# FUNCTIONS TO PLOT THE DIFFERENCES BETWEEN TWO MODELS 
# =============================================================================
    
def distributed_mass_balance_difference(model_name,avg_years,all_years,massbal_dist_ref,massbal_dist_alt,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY_REF = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_REF[year-avg_years[0],:,:] = massbal_dist_ref[year-all_years[0]]
    MB_REF = np.nanmean(MB_ARRAY_REF,axis=0)
        
    MB_ARRAY_ALT = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_ALT[year-avg_years[0],:,:] = massbal_dist_alt[year-all_years[0]]
    MB_ALT = np.nanmean(MB_ARRAY_ALT,axis=0)
    
    DIFF = np.array(MB_REF-MB_ALT)
    DIFF[np.where(DIFF==0)] = np.nan
        
    print(np.nanmin(DIFF),np.nanmax(DIFF))
    
    plt.figure(figsize=(9,5))
    plt.title('Mass Balance ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='massbal_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,20,1))
    legend.ax.set_ylabel('Difference (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=1)
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    #plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nmass balance\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    plt.tight_layout()
    
def distributed_runoff_difference(model_name,avg_years,all_years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,netsnowmelt_dist_alt,gl_icemelt_dist_alt,superimp_icemelt_dist_alt,rain_runoff_dist_alt,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY_REF = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_REF[year-avg_years[0],:,:] = netsnowmelt_dist_ref[year-all_years[0]] + gl_icemelt_dist_ref[year-all_years[0]] + superimp_icemelt_dist_ref[year-all_years[0]] + rain_runoff_dist_ref[year-all_years[0]]
    MB_REF = np.nanmean(MB_ARRAY_REF,axis=0)
        
    MB_ARRAY_ALT = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_ALT[year-avg_years[0],:,:] = netsnowmelt_dist_alt[year-all_years[0]] + gl_icemelt_dist_alt[year-all_years[0]] + superimp_icemelt_dist_alt[year-all_years[0]] + rain_runoff_dist_alt[year-all_years[0]]
    MB_ALT = np.nanmean(MB_ARRAY_ALT,axis=0)
    
    DIFF = np.array(MB_REF-MB_ALT)
    DIFF[np.where(DIFF==0)] = np.nan
        
    print(np.nanmin(DIFF),np.nanmax(DIFF))    
    
    fig = plt.figure(figsize=(9,5))
    plt.title('Total Runoff ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='runoff_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,20,1))
    legend.ax.set_ylabel('Difference (m w.e. a$^{-1}$)', rotation=270,fontsize=16,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=16)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    fig.patch.set_facecolor('#dce3e4')
    plt.axis('off')
    plt.gca().set_facecolor('#dce3e4')
    plt.tight_layout()
    
def distributed_glaciermelt_difference(model_name,avg_years,all_years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,netsnowmelt_dist_alt,gl_icemelt_dist_alt,superimp_icemelt_dist_alt,rain_runoff_dist_alt,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY_REF = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_REF[year-avg_years[0],:,:] = gl_icemelt_dist_ref[year-all_years[0]] 
    MB_REF = np.nanmean(MB_ARRAY_REF,axis=0)
        
    MB_ARRAY_ALT = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_ALT[year-avg_years[0],:,:] = gl_icemelt_dist_alt[year-all_years[0]] 
    MB_ALT = np.nanmean(MB_ARRAY_ALT,axis=0)
    
    DIFF = np.array(MB_REF-MB_ALT)
    DIFF[np.where(DIFF==0)] = np.nan
        
    print(np.nanmin(DIFF),np.nanmax(DIFF))    
    
    plt.figure(figsize=(9,5))
    plt.title('Total glacier ice melt ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='runoff_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,20,0.5))
    legend.ax.set_ylabel('Difference (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.tight_layout()
    
def distributed_snowrunoff_difference(model_name,avg_years,all_years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,netsnowmelt_dist_alt,gl_icemelt_dist_alt,superimp_icemelt_dist_alt,rain_runoff_dist_alt,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY_REF = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_REF[year-avg_years[0],:,:] = netsnowmelt_dist_ref[year-all_years[0]]
    MB_REF = np.nanmean(MB_ARRAY_REF,axis=0)
        
    MB_ARRAY_ALT = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_ALT[year-avg_years[0],:,:] = netsnowmelt_dist_alt[year-all_years[0]] 
    MB_ALT = np.nanmean(MB_ARRAY_ALT,axis=0)
    
    DIFF = np.array(MB_REF-MB_ALT)
    DIFF[np.where(DIFF==0)] = np.nan
        
    print(np.nanmin(DIFF),np.nanmax(DIFF))    
    
    plt.figure(figsize=(9,5))
    plt.title('Runoff from snow melt ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='snowmelt_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-1,1,0.1))
    legend.ax.set_ylabel('Difference (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.tight_layout()
    
def distributed_rainrunoff_difference(model_name,avg_years,all_years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,netsnowmelt_dist_alt,gl_icemelt_dist_alt,superimp_icemelt_dist_alt,rain_runoff_dist_alt,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY_REF = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_REF[year-avg_years[0],:,:] =  rain_runoff_dist_ref[year-all_years[0]]
    MB_REF = np.nanmean(MB_ARRAY_REF,axis=0)
        
    MB_ARRAY_ALT = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_ALT[year-avg_years[0],:,:] =  rain_runoff_dist_alt[year-all_years[0]]
    MB_ALT = np.nanmean(MB_ARRAY_ALT,axis=0)
    
    DIFF = np.array(MB_REF-MB_ALT)
    DIFF[np.where(DIFF==0)] = np.nan
        
    print(np.nanmin(DIFF),np.nanmax(DIFF))    
    
    plt.figure(figsize=(9,5))
    plt.title('Runoff from rainfall ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='rainrunoff_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-0.06,0.004,0.01))
    legend.ax.set_ylabel('Difference (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.tight_layout()
    
def distributed_SImelt_difference(model_name,avg_years,all_years,netsnowmelt_dist_ref,gl_icemelt_dist_ref,superimp_icemelt_dist_ref,rain_runoff_dist_ref,netsnowmelt_dist_alt,gl_icemelt_dist_alt,superimp_icemelt_dist_alt,rain_runoff_dist_alt,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY_REF = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_REF[year-avg_years[0],:,:] =  superimp_icemelt_dist_ref[year-all_years[0]]
    MB_REF = np.nanmean(MB_ARRAY_REF,axis=0)
        
    MB_ARRAY_ALT = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        MB_ARRAY_ALT[year-avg_years[0],:,:] = superimp_icemelt_dist_alt[year-all_years[0]]
    MB_ALT = np.nanmean(MB_ARRAY_ALT,axis=0)
    
    DIFF = np.array(MB_REF-MB_ALT)
    DIFF[np.where(DIFF==0)] = np.nan
        
    print(np.nanmin(DIFF),np.nanmax(DIFF))    
    
    plt.figure(figsize=(9,5))
    plt.title('Superimposed ice melt ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='SImelt_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-1,1,0.01))
    legend.ax.set_ylabel('Difference (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.tight_layout()




def compare_hydrographs(units,title,alt_model_name,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt_ref, snowmelt_runoff_ref, rain_runoff_ref, superimp_icemelt_ref,gl_icemelt_std_ref, snowmelt_runoff_std_ref, rain_runoff_std_ref, superimp_icemelt_std_ref, \
                                        gl_icemelt_alt, snowmelt_runoff_alt, rain_runoff_alt, superimp_icemelt_alt,gl_icemelt_std_alt, snowmelt_runoff_std_alt, rain_runoff_std_alt, superimp_icemelt_std_alt,area_map):
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt_ref[i][:365])
        snow_sum += np.array(snowmelt_runoff_ref[i][:365])
        rain_sum += np.array(rain_runoff_ref[i][:365])
        SI_sum += np.array(superimp_icemelt_ref[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std_ref[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std_ref[i][:365])
        rainstd_sum += np.array(rain_runoff_std_ref[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std_ref[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean_ref = np.array(ice_sum/len(avg_years))
    snow_mean_ref = np.array(snow_sum/len(avg_years))
    rain_mean_ref = np.array(rain_sum/len(avg_years))
    SI_mean_ref = np.array(SI_sum/len(avg_years))
    total_runoff_ref = ice_mean_ref + snow_mean_ref + rain_mean_ref + SI_mean_ref
    
    icestd_mean_ref = np.array(icestd_sum/len(avg_years))
    snowstd_mean_ref = np.array(snowstd_sum/len(avg_years))
    rainstd_mean_ref = np.array(rainstd_sum/len(avg_years))
    SIstd_mean_ref = np.array(SIstd_sum/len(avg_years))
    total_std_ref = icestd_mean_ref + snowstd_mean_ref + rainstd_mean_ref + SIstd_mean_ref
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt_alt[i][:365])
        snow_sum += np.array(snowmelt_runoff_alt[i][:365])
        rain_sum += np.array(rain_runoff_alt[i][:365])
        SI_sum += np.array(superimp_icemelt_alt[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std_alt[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std_alt[i][:365])
        rainstd_sum += np.array(rain_runoff_std_alt[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std_alt[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean_alt = np.array(ice_sum/len(avg_years))
    snow_mean_alt = np.array(snow_sum/len(avg_years))
    rain_mean_alt = np.array(rain_sum/len(avg_years))
    SI_mean_alt = np.array(SI_sum/len(avg_years))
    total_runoff_alt = ice_mean_alt + snow_mean_alt + rain_mean_alt + SI_mean_alt
    
    icestd_mean_alt = np.array(icestd_sum/len(avg_years))
    snowstd_mean_alt = np.array(snowstd_sum/len(avg_years))
    rainstd_mean_alt = np.array(rainstd_sum/len(avg_years))
    SIstd_mean_alt = np.array(SIstd_sum/len(avg_years))
    total_std_alt = icestd_mean_alt + snowstd_mean_alt + rainstd_mean_alt + SIstd_mean_alt
    
    time = np.arange(0,len(total_runoff_ref))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2
   
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(9,5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    #ax.plot(time,total_runoff_ref*dc,c='k',label='Total runoff')
    ax.plot(time,ice_mean_ref*dc,c='turquoise',label='Glacier ice melt')    
    ax.plot(time,snow_mean_ref*dc,c='royalblue',label='Snow melt')
    ax.plot(time,rain_mean_ref*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean_ref*dc,c='darkorange',label='Superimposed ice melt')
    
    ax.plot(time,ice_mean_alt*dc,c='#2d9d92',linestyle='--')    
    ax.plot(time,snow_mean_alt*dc,c='#2d4a9d',linestyle='--') 
    ax.plot(time,rain_mean_alt*dc,c='#b20e67',linestyle='--') 
    ax.plot(time,SI_mean_alt*dc,c='#b26200',linestyle='--') 
    
    plt.fill_between(time,(ice_mean_ref-icestd_mean_ref)*dc,(ice_mean_ref+icestd_mean_ref)*dc,color='turquoise',alpha=0.35)    
    plt.fill_between(time,(snow_mean_ref-snowstd_mean_ref)*dc,(snow_mean_ref+snowstd_mean_ref)*dc,color='royalblue',alpha=0.35)
    plt.fill_between(time,(rain_mean_ref-rainstd_mean_ref)*dc,(rain_mean_ref+rainstd_mean_ref)*dc,color='deeppink',alpha=0.35)
    plt.fill_between(time,(SI_mean_ref-SIstd_mean_ref)*dc,(SI_mean_ref+SIstd_mean_ref)*dc,color='darkorange',alpha=0.35)
    
    plt.fill_between(time,(ice_mean_alt-icestd_mean_alt)*dc,(ice_mean_alt+icestd_mean_alt)*dc,color='#2d9d92',alpha=0.35)    
    plt.fill_between(time,(snow_mean_alt-snowstd_mean_alt)*dc,(snow_mean_alt+snowstd_mean_alt)*dc,color='#2d4a9d',alpha=0.35)
    plt.fill_between(time,(rain_mean_alt-rainstd_mean_alt)*dc,(rain_mean_alt+rainstd_mean_alt)*dc,color='#b20e67',alpha=0.35)
    plt.fill_between(time,(SI_mean_alt-SIstd_mean_alt)*dc,(SI_mean_alt+SIstd_mean_alt)*dc,color='#b26200',alpha=0.35)
    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335,365,395,425])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','REF_MODEL',alt_model_name],rotation=45,fontsize=14)
    ax.grid()
    plt.xlim(0,455)
    ax.set_ylim(0,daily_runoff_upperlim)

    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)

    ax0 = ax.twinx()
    ax0.set_ylim(0,cumu_runoff_upperlim)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.margins(x=0)
    #ax0.set_xticks(ticks=[395,425])
    #ax0.set_xticklabels(['REF_MODEL','ROUNCE_DEBRIS'],rotation=45,fontsize=14)
    
    if units=='mwe':
        ax.set_ylabel('Runoff (m w.e. day$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (m w.e. a$^{-1}$)',fontsize=14)
    elif units=='m3':
        ax.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (km$^3$ a$^{-1}$)',fontsize=14)
    
    width=13
    plt.vlines(395, ymin=(np.cumsum(total_runoff_ref-total_std_ref)*yc)[-1], ymax=(np.cumsum(total_runoff_ref+total_std_ref)*yc)[-1], linewidth=width, color='k',alpha=0.3)
    plt.hlines(y=(np.cumsum(total_runoff_ref)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='k',linewidth=4,label='Cumulative total runoff')
    
    plt.vlines(395, ymin=(np.cumsum(ice_mean_ref-icestd_mean_ref)*yc)[-1], ymax=(np.cumsum(ice_mean_ref+icestd_mean_ref)*yc)[-1], linewidth=width, color='turquoise',alpha=0.3)
    plt.hlines(y=(np.cumsum(ice_mean_ref)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='turquoise',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(snow_mean_ref-snowstd_mean_ref)*yc)[-1], ymax=(np.cumsum(snow_mean_ref+snowstd_mean_ref)*yc)[-1], linewidth=width, color='royalblue',alpha=0.3)
    plt.hlines(y=(np.cumsum(snow_mean_ref)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='royalblue',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(rain_mean_ref-rainstd_mean_ref)*yc)[-1], ymax=(np.cumsum(rain_mean_ref+rainstd_mean_ref)*yc)[-1], linewidth=width, color='deeppink',alpha=0.3)
    plt.hlines(y=(np.cumsum(rain_mean_ref)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='deeppink',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(SI_mean_ref-SIstd_mean_ref)*yc)[-1], ymax=(np.cumsum(SI_mean_ref+SIstd_mean_ref)*yc)[-1], linewidth=width, color='darkorange',alpha=0.3)
    plt.hlines(y=(np.cumsum(SI_mean_ref)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='darkorange',linewidth=4)

    plt.vlines(425, ymin=(np.cumsum(total_runoff_alt-total_std_alt)*yc)[-1], ymax=(np.cumsum(total_runoff_alt+total_std_alt)*yc)[-1], linewidth=width, color='k',alpha=0.3)
    plt.hlines(y=(np.cumsum(total_runoff_alt)*yc)[-1], xmin=425 - width / 2, xmax=425 + width / 2,colors='k',linewidth=4,label='Cumulative total runoff')
    
    plt.vlines(425, ymin=(np.cumsum(ice_mean_alt-icestd_mean_alt)*yc)[-1], ymax=(np.cumsum(ice_mean_alt+icestd_mean_alt)*yc)[-1], linewidth=width, color='#2d9d92',alpha=0.3)
    plt.hlines(y=(np.cumsum(ice_mean_alt)*yc)[-1], xmin=425 - width / 2, xmax=425 + width / 2,colors='#2d9d92',linewidth=4)
    
    plt.vlines(425, ymin=(np.cumsum(snow_mean_alt-snowstd_mean_alt)*yc)[-1], ymax=(np.cumsum(snow_mean_alt+snowstd_mean_alt)*yc)[-1], linewidth=width, color='#2d4a9d',alpha=0.3)
    plt.hlines(y=(np.cumsum(snow_mean_alt)*yc)[-1], xmin=425 - width / 2, xmax=425 + width / 2,colors='#2d4a9d',linewidth=4)
    
    plt.vlines(425, ymin=(np.cumsum(rain_mean_alt-rainstd_mean_alt)*yc)[-1], ymax=(np.cumsum(rain_mean_alt+rainstd_mean_alt)*yc)[-1], linewidth=width, color='#b20e67',alpha=0.3)
    plt.hlines(y=(np.cumsum(rain_mean_alt)*yc)[-1], xmin=425 - width / 2, xmax=425 + width / 2,colors='#b20e67',linewidth=4)
    
    plt.vlines(425, ymin=(np.cumsum(SI_mean_alt-SIstd_mean_alt)*yc)[-1], ymax=(np.cumsum(SI_mean_alt+SIstd_mean_alt)*yc)[-1], linewidth=width, color='#b26200',alpha=0.3)
    plt.hlines(y=(np.cumsum(SI_mean_alt)*yc)[-1], xmin=425 - width / 2, xmax=425 + width / 2,colors='#b26200',linewidth=4)

    ax.hlines(y=-10, xmin=375 - width / 2, xmax=375 + width / 2,colors='k',linewidth=4,label='Cumulative total runoff')
    ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='Standard deviation')
    
    #'REF_MODEL','ROUNCE_DEBRIS'
    #ax.text(395,250,s='REF\nMODEL',horizontalalignment='center',verticalalignment='center')
    
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),loc='upper left',fontsize=14, ncol=1, borderaxespad=0.19)
    
def massbalance_timeseries_comparison(title,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim,  \
                                      accumulation_ref, refrozen_rain_ref, netsnowmelt_ref, superimp_icemelt_ref, gl_icemelt_ref, accumulation_std_ref, refrozen_rain_std_ref, netsnowmelt_std_ref, superimp_icemelt_std_ref, gl_icemelt_std_ref, \
                                      accumulation_alt, refrozen_rain_alt, netsnowmelt_alt, superimp_icemelt_alt, gl_icemelt_alt, accumulation_std_alt, refrozen_rain_std_alt, netsnowmelt_std_alt, superimp_icemelt_std_alt, gl_icemelt_std_alt):
    
    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation_alt[i][:365])
        refrain_sum += np.array(refrozen_rain_alt[i][:365])
        snowmelt_sum += np.array(netsnowmelt_alt[i][:365])
        SImelt_sum += np.array(superimp_icemelt_alt[i][:365])
        gl_melt_sum += np.array(gl_icemelt_alt[i][:365]) 
        
        snowfallstd_sum += np.array(accumulation_std_alt[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std_alt[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std_alt[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std_alt[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std_alt[i][:365]) 
        
    snowfall_mean_alt = np.array(snowfall_sum/len(avg_years))
    refrain_mean_alt = np.array(refrain_sum/len(avg_years))
    snowmelt_mean_alt = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean_alt = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean_alt = np.array(gl_melt_sum/len(avg_years))
      
    snowfallstd_mean_alt = np.array(snowfallstd_sum/len(avg_years))
    refrainstd_mean_alt = np.array(refrainstd_sum/len(avg_years))
    snowmeltstd_mean_alt = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean_alt = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean_alt = np.array(gl_meltstd_sum/len(avg_years))

    total_accumulation_alt = snowfall_mean_alt + refrain_mean_alt
    total_accumulation_std_alt = snowfallstd_mean_alt + refrainstd_mean_alt
    
    total_ablation_alt = snowmelt_mean_alt + SImelt_mean_alt + gl_melt_mean_alt
    total_ablation_std_alt = snowmeltstd_mean_alt + SImeltstd_mean_alt + gl_meltstd_mean_alt
    

    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation_ref[i][:365])
        refrain_sum += np.array(refrozen_rain_ref[i][:365])
        snowmelt_sum += np.array(netsnowmelt_ref[i][:365])
        SImelt_sum += np.array(superimp_icemelt_ref[i][:365])
        gl_melt_sum += np.array(gl_icemelt_ref[i][:365]) 
        
        snowfallstd_sum += np.array(accumulation_std_ref[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std_ref[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std_ref[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std_ref[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std_ref[i][:365]) 
        
    snowfall_mean_ref = np.array(snowfall_sum/len(avg_years))
    refrain_mean_ref = np.array(refrain_sum/len(avg_years))
    snowmelt_mean_ref = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean_ref = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean_ref = np.array(gl_melt_sum/len(avg_years))
      
    snowfallstd_mean_ref = np.array(snowfallstd_sum/len(avg_years))
    refrainstd_mean_ref = np.array(refrainstd_sum/len(avg_years))
    snowmeltstd_mean_ref = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean_ref = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean_ref = np.array(gl_meltstd_sum/len(avg_years))

    total_accumulation_ref = snowfall_mean_ref + refrain_mean_ref
    total_accumulation_std_ref = snowfallstd_mean_ref + refrainstd_mean_ref
    
    total_ablation_ref = snowmelt_mean_ref + SImelt_mean_ref + gl_melt_mean_ref
    total_ablation_std_ref = snowmeltstd_mean_ref + SImeltstd_mean_ref + gl_meltstd_mean_ref


    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')   
    time = np.arange(0,len(dates))

    # Plot the average
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=14)
    ax.plot(time,total_accumulation_alt,c='mediumblue',label='Accumulation',linestyle='--')    
    ax.plot(time,-total_ablation_alt,c='red',label='Ablation',linestyle='--')    
    plt.fill_between(time,(total_accumulation_alt-total_accumulation_std_alt),(total_accumulation_alt+total_accumulation_std_alt),color='mediumblue',alpha=0.35)    
    plt.fill_between(time,(-total_ablation_alt-total_ablation_std_alt),(-total_ablation_alt+total_ablation_std_alt),color='red',alpha=0.35)
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax.grid()
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='grey')
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    
    #ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    #ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=14)
    ax.plot(time,total_accumulation_ref,c='mediumblue',label='Accumulation')    
    ax.plot(time,-total_ablation_ref,c='red',label='Ablation')    
    plt.fill_between(time,(total_accumulation_ref-total_accumulation_std_ref),(total_accumulation_ref+total_accumulation_std_ref),color='mediumblue',alpha=0.35)    
    plt.fill_between(time,(-total_ablation_ref-total_ablation_std_ref),(-total_ablation_ref+total_ablation_std_ref),color='red',alpha=0.35)
    #ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    #ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    #ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    #ax.grid()
    #ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    #ax.tick_params(axis='y',labelsize=14)
    #ax.margins(x=0)
    
    massbal_alt = total_accumulation_alt - total_ablation_alt
    max_mb_alt = (total_accumulation_alt + total_accumulation_std_alt) - (total_ablation_alt - total_ablation_std_alt)
    min_mb_alt = (total_accumulation_alt - total_accumulation_std_alt) - (total_ablation_alt + total_ablation_std_alt)

    massbal_ref = total_accumulation_ref - total_ablation_ref

    #transition_date = dates[50:][np.where(np.cumsum(massbal_alt)[50:] <=0)[0][0]]
    ax0 = ax.twinx()
    ax0.plot(time,np.cumsum(massbal_alt),c='k',linewidth=3,linestyle='--')
    ax0.plot(time,np.cumsum(massbal_ref),c='k',linewidth=3)
    #plt.fill_between(time,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='std. dev.')

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass Balance (m w.e.)',fontsize=14)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    #ax0.legend(fontsize=14,loc='lower left')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
    fig.tight_layout()
    
    
def compare_date_of_zero_balance(years,mb_ref,mb_alt):
    # transition date
    transition_dates_alt = []
    transition_DOYs_alt = []
    for year in years:
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
        massbal_alt = mb_alt[year-years[0]]
        
        d_tr_alt = [np.where(np.cumsum(massbal_alt)[50:] <=0)][0][0]
        if len(d_tr_alt) == 0:
            transition_date = 'N/A'
            transition_DOY = np.nan
        else:
            transition_date = (dates[50:][np.where(np.cumsum(massbal_alt)[50:] <=0)[0][0]])
            transition_DOY = np.where(np.cumsum(massbal_alt)[50:] <=0)[0][0]
        
        transition_dates_alt.append(transition_date)
        transition_DOYs_alt.append(transition_DOY+50)
        
    monthday = []
    for i in np.arange(min(transition_DOYs_alt),max(transition_DOYs_alt)+1,10):
        monthday.append( str(dates[i])[5:10])
        
        
    m, b = np.polyfit(years[np.isfinite(transition_DOYs_alt)],np.array(transition_DOYs_alt)[np.isfinite(transition_DOYs_alt)],deg=1)
    
    plt.figure(figsize=(12,5))
    plt.title('Date of zero balance',fontsize=14)
    plt.scatter(years,transition_DOYs_alt,c='darkblue',s=90,label='ROUNCE_DEBRIS')
    #plt.yticks(np.arange(min(transition_DOYs_alt),max(transition_DOYs_alt)+1,10),monthday,fontsize=14)
    #plt.xticks(fontsize=14)
    #plt.ylabel('date of $\dot{B}$ = 0',fontsize=14)
    #plt.grid()
    plt.plot(years,m*years + b,linestyle='--',c='darkblue',linewidth=2)
    
    
    # transition date
    transition_dates_ref = []
    transition_DOYs_ref = []
    for year in years:
        dates = pd.date_range(start= str(year) + '-10-01 00:00:00',end= str(year+1) + '-09-30 21:00:00',freq='1D')
        massbal_ref = mb_ref[year-years[0]]
        
        d_tr_ref = [np.where(np.cumsum(massbal_ref)[50:] <=0)][0][0]
        if len(d_tr_ref) == 0:
            transition_date = 'N/A'
            transition_DOY = np.nan
        else:
            transition_date = (dates[50:][np.where(np.cumsum(massbal_ref)[50:] <=0)[0][0]])
            transition_DOY = np.where(np.cumsum(massbal_ref)[50:] <=0)[0][0]
        
        transition_dates_ref.append(transition_date)
        transition_DOYs_ref.append(transition_DOY+50)
        
    monthday = []
    for i in np.arange(min(transition_DOYs_ref),max(transition_DOYs_ref)+1,10):
        monthday.append( str(dates[i])[5:10])
        
        
    m, b = np.polyfit(years[np.isfinite(transition_DOYs_ref)],np.array(transition_DOYs_ref)[np.isfinite(transition_DOYs_ref)],deg=1)
    
    plt.scatter(years,transition_DOYs_ref,c='red',s=60,label='REF_MODEL')
    plt.yticks(np.arange(min(transition_DOYs_alt),max(transition_DOYs_alt)+1,10),monthday,fontsize=14)
    plt.xticks(fontsize=14)
    plt.ylabel('date of $\dot{B}$ = 0',fontsize=14)
    plt.grid()
    plt.plot(years,m*years + b,linestyle='--',c='red',linewidth=2)
    plt.legend(fontsize=14)
    
def massbalance_timeseries_difference(title,alt_model_name,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim,  \
                                      accumulation_ref, refrozen_rain_ref, netsnowmelt_ref, superimp_icemelt_ref, gl_icemelt_ref, accumulation_std_ref, refrozen_rain_std_ref, netsnowmelt_std_ref, superimp_icemelt_std_ref, gl_icemelt_std_ref, \
                                      accumulation_alt, refrozen_rain_alt, netsnowmelt_alt, superimp_icemelt_alt, gl_icemelt_alt, accumulation_std_alt, refrozen_rain_std_alt, netsnowmelt_std_alt, superimp_icemelt_std_alt, gl_icemelt_std_alt):
    
    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation_alt[i][:365])
        refrain_sum += np.array(refrozen_rain_alt[i][:365])
        snowmelt_sum += np.array(netsnowmelt_alt[i][:365])
        SImelt_sum += np.array(superimp_icemelt_alt[i][:365])
        gl_melt_sum += np.array(gl_icemelt_alt[i][:365]) 
        
        snowfallstd_sum += np.array(accumulation_std_alt[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std_alt[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std_alt[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std_alt[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std_alt[i][:365]) 
        
    snowfall_mean_alt = np.array(snowfall_sum/len(avg_years))
    refrain_mean_alt = np.array(refrain_sum/len(avg_years))
    snowmelt_mean_alt = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean_alt = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean_alt = np.array(gl_melt_sum/len(avg_years))
      
    snowfallstd_mean_alt = np.array(snowfallstd_sum/len(avg_years))
    refrainstd_mean_alt = np.array(refrainstd_sum/len(avg_years))
    snowmeltstd_mean_alt = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean_alt = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean_alt = np.array(gl_meltstd_sum/len(avg_years))

    total_accumulation_alt = snowfall_mean_alt + refrain_mean_alt
    total_accumulation_std_alt = snowfallstd_mean_alt + refrainstd_mean_alt
    
    total_ablation_alt = snowmelt_mean_alt + SImelt_mean_alt + gl_melt_mean_alt
    total_ablation_std_alt = snowmeltstd_mean_alt + SImeltstd_mean_alt + gl_meltstd_mean_alt
    

    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation_ref[i][:365])
        refrain_sum += np.array(refrozen_rain_ref[i][:365])
        snowmelt_sum += np.array(netsnowmelt_ref[i][:365])
        SImelt_sum += np.array(superimp_icemelt_ref[i][:365])
        gl_melt_sum += np.array(gl_icemelt_ref[i][:365]) 
        
        snowfallstd_sum += np.array(accumulation_std_ref[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std_ref[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std_ref[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std_ref[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std_ref[i][:365]) 
        
    snowfall_mean_ref = np.array(snowfall_sum/len(avg_years))
    refrain_mean_ref = np.array(refrain_sum/len(avg_years))
    snowmelt_mean_ref = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean_ref = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean_ref = np.array(gl_melt_sum/len(avg_years))
      
    snowfallstd_mean_ref = np.array(snowfallstd_sum/len(avg_years))
    refrainstd_mean_ref = np.array(refrainstd_sum/len(avg_years))
    snowmeltstd_mean_ref = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean_ref = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean_ref = np.array(gl_meltstd_sum/len(avg_years))

    total_accumulation_ref = snowfall_mean_ref + refrain_mean_ref
    total_accumulation_std_ref = snowfallstd_mean_ref + refrainstd_mean_ref
    
    total_ablation_ref = snowmelt_mean_ref + SImelt_mean_ref + gl_melt_mean_ref
    total_ablation_std_ref = snowmeltstd_mean_ref + SImeltstd_mean_ref + gl_meltstd_mean_ref


    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')   
    time = np.arange(0,len(dates))

    # Plot the average
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1) + '\nREF_MODEL - ' + alt_model_name,fontsize=14)
    ax.set_ylabel('Difference (m w.e. d$^{-1}$)',fontsize=14)
    ax.plot(time,total_accumulation_ref-total_accumulation_alt,c='mediumblue',label='Accumulation',zorder=10,linewidth=3)    
    ax.plot(time,(total_ablation_ref)-(total_ablation_alt),c='red',label='Ablation',zorder=10,linewidth=3)    
    #plt.fill_between(time,(total_accumulation_alt-total_accumulation_std_alt),(total_accumulation_alt+total_accumulation_std_alt),color='mediumblue',alpha=0.35)    
    #plt.fill_between(time,(-total_ablation_alt-total_ablation_std_alt),(-total_ablation_alt+total_ablation_std_alt),color='red',alpha=0.35)
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax.grid()
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='grey')
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    ax.plot(time,np.linspace(-100,-200,len(time)),c='k',linewidth=3,label='Cumulative difference')
    
    #ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    #ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=14)
    #ax.plot(time,total_accumulation_ref,c='mediumblue',label='Accumulation')    
    #ax.plot(time,-total_ablation_ref,c='red',label='Ablation')    
    #plt.fill_between(time,(total_accumulation_ref-total_accumulation_std_ref),(total_accumulation_ref+total_accumulation_std_ref),color='mediumblue',alpha=0.35)    
    #plt.fill_between(time,(-total_ablation_ref-total_ablation_std_ref),(-total_ablation_ref+total_ablation_std_ref),color='red',alpha=0.35)

    
    massbal_alt = total_accumulation_alt - total_ablation_alt
    max_mb_alt = (total_accumulation_alt + total_accumulation_std_alt) - (total_ablation_alt - total_ablation_std_alt)
    min_mb_alt = (total_accumulation_alt - total_accumulation_std_alt) - (total_ablation_alt + total_ablation_std_alt)

    massbal_ref = total_accumulation_ref - total_ablation_ref

    #transition_date = dates[50:][np.where(np.cumsum(massbal_alt)[50:] <=0)[0][0]]
    ax0 = ax.twinx()
    ax0.plot(time,np.cumsum(massbal_ref)-np.cumsum(massbal_alt),c='k',linewidth=3,label='Cumulative difference')
    #ax0.plot(time,np.cumsum(massbal_ref),c='k',linewidth=3)
    #plt.fill_between(time,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='std. dev.')

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Difference (m w.e.)',fontsize=14)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    #ax0.legend(fontsize=14,loc='lower left')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='best', ncol=1, borderaxespad=0.19)
    fig.tight_layout()

def percent_difference_with_direction(A, B):
    abs_diff = abs(A - B)
    avg_magnitude = 0.5 * (abs(A) + abs(B))
    percent_diff = float('{:.2g}'.format((abs_diff / avg_magnitude) * 100))

    if A > B:
        direction = "Decrease"
    elif A < B:
        direction = "Increase"
    else:
        direction = "No Change"

    return percent_diff if direction == "Increase" else -percent_diff

def compare_hydrographs_differences(units,title,avg_years,all_years,lowerlim,upperlim,gl_icemelt_ref, snowmelt_runoff_ref, rain_runoff_ref, superimp_icemelt_ref,gl_icemelt_std_ref, snowmelt_runoff_std_ref, rain_runoff_std_ref, superimp_icemelt_std_ref, \
                                        gl_icemelt_alt, snowmelt_runoff_alt, rain_runoff_alt, superimp_icemelt_alt,gl_icemelt_std_alt, snowmelt_runoff_std_alt, rain_runoff_std_alt, superimp_icemelt_std_alt,area_map):
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt_ref[i][:365])
        snow_sum += np.array(snowmelt_runoff_ref[i][:365])
        rain_sum += np.array(rain_runoff_ref[i][:365])
        SI_sum += np.array(superimp_icemelt_ref[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std_ref[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std_ref[i][:365])
        rainstd_sum += np.array(rain_runoff_std_ref[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std_ref[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean_ref = np.array(ice_sum/len(avg_years))
    snow_mean_ref = np.array(snow_sum/len(avg_years))
    rain_mean_ref = np.array(rain_sum/len(avg_years))
    SI_mean_ref = np.array(SI_sum/len(avg_years))
    total_runoff_ref = ice_mean_ref + snow_mean_ref + rain_mean_ref + SI_mean_ref
    
    #icestd_mean_ref = np.array(icestd_sum/len(avg_years))
    #snowstd_mean_ref = np.array(snowstd_sum/len(avg_years))
    #rainstd_mean_ref = np.array(rainstd_sum/len(avg_years))
    #SIstd_mean_ref = np.array(SIstd_sum/len(avg_years))
    #total_std_ref = icestd_mean_ref + snowstd_mean_ref + rainstd_mean_ref + SIstd_mean_ref
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt_alt[i][:365])
        snow_sum += np.array(snowmelt_runoff_alt[i][:365])
        rain_sum += np.array(rain_runoff_alt[i][:365])
        SI_sum += np.array(superimp_icemelt_alt[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std_alt[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std_alt[i][:365])
        rainstd_sum += np.array(rain_runoff_std_alt[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std_alt[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean_alt = np.array(ice_sum/len(avg_years))
    snow_mean_alt = np.array(snow_sum/len(avg_years))
    rain_mean_alt = np.array(rain_sum/len(avg_years))
    SI_mean_alt = np.array(SI_sum/len(avg_years))
    total_runoff_alt = ice_mean_alt + snow_mean_alt + rain_mean_alt + SI_mean_alt
    
    #icestd_mean_alt = np.array(icestd_sum/len(avg_years))
    #snowstd_mean_alt = np.array(snowstd_sum/len(avg_years))
    #rainstd_mean_alt = np.array(rainstd_sum/len(avg_years))
    #SIstd_mean_alt = np.array(SIstd_sum/len(avg_years))
    #total_std_alt = icestd_mean_alt + snowstd_mean_alt + rainstd_mean_alt + SIstd_mean_alt
    
    time = np.arange(0,len(total_runoff_ref))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2
   
    icediff = percent_difference_with_direction(np.sum(ice_mean_ref),np.sum(ice_mean_alt))
    snowdiff = percent_difference_with_direction(np.sum(snow_mean_ref),np.sum(snow_mean_alt))
    raindiff = percent_difference_with_direction(np.sum(rain_mean_ref),np.sum(rain_mean_alt))
    SIdiff = percent_difference_with_direction(np.sum(SI_mean_ref),np.sum(SI_mean_alt))
    totaldiff = percent_difference_with_direction(np.sum(total_runoff_ref),np.sum(total_runoff_alt))
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1) + '\nREF_MODEL - ROUNCE_DEBRIS',fontsize=14)
    ax.plot(time,(ice_mean_ref-ice_mean_alt)*dc,c='turquoise',label='Glacier ice melt (' + str(icediff) + ' %)',linewidth=3)    
    ax.plot(time,(snow_mean_ref-snow_mean_alt)*dc,c='royalblue',label='Snow melt ('+ str(snowdiff) + ' %)',linewidth=3)   
    ax.plot(time,(rain_mean_ref-rain_mean_alt)*dc,c='deeppink',label='Rain (' + str(raindiff) + ' %)',linewidth=3)   
    ax.plot(time,(SI_mean_ref-SI_mean_alt)*dc,c='darkorange',label='Superimposed ice melt (' + str(SIdiff) + ' %)',linewidth=3)    
    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335,365,395])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',''],rotation=45,fontsize=14)
    ax.grid()
    #plt.xlim(0,395)
    plt.xlim(0,365)
    ax.set_ylim(lowerlim,upperlim)

    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)

    #ax0 = ax.twinx()
    #ax0.set_ylim(-cumu_runoff_upperlim,cumu_runoff_upperlim)
    #ax0.tick_params(axis='y',labelsize=14)
    #ax0.margins(x=0)
    
    if units=='mwe':
        ax.set_ylabel('Difference (m w.e. day$^{-1}$)',fontsize=14)
    #    ax0.set_ylabel('Cumulative difference (m w.e. a$^{-1}$)',fontsize=14)
    elif units=='m3':
        ax.set_ylabel('Difference (m$^3$ s$^{-1}$)',fontsize=14)
    #    ax0.set_ylabel('Cumulative difference (km$^3$ a$^{-1}$)',fontsize=14)
    
    #width=10
# ================================================================================================================
     #plt.vlines(395, ymin=(np.cumsum(total_runoff_ref-total_std_ref)*yc)[-1], ymax=(np.cumsum(total_runoff_ref+total_std_ref)*yc)[-1], linewidth=width, color='k',alpha=0.3)
    #plt.hlines(y=((np.cumsum(total_runoff_ref)-np.cumsum(total_runoff_alt))*yc)[-1], xmin=380 - width / 2, xmax=380 + width / 2,colors='k',linewidth=4,label='Cumulative total runoff')
     
     #plt.vlines(395, ymin=(np.cumsum(ice_mean_ref-icestd_mean_ref)*yc)[-1], ymax=(np.cumsum(ice_mean_ref+icestd_mean_ref)*yc)[-1], linewidth=width, color='turquoise',alpha=0.3)
    #plt.hlines(y=((np.cumsum(ice_mean_ref)-np.cumsum(ice_mean_alt))*yc)[-1], xmin=380 - width / 2, xmax=380 + width / 2,colors='turquoise',linewidth=4)
     
     #plt.vlines(395, ymin=(np.cumsum(snow_mean_ref-snowstd_mean_ref)*yc)[-1], ymax=(np.cumsum(snow_mean_ref+snowstd_mean_ref)*yc)[-1], linewidth=width, color='royalblue',alpha=0.3)
    #plt.hlines(y=((np.cumsum(snow_mean_ref)-np.cumsum(snow_mean_alt))*yc)[-1], xmin=380 - width / 2, xmax=380 + width / 2,colors='royalblue',linewidth=4)
     
     #plt.vlines(395, ymin=(np.cumsum(rain_mean_ref-rainstd_mean_ref)*yc)[-1], ymax=(np.cumsum(rain_mean_ref+rainstd_mean_ref)*yc)[-1], linewidth=width, color='deeppink',alpha=0.3)
    #plt.hlines(y=((np.cumsum(rain_mean_ref)-np.cumsum(rain_mean_alt))*yc)[-1], xmin=380 - width / 2, xmax=380 + width / 2,colors='deeppink',linewidth=4)
     
     #plt.vlines(395, ymin=(np.cumsum(SI_mean_ref-SIstd_mean_ref)*yc)[-1], ymax=(np.cumsum(SI_mean_ref+SIstd_mean_ref)*yc)[-1], linewidth=width, color='darkorange',alpha=0.3)
    #plt.hlines(y=((np.cumsum(SI_mean_ref)-np.cumsum(SI_mean_alt))*yc)[-1], xmin=380 - width / 2, xmax=380 + width / 2,colors='darkorange',linewidth=4)
 
    ax.hlines(y=-1000, xmin=375 - 10 / 2, xmax=375 + 10 / 2,colors='k',linewidth=4,label='Total runoff (' + str(totaldiff) + ' %)')    
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),loc='upper left',fontsize=14, ncol=1, borderaxespad=0.19)
    

def surface_elevation_profiles():
    # profiles of surface elevation across the medial moraines
    
    # Load initial surface elev profile (1979 DEM)
    z_initial = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Zgrids/DEM_KRH_1979.txt')
    debris_thickness_map = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/DebrisThickness/KRH_debrismap.txt')
    
    # select transect of interest
    plt.figure(figsize=(8,5))
    #plt.contourf(Xgrid,Ygrid,debris_map)
    plt.contourf(Xgrid[:,250:350],Ygrid[:,250:350],debris_map[:,250:350],levels=np.linspace(0,1,21),cmap = 'PuOr')
    legend = plt.colorbar()
    plt.plot(Xgrid[:,295][np.where(KRH_tributaries[:,295] ==5 )],Ygrid[:,295][np.where(KRH_tributaries[:,295] ==5 )],linewidth=3,color='gold')
    plt.plot(Xgrid[:,315][np.where(KRH_tributaries[:,315] ==5 )],Ygrid[:,315][np.where(KRH_tributaries[:,315] ==5 )],linewidth=3,color='cyan')
    plt.plot(Xgrid[26,:][np.where(KRH_tributaries[26,:] ==5 )],Ygrid[26,:][np.where(KRH_tributaries[26,:] ==5 )],linewidth=3,color='deeppink')
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.xlim(610000,635000)
    plt.ylim(6730000,6750000)
    legend.ax.set_ylabel('Debris thickness (m)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    legend.ax.tick_params(labelsize=14)
    #plt.axis('equal')
    
    # Load melt scaling map for both models
    debris_m_ref = debris('Sub-debris melt scaling','F:/Mass Balance Model/Kaskawulsh-Mass-Balance/DebrisThickness/KRH_debrismap.txt',Sfc,2.0277,2.1717,0.006,0.019,11.0349260206858,1.98717418666925)
    debris_m_rounce = debris('Sub-debris melt scaling','F:/Mass Balance Model/Kaskawulsh-Mass-Balance/DebrisThickness/KRH_debrismap.txt',Sfc,2.9200183038347,6.62306028549649,0.02,0.1268126812681268,11.0349260206858,1.98717418666925)
    
    
    plt.figure(figsize=(8,8))
    plt.subplot(3,2,1)
    plt.title('Transect 1',fontsize=14,color='gold')
    debris_thickness = debris_thickness_map[:,295][np.where(KRH_tributaries[:,295] ==5 )]
    plt.plot(np.arange(0,len(debris_thickness)*0.2,0.2),debris_thickness,label='1979 surface',color='k',linewidth=3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Debris thickness (m)',fontsize=14)
    plt.grid()
    plt.ylim(0,0.6)
    plt.subplot(3,2,3)
    plt.title('Transect 2',fontsize=14,color='cyan')
    debris_thickness = debris_thickness_map[:,315][np.where(KRH_tributaries[:,315] ==5 )]
    plt.plot(np.arange(0,len(debris_thickness)*0.2,0.2),debris_thickness,label='1979 surface',color='k',linewidth=3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Debris thickness (m)',fontsize=14)
    plt.grid()
    plt.ylim(0,0.6)
    plt.subplot(3,2,5)
    plt.title('Transect 3',fontsize=14,color='deeppink')
    debris_thickness = debris_thickness_map[26,:][np.where(KRH_tributaries[26,:] ==5 )]
    plt.plot(np.arange(0,len(debris_thickness)*0.2,0.2),debris_thickness,label='1979 surface',color='k',linewidth=3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Debris thickness (m)',fontsize=14)
    plt.grid()
    plt.ylim(0,0.6)
    
    plt.subplot(3,2,2)
    plt.title('Transect 1',fontsize=14,color='gold')
    debris = debris_m_ref[:,295][np.where(KRH_tributaries[:,295] ==5 )]
    debris2 = debris_m_rounce[:,295][np.where(KRH_tributaries[:,295] ==5 )]
    plt.plot(np.arange(0,len(debris)*0.2,0.2),debris,label='REF_MODEL',color='cornflowerblue',linewidth=3)
    plt.plot(np.arange(0,len(debris2)*0.2,0.2),debris2,label='ROUNCE_DEBRIS',color='firebrick',linewidth=3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Melt scaling',fontsize=14)
    plt.grid()
    plt.ylim(0,2.5)
    plt.legend(fontsize=12)
    plt.subplot(3,2,4)
    plt.title('Transect 2',fontsize=14,color='cyan')
    debris = debris_m_ref[:,315][np.where(KRH_tributaries[:,315] ==5 )]
    debris2 = debris_m_rounce[:,315][np.where(KRH_tributaries[:,315] ==5 )]
    plt.plot(np.arange(0,len(debris)*0.2,0.2),debris,label='REF_MODEL',color='cornflowerblue',linewidth=3)
    plt.plot(np.arange(0,len(debris2)*0.2,0.2),debris2,label='ROUNCE_DEBRIS',color='firebrick',linewidth=3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Melt scaling',fontsize=14)
    plt.grid()
    plt.ylim(0,2.5)
    plt.subplot(3,2,6)
    plt.title('Transect 3',fontsize=14,color='deeppink')
    debris = debris_m_ref[26,:][np.where(KRH_tributaries[26,:] ==5 )]
    debris2 = debris_m_rounce[26,:][np.where(KRH_tributaries[26,:] ==5 )]
    plt.plot(np.arange(0,len(debris)*0.2,0.2),debris,label='REF_MODEL',color='cornflowerblue',linewidth=3)
    plt.plot(np.arange(0,len(debris2)*0.2,0.2),debris2,label='ROUNCE_DEBRIS',color='firebrick',linewidth=3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Melt scaling',fontsize=14)
    plt.grid()
    plt.ylim(0,2.5)
    
    plt.tight_layout()
    
    plt.figure(figsize=(9,10))
    plt.subplot(3,1,1)
    plt.title('Transect 1',fontsize=14,color='gold')
    initial_surface = z_initial[:,295][np.where(KRH_tributaries[:,295] ==5 )]
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),initial_surface,label='1979 surface',color='k',linewidth=3)
    mb = np.zeros((z_initial[:,295][np.where(KRH_tributaries[:,295] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_ref[year-years[0]][:,295][np.where(KRH_tributaries[:,295] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z,label='2022 surface (REF_MODEL)',color='cornflowerblue',linewidth=3)
    mb = np.zeros((z_initial[:,295][np.where(KRH_tributaries[:,295] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_rounce[year-years[0]][:,295][np.where(KRH_tributaries[:,295] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z,label='2022 surface (ROUNCE_DEBRIS)',color='firebrick',linewidth=3)
    plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Surface Elevation (m a.s.l.)',fontsize=14)
    plt.grid()
    
    plt.subplot(3,1,2)
    plt.title('Transect 2',fontsize=14,color='cyan')
    initial_surface = z_initial[:,315][np.where(KRH_tributaries[:,315] ==5 )]
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),initial_surface,label='1979 surface',color='k',linewidth=3)
    mb = np.zeros((z_initial[:,315][np.where(KRH_tributaries[:,315] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_ref[year-years[0]][:,315][np.where(KRH_tributaries[:,315] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z,label='2022 surface (REF_MODEL)',color='cornflowerblue',linewidth=3)
    mb = np.zeros((z_initial[:,315][np.where(KRH_tributaries[:,315] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_rounce[year-years[0]][:,315][np.where(KRH_tributaries[:,315] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z,label='2022 surface (ROUNCE_DEBRIS)',color='firebrick',linewidth=3)
    #plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Surface Elevation (m a.s.l.)',fontsize=14)
    plt.grid()
    
    plt.subplot(3,1,3)
    plt.title('Transect 3',fontsize=14,color='deeppink')
    initial_surface = z_initial[26,:][np.where(KRH_tributaries[26,:] ==5 )]
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),initial_surface,label='1979 surface',color='k',linewidth=3)
    mb = np.zeros((z_initial[26,:][np.where(KRH_tributaries[26,:] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_ref[year-years[0]][26,:][np.where(KRH_tributaries[26,:] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z,label='2022 surface (REF_MODEL)',color='cornflowerblue',linewidth=3)
    mb = np.zeros((z_initial[26,:][np.where(KRH_tributaries[26,:] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_rounce[year-years[0]][26,:][np.where(KRH_tributaries[26,:] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z,label='2022 surface (ROUNCE_DEBRIS)',color='firebrick',linewidth=3)
    #plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Surface Elevation (m a.s.l.)',fontsize=14)
    plt.grid()
    
    plt.tight_layout()
    
    plt.figure(figsize=(9,10))
    plt.subplot(3,1,1)
    plt.title('Transect 1',fontsize=14,color='gold')
    initial_surface = z_initial[:,295][np.where(KRH_tributaries[:,295] ==5 )]
    #plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),initial_surface,label='1979 surface',color='k',linewidth=3)
    mb = np.zeros((z_initial[:,295][np.where(KRH_tributaries[:,295] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_ref[year-years[0]][:,295][np.where(KRH_tributaries[:,295] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z-initial_surface,label='2022 surface (REF_MODEL)',color='cornflowerblue',linewidth=3)
    mb = np.zeros((z_initial[:,295][np.where(KRH_tributaries[:,295] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_rounce[year-years[0]][:,295][np.where(KRH_tributaries[:,295] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z-initial_surface,label='2022 surface (ROUNCE_DEBRIS)',color='firebrick',linewidth=3)
    plt.legend(fontsize=14,loc='upper left')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Elevation Change (m a.s.l.)',fontsize=14)
    plt.grid()
    plt.ylim(-16,0)
    
    plt.subplot(3,1,2)
    plt.title('Transect 2',fontsize=14,color='cyan')
    initial_surface = z_initial[:,315][np.where(KRH_tributaries[:,315] ==5 )]
    #plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),initial_surface,label='1979 surface',color='k',linewidth=3)
    mb = np.zeros((z_initial[:,315][np.where(KRH_tributaries[:,315] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_ref[year-years[0]][:,315][np.where(KRH_tributaries[:,315] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z-initial_surface,label='2022 surface (REF_MODEL)',color='cornflowerblue',linewidth=3)
    mb = np.zeros((z_initial[:,315][np.where(KRH_tributaries[:,315] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_rounce[year-years[0]][:,315][np.where(KRH_tributaries[:,315] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z-initial_surface,label='2022 surface (ROUNCE_DEBRIS)',color='firebrick',linewidth=3)
    #plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Elevation Change (m a.s.l.)',fontsize=14)
    plt.grid()
    plt.ylim(-16,0)
    
    plt.subplot(3,1,3)
    plt.title('Transect 3',fontsize=14,color='deeppink')
    initial_surface = z_initial[26,:][np.where(KRH_tributaries[26,:] ==5 )]
    #plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),initial_surface,label='1979 surface',color='k',linewidth=3)
    mb = np.zeros((z_initial[26,:][np.where(KRH_tributaries[26,:] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_ref[year-years[0]][26,:][np.where(KRH_tributaries[26,:] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z-initial_surface,label='2022 surface (REF_MODEL)',color='cornflowerblue',linewidth=3)
    mb = np.zeros((z_initial[26,:][np.where(KRH_tributaries[26,:] ==5 )].shape))
    for year in years:
        mb =+ massbal_dist_rounce[year-years[0]][26,:][np.where(KRH_tributaries[26,:] ==5 )]
        new_z = initial_surface + mb
    plt.plot(np.arange(0,len(initial_surface)*0.2,0.2),new_z-initial_surface,label='2022 surface (ROUNCE_DEBRIS)',color='firebrick',linewidth=3)
    #plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Distance (km)',fontsize=14)
    plt.ylabel('Elevation Change (m a.s.l.)',fontsize=14)
    plt.grid()
    plt.ylim(-16,0)
    
    plt.tight_layout()


def Bmod_vs_Bcal(avg_years,KRH_tributaries,KW_fluxgates,dist_massbal):
    all_years = np.arange(1979,2022)
    MB_ARRAY = np.empty((len(avg_years),230,329))
    for year in avg_years:
        print(year)
        MB_ARRAY[year-avg_years[0],:,:] = dist_massbal[year-all_years[0]]
        
    bmod = np.nanmean(MB_ARRAY,axis=0)
    
    
    print('Trunk',np.round(np.mean(bmod[np.where(KRH_tributaries==5)]),2))
    print('NA',np.round(np.mean(bmod[np.where(KRH_tributaries==4)]),2))
    print('CA',np.round(np.mean(bmod[np.where(KRH_tributaries==3)]),2))
    print('SW',np.round(np.mean(bmod[np.where(KRH_tributaries==2)]),2))
    print('SA',np.round(np.mean(bmod[np.where(KRH_tributaries==1)]),2))
    print('All tribs',np.round(np.mean(bmod[np.where(KRH_tributaries<=4)]),2))
    print('All KW', np.round(np.mean(bmod[np.isfinite(KRH_tributaries)]),2))
    #    0 = KW0
    #    1 = NA
    #    2 = CA
    #    3 = SW
    #    4 = SA
    #    5 = KW5
    #    6 = KW4
    #    7 = KW3
    #    8 = KW2
    #    9 = KW1
    print('KW0',np.round(np.mean(bmod[np.where(KW_fluxgates==0)]),2))
    print('KW1',np.round(np.mean(bmod[np.where(KW_fluxgates==9)]),2))
    print('KW2',np.round(np.mean(bmod[np.where(KW_fluxgates==8)]),2))
    print('KW3',np.round(np.mean(bmod[np.where(KW_fluxgates==7)]),2))
    print('KW4',np.round(np.mean(bmod[np.where(KW_fluxgates==6)]),2))
    print('KW5',np.round(np.mean(bmod[np.where(KW_fluxgates==5)]),2))
    print('KW1 + KW2',np.round(np.mean(bmod[np.where(KW_fluxgates>=8)]),2))
    print('KW4 + KW5',np.round(np.mean(bmod[np.where((KW_fluxgates==6) | (KW_fluxgates==5))]),2))
    
    Zones = ['NA','CA','SW','SA','Main\nTrunk','KW0','KW1','KW2','KW3','KW4','KW5','KW\n1-2','KW\n4-5','All\nTs','Glacier-\nwide']
    bmod_KR = [np.round(np.mean(bmod[np.where(KRH_tributaries==4)]),2),\
        np.round(np.mean(bmod[np.where(KRH_tributaries==3)]),2),\
        np.round(np.mean(bmod[np.where(KRH_tributaries==2)]),2),\
        np.round(np.mean(bmod[np.where(KRH_tributaries==1)]),2),\
        np.round(np.mean(bmod[np.where(KRH_tributaries==5)]),2),\
        np.round(np.mean(bmod[np.where(KW_fluxgates==0)]),2),\
        np.round(np.mean(bmod[np.where(KW_fluxgates==9)]),2),\
        np.round(np.mean(bmod[np.where(KW_fluxgates==8)]),2),\
        np.round(np.mean(bmod[np.where(KW_fluxgates==7)]),2),\
        np.round(np.mean(bmod[np.where(KW_fluxgates==6)]),2),\
        np.round(np.mean(bmod[np.where(KW_fluxgates==5)]),2),\
        np.round(np.mean(bmod[np.where(KW_fluxgates>=8)]),2),\
        np.round(np.mean(bmod[np.where((KW_fluxgates==6) | (KW_fluxgates==5))]),2),\
        np.round(np.mean(bmod[np.where(KRH_tributaries<=4)]),2),\
        np.round(np.mean(bmod[np.isfinite(KRH_tributaries)]),2)]
        
    bmod_EMY = [0.50,0.25,1.04,0.28,-4.01,-5.62,-5.23,-3.92,-4.26,-3.84,-2.70,-4.64,-3.16,0.41,-0.42]
    bcal = [0.74,0.31,0.03,-0.14,-3.41,-4.82,-2.89,-5.81,-2.93,-1.69,-3.28,-4.22,-2.64,0.25,-0.46]
        
    # Load text files with std devs for each zone:
    bmod_KR_stddev = np.loadtxt('D:/Model Runs/REF_MODEL/Processing/Bmod_Zones_REFMODEL_stddev.txt')
    
    # Bar chart
    bar_width = 0.25
    index = np.arange(len(Zones))
    
    fig, ax = plt.subplots(figsize=(12, 5.5))
    
    bar1 = ax.bar(index - bar_width, bmod_EMY, bar_width, label='$\dot{B}_{mod}$\n(Young et al. 2021)',color='lightsteelblue',zorder=10)
    bar2 = ax.bar(index, bcal, bar_width, label='$\dot{B}_{cal}$\n(Young et al. 2021)',color='mediumpurple',zorder=10)
    bar3 = ax.bar(index + bar_width, bmod_KR, bar_width, label='$\dot{B}_{mod}$\n(This study)',color='navy',zorder=10)
    ax.errorbar(index + bar_width, bmod_KR,bmod_KR_stddev,fmt='none',zorder=20,capsize=5,color='dimgrey')
    
    # Customize the plot
    #ax.set_xlabel('Categories')
    ax.set_ylabel('$\dot{B}_{sfc}$ (m w.e. a$^{-1}$)',fontsize=14)
    ax.tick_params(axis='both',labelsize=14)
    #ax.set_title('Bar Chart with 14 Categories and 3 Values Each')
    ax.set_xticks(index)
    ax.set_xticklabels(Zones, rotation=0)
    ax.legend(fontsize=14,loc='lower right',ncol=3,bbox_to_anchor=(1,1))
    ax.grid()
    ax.set_xlim(-2.5,15)
    #plt.text(-2.2,1.5,'a)',fontsize=18,weight='bold')
    #ax.set_ylim(-7,1.25)
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(9,5))
    plt.axis('equal')
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.7,alpha=0.8)
    plt.contour(Xgrid,Ygrid,KW_fluxgates,levels=(0,4,5,6,7,8,9,10),colors='navy',linewidths=3,alpha=1)
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.axis('off')
    #plt.text(5.8e5,6.71e6,'b)',fontsize=18,weight='bold')
    plt.tight_layout()
    
    # Plot the MAE between Bmod(this study) and Bcal, and Bmod(young et al) and Bcal
    Bmod_Young_MAE = np.abs(np.array(bcal)- np.array(bmod_EMY))
    Bmod_Robinson_MAE = np.abs(np.array(bcal)- np.array(bmod_KR))
    
    bar_width = 0.3
    fig, ax = plt.subplots(figsize=(12, 2.6))
    bar1 = ax.bar(index - 0.5*bar_width, Bmod_Young_MAE, bar_width, label='$\dot{B}_{mod}$\n(Young et al. 2021)',color='lightsteelblue',zorder=10)
    bar2 = ax.bar(index + 0.5*bar_width, Bmod_Robinson_MAE, bar_width, label='$\dot{B}_{mod}$\n(This study)',color='navy',zorder=10)
    
    ax.set_ylabel('Abs. diff. (m w.e. a$^{-1}$)',fontsize=14)
    ax.tick_params(axis='both',labelsize=14)
    #ax.set_title('Bar Chart with 14 Categories and 3 Values Each')
    ax.set_xticks(index)
    ax.set_xticklabels([' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '], rotation=0)
    ax.xaxis.tick_top()
    #ax.legend(fontsize=14,loc='lower right',ncol=3,bbox_to_anchor=(1,1))
    ax.grid()
    ax.set_xlim(-2.5,15)
    
    ax.set_yticks(np.arange(0,2.6,0.5))
    
    #plt.text(-2.2,1.5,'a)',fontsize=18,weight='bold')
    #ax.set_ylim(-7,1.25)
    plt.tight_layout()
    plt.show()    
    
    return bmod_EMY, bmod_KR
    
    #print(Bmod_Young_MAE)
    #print(Bmod_Robinson_MAE)
def Calculate_mass_loss_rates():  
# Kaskawulsh mass loss rates
    eighties = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,1990),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref,KRH_tributaries)/10
    nineties = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1990,2000),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref,KRH_tributaries)/10
    twothousands = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2000,2010),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref,KRH_tributaries)/10
    tens = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2010,2020),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref,KRH_tributaries)/10
    twenties = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2020,2022),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref,KRH_tributaries)/2
    allyears = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref,KRH_tributaries)/42
    
    Otherice = np.array((Sfc))
    Otherice[np.isfinite(KRH_tributaries)] = np.nan
    Otherice[Otherice==1] = np.nan
    eighties = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,1990),0.07,17,accumulation_nonKWice_ref, refrozen_rain_nonKWice_ref, gl_icemelt_nonKWice_ref, netsnowmelt_nonKWice_ref, superimp_icemelt_nonKWice_ref,Otherice)/10
    nineties = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1990,2000),0.07,17,accumulation_nonKWice_ref, refrozen_rain_nonKWice_ref, gl_icemelt_nonKWice_ref, netsnowmelt_nonKWice_ref, superimp_icemelt_nonKWice_ref,Otherice)/10
    twothousands = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2000,2010),0.07,17,accumulation_nonKWice_ref, refrozen_rain_nonKWice_ref, gl_icemelt_nonKWice_ref, netsnowmelt_nonKWice_ref, superimp_icemelt_nonKWice_ref,Otherice)/10
    tens = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2010,2020),0.07,17,accumulation_nonKWice_ref, refrozen_rain_nonKWice_ref, gl_icemelt_nonKWice_ref, netsnowmelt_nonKWice_ref, superimp_icemelt_nonKWice_ref,Otherice)/10
    twenties = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2020,2022),0.07,17,accumulation_nonKWice_ref, refrozen_rain_nonKWice_ref, gl_icemelt_nonKWice_ref, netsnowmelt_nonKWice_ref, superimp_icemelt_nonKWice_ref,Otherice)/2
    allyears = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_nonKWice_ref, refrozen_rain_nonKWice_ref, gl_icemelt_nonKWice_ref, netsnowmelt_nonKWice_ref, superimp_icemelt_nonKWice_ref,Otherice)/42
    
    # All glacierized area mass loss rates
    eighties2 = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,1990),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area)
    nineties2 = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1990,2000),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area)
    twothousands2 = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2000,2010),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area)
    tens2 = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2010,2020),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area)
    twenties2 = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(2020,2022),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area)
    allyears2 = cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area)
      
    print('total mass loss from 1980s, 1990s, 2000s, 2010s, 2020-2022, 1980-2022:',\
          eighties2,nineties2,twothousands2,tens2,twenties2,allyears2)
    
def Calculate_KRH_contribution_to_Alaska_massloss():
    KRH_ice_area = np.where(Sfc==0)[0].shape[0]*(0.2*0.2) # km2
    Alaska_ice_area = 86734  #km2, from https://www.glims.org/RGI/00_rgi60_TechnicalNote.pdf
    frac_KRH_in_Alaska = KRH_ice_area/Alaska_ice_area*100
    
    # what is the modelled mass loss in KRH for 2000--2019:
    KRH_2000_19_massloss = np.array(cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1999,2019),0.07,17,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area))[0]/20 #Gt
    frac_KRH_massloss = KRH_2000_19_massloss/-66.7*100 #Gt per yr
    print('KRH contributed',np.round(frac_KRH_massloss,2),'% to Alaska mass loss 2000-2019')
    
    KW_2000_19_massloss = np.array(cumulative_massbalance('Glacierized Area Mass Balance\n',years,np.arange(1999,2019),0.07,17,accumulation_kw_ref, refrozen_rain_kw_ref, gl_icemelt_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref,KRH_tributaries))[1]/20 # m w.e.
    print('Modelled mass loss from Kaskawulsh Glacier (2000-2019) = ',np.round(KW_2000_19_massloss,2),' m w.e. a$^{-1}$')

    
    alaska_mb = pd.read_csv('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/hugonnet_mb_Alaska.csv')
    AK_gl_subregions = np.array(alaska_mb['O2Region'])
    
    AK_gl_areas = np.array(alaska_mb['area'])
    stelias_gl_area = np.sum(AK_gl_areas[np.where(AK_gl_subregions==5)])
    
    AK_gl_mb_mwe = np.array(alaska_mb['mb_mwea'])
    AK_gl_mb_Gta = ((AK_gl_areas*1e6)*AK_gl_mb_mwe)/1e9 # area in m2 x mb in m w e = volume in m3 
    
    stelias_gl_mb_Gta = np.sum(AK_gl_mb_Gta[np.where(AK_gl_subregions==5)])
    frac_KRH_massloss_StElias = KRH_2000_19_massloss/stelias_gl_mb_Gta*100 #Gt per yr
    print('KRH contributed',np.round(frac_KRH_massloss_StElias,2),'% to St. Elias mass loss 2000-2019')
        
    frac_KRH_in_StElias = KRH_ice_area/stelias_gl_area*100
    print('KRH represents ',np.round(frac_KRH_in_StElias,2),' % of ice in the St. Elias mountains')

        
    
    