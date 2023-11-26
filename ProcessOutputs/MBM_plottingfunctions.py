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


def load_hydrologic_year(sim,year,varname,NARR_INPUTS,MODEL_OUTPUTS,Glacier_ID,stddev=False,NARRvar=False):
    '''
    Returns a given variable from Oct 1 -- Sept 30 (hydrological year)
    '''
    
    dates_yr1 = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')
    dates_yr2 = pd.date_range(start= str(year+1) + '-01-01 00:00:00',end= str(year+1) + '-12-31 21:00:00',freq='D')
        
    # Concatenate variable from Oct 1 of current year to Sept 30 of following year    
    if NARRvar == True:
        inMB1 = Dataset(os.path.join(NARR_INPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
        inMB2 = Dataset(os.path.join(NARR_INPUTS,varname + '_' + str(Glacier_ID) + '_' + str(year+1) + '.nc'),'r')
    else:
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
    domain (str): what area these outputs cover (e.g. krh, allgl, kw)
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
        
        ax.set_title(str(year)+'-'+str(year+1))
    
        ax.plot(np.arange(0,len(dates)),gl_icemelt[year-years[0]][:365]*dc,c='turquoise',label='Glacier ice melt')    
        ax.plot(np.arange(0,len(dates)),snowmelt_runoff[year-years[0]][:365]*dc,c='royalblue',label='Snow melt')
        ax.plot(np.arange(0,len(dates)),rain_runoff[year-years[0]][:365]*dc,c='deeppink',label='Rain')
        ax.plot(np.arange(0,len(dates)),superimp_icemelt[year-years[0]][:365]*dc,c='darkorange',label='Superimposed ice melt')
        ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
        ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45)
        ax.set_ylim(0,daily_runoff_upperlim)
        ax.grid()
        ax.margins(x=0)
    
        total_runoff = gl_icemelt[year-years[0]][:365] + snowmelt_runoff[year-years[0]][:365] + rain_runoff[year-years[0]][:365] + superimp_icemelt[year-years[0]][:365]
    
        ax0 = ax.twinx()
        ax0.plot(np.arange(0,len(dates)),np.cumsum(total_runoff)*yc,c='k',label='cumulative runoff')
        ax0.set_ylim(0,cumu_runoff_upperlim)
        ax0.margins(x=0)
        
        if units=='mwe':
            ax.set_ylabel('Runoff (m w.e. day$^{-1}$)',fontsize=11)
            ax0.set_ylabel('Cumulative Runoff (m w.e. a$^{-1}$)',fontsize=11)
        elif units=='m3':
            ax.set_ylabel('Runoff (m$^3$ s$^{-1}$)',fontsize=11)
            ax0.set_ylabel('Cumulative Runoff (km$^3$ a$^{-1}$)',fontsize=11)        
        
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
    
    
def distributed_average_mass_balance(model_name,avg_years,all_years,dist_massbal,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        print(year)
        MB_ARRAY[year-avg_years[0],:,:] = dist_massbal[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.mean((np.nanmean(MB_ARRAY,axis=0))[np.isfinite(KRH_tributaries)]),2)
    print(np.nanmin(np.nanmean(MB_ARRAY,axis=0)),np.nanmax(np.nanmean(MB_ARRAY,axis=0)))  
    
    plt.figure(figsize=(9,5))
    plt.title(model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),cmap='massbal',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,5,2))
    legend.ax.set_ylabel('Mass Balance (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=1)
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nmass balance\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
    plt.tight_layout()
    

def distributed_runoff(avg_years,all_years,netsnowmelt_dist,gl_icemelt_dist,superimp_icemelt_dist,rain_runoff_dist,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KRH_tributaries):

    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] = netsnowmelt_dist[year-all_years[0]] + gl_icemelt_dist[year-all_years[0]] + superimp_icemelt_dist[year-all_years[0]] + rain_runoff_dist[year-all_years[0]]
        
    kaskawulsh_mb = np.round(np.nanmean(np.nanmean(MB_ARRAY,axis=0)),2)
    
    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))
    
    plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,Runoff,cmap=cmocean.cm.ice_r,levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.linspace(0,20,11))
    legend.ax.set_ylabel('Total runoff (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=14)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=3,alpha=1)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=0.8,linestyles = 'dashed')
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nMean runoff\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)
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
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),fontsize=14,bbox_to_anchor=(0.98,1.04), ncol=2, borderaxespad=0.19)
    fig.suptitle(title,fontsize=14,y=1.01) 
    fig.tight_layout()
    
    
def massbalance_timeseries_average(title,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt, accumulation_std, refrozen_rain_std, netsnowmelt_std, superimp_icemelt_std, gl_icemelt_std):
    
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

    # Plot the average
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=14)
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
    ax0.plot(time,np.cumsum(massbal),c='k',label='Cumulative balance\n$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    #plt.fill_between(time,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='std. dev.')

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

def cumulative_massbalance(title,years,plotting_years,daily_mb_abs_lim,cumu_mb_abs_lim,accumulation,refrozen_rain,gl_icemelt,netsnowmelt,SImelt):
    
    year_i_index = plotting_years[0] - years[0]
    year_f_index = len(years) - (years[-1] - plotting_years[-1])
    #year_i_index = 0
    #year_f_index = 43
    
    all_accumulation = np.concatenate(accumulation[year_i_index:year_f_index], axis=0) + np.concatenate(refrozen_rain[year_i_index:year_f_index], axis=0)
    all_ablation = np.concatenate(gl_icemelt[year_i_index:year_f_index], axis=0) + np.concatenate(netsnowmelt[year_i_index:year_f_index], axis=0) + np.concatenate(SImelt[year_i_index:year_f_index], axis=0)
    
    acc_final = all_accumulation[~np.isnan(all_accumulation)]
    abl_final = all_ablation[~np.isnan(all_ablation)]
    
    massbal = acc_final - abl_final    
    
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
    ax0.plot(dates,np.cumsum(massbal),c='k',linewidth=3,label='Cumulative Mass Balance')
    #plt.fill_between(time,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='std. dev.')

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
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.tick_params(axis='x',labelsize=14)
    ax.margins(x=0.01)
    
    ax.bar(plotting_years,massbal,color='k',linewidth=3,label='Mass Balance',zorder=20)
        
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
    
def annualrunoff_stackedbar(title,years,netsnowmelt,gl_icemelt,superimp_icemelt,rain_runoff,area_map):
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
    m, b = np.polyfit(years,total_runoff,deg=1)
    
    width = 0.75       # the width of the bars: can also be len(x) sequence
    
    fig, ax = plt.subplots(figsize=(12,5))
    ax.bar(years, np.array(yrly_rain)*area, width, label='Rain',color='deeppink',zorder=5)
    ax.bar(years, np.array(yrly_SImelt)*area, width, bottom=np.array(yrly_rain)*area,label='Superimposed ice melt',color='darkorange',zorder=5)
    ax.bar(years, np.array(yrly_snowmelt)*area, width, bottom=(np.array(yrly_rain) + np.array(yrly_SImelt))*area,label='Snow melt',color='royalblue',zorder=5)
    ax.bar(years, np.array(yrly_icemelt)*area, width, bottom=(np.array(yrly_rain) + np.array(yrly_SImelt) + np.array(yrly_snowmelt))*area,label='Glacier ice melt',color='turquoise',zorder=5)
    plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2,zorder=20,label='Trend = +' + str(np.round(m,2))+ ' km$^{3}$ a$^{-1}$')
    
    ax.set_ylabel('Runoff (km$^{3}$ a$^{-1}$)',fontsize=14)
    plt.xticks(np.arange(years[0],years[-1]+1,5),fontsize=14,rotation=45)
    plt.yticks(fontsize=14)
    plt.ylim(0,2.75)
    #ax.set_title('Monthly Runoff (2007-2018) \n Debris Case',fontsize = 15)
    ax.legend(fontsize=14)
    plt.grid(zorder=-30)
    plt.title(title,fontsize=14)
    
    plt.tight_layout() 

    
def date_of_zero_balance(title,years,mb):
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
    for i in np.arange(int(np.nanmin(transition_DOYs)),int(np.nanmax(transition_DOYs))+1,10):
        monthday.append( str(dates[i])[5:10])
          
    m, b = np.polyfit(years[np.isfinite(transition_DOYs)],np.array(transition_DOYs)[np.isfinite(transition_DOYs)],deg=1)
    
    plt.figure(figsize=(7.2,3))
    plt.title(title + '\nDate of zero balance',fontsize=14)
    plt.scatter(years,transition_DOYs,c='darkblue',s=60)
    plt.yticks(np.arange(int(np.nanmin(transition_DOYs)),int(np.nanmax(transition_DOYs))+1,10),monthday,fontsize=14)
    plt.ylim(int(np.nanmin(transition_DOYs))-10,int(np.nanmax(transition_DOYs))+10)
    plt.xticks(fontsize=14)
    plt.ylabel('date of $\dot{B}$ = 0',fontsize=14)
    plt.grid()
    plt.plot(years,m*years + b,linestyle='--',c='k',linewidth=2,label='Trend = ' + str(np.round(m,2)) + ' days a$^{-1}$')
    plt.legend(fontsize=14,loc='upper right')
    
    width=8.5
    for year in years:
        #print(year,transition_DOYs[year-years[0]])
        if np.isnan(transition_DOYs[year-years[0]]):
            plt.vlines(year, ymin=0, ymax=365, linewidth=width, color='k',alpha=0.3)

    
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
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,5))
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    #ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    ax.plot(time,ice_mean*dc,c='turquoise',label='Glacier ice melt')    
    ax.plot(time,snow_mean*dc,c='royalblue',label='Snow melt')
    ax.plot(time,rain_mean*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean*dc,c='darkorange',label='Superimposed ice melt')
    
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
        ax0.set_ylabel('Cumulative Runoff (km$^3$ a$^{-1}$)',fontsize=14)
    
    width=10
    plt.vlines(375, ymin=(np.cumsum(total_runoff-total_std)*yc)[-1], ymax=(np.cumsum(total_runoff+total_std)*yc)[-1], linewidth=width, color='k',alpha=0.3)
    plt.hlines(y=(np.cumsum(total_runoff)*yc)[-1], xmin=375 - width / 2, xmax=375 + width / 2,colors='k',linewidth=4,label='Cumulative total runoff')
    
    plt.vlines(385, ymin=(np.cumsum(ice_mean-icestd_mean)*yc)[-1], ymax=(np.cumsum(ice_mean+icestd_mean)*yc)[-1], linewidth=width, color='turquoise',alpha=0.3)
    plt.hlines(y=(np.cumsum(ice_mean)*yc)[-1], xmin=385 - width / 2, xmax=385 + width / 2,colors='turquoise',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(snow_mean-snowstd_mean)*yc)[-1], ymax=(np.cumsum(snow_mean+snowstd_mean)*yc)[-1], linewidth=width, color='royalblue',alpha=0.3)
    plt.hlines(y=(np.cumsum(snow_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='royalblue',linewidth=4)
    
    plt.vlines(405, ymin=(np.cumsum(rain_mean-rainstd_mean)*yc)[-1], ymax=(np.cumsum(rain_mean+rainstd_mean)*yc)[-1], linewidth=width, color='deeppink',alpha=0.3)
    plt.hlines(y=(np.cumsum(rain_mean)*yc)[-1], xmin=405 - width / 2, xmax=405 + width / 2,colors='deeppink',linewidth=4)
    
    plt.vlines(415, ymin=(np.cumsum(SI_mean-SIstd_mean)*yc)[-1], ymax=(np.cumsum(SI_mean+SIstd_mean)*yc)[-1], linewidth=width, color='darkorange',alpha=0.3)
    plt.hlines(y=(np.cumsum(SI_mean)*yc)[-1], xmin=415 - width / 2, xmax=415 + width / 2,colors='darkorange',linewidth=4)
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
    ax1 = plt.axes([0.05, 0.1, 0.45, 0.45])
    #box_plot_data = [cum_ice_melt, cum_snow_melt, cum_superimposed_ice_melt, cum_rain]
    #box_colors = ['blue', 'orange', 'green', 'red']
    #box = ax2.boxplot(box_plot_data, positions=[1, 2, 3, 4], widths=0.15, patch_artist=True, showfliers=False)

    piechart_colours = ['turquoise','royalblue','deeppink','darkorange']
    fig = plt.figure(figsize=(5,5))
    
    total_runoff = np.sum(ice_mean + snow_mean + SI_mean + rain_mean)
    gl_ice_percent = round((np.sum(ice_mean)/total_runoff)*100,1)
    snow_percent = round((np.sum(snow_mean)/total_runoff)*100,1)
    SIice_percent = round((np.sum(SI_mean)/total_runoff)*100,1)
    rain_percent = round((np.sum(rain_mean)/total_runoff)*100,1)
    
    ax1.pie([gl_ice_percent,snow_percent,rain_percent,SIice_percent],autopct='%.1f%%',colors=piechart_colours, textprops={'fontsize': 14})
    #fig.suptitle(title,fontsize=14,y=1.01) 
    #plt.legend(['Glacier ice melt','Snow melt','Rain','Superimposed ice melt'],fontsize=14,ncol=1,bbox_to_anchor=(1.3,0.8))
    plt.tight_layout()
    
    fig.tight_layout()
    
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
    
    plt.figure(figsize=(9,5))
    plt.title('Total Runoff ' + str(avg_years[0]) + '--' + str(avg_years[-1]+1) +'\nREF_MODEL minus ' + model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,DIFF,cmap='runoff_diff',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,20))
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
    legend = plt.colorbar(ticks=np.arange(-20,20))
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
    legend = plt.colorbar(ticks=np.arange(-0.004,0.00001,0.001))
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




def compare_hydrographs(units,title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt_ref, snowmelt_runoff_ref, rain_runoff_ref, superimp_icemelt_ref,gl_icemelt_std_ref, snowmelt_runoff_std_ref, rain_runoff_std_ref, superimp_icemelt_std_ref, \
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
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','REF_MODEL','ROUNCE_DEBRIS'],rotation=45,fontsize=14)
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
    
    
def compare_date_of_zero_balance(years,mb_alt,mb_ref):
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
    
def massbalance_timeseries_difference(title,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim,  \
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
    ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1) + '\nREF_MODEL - ROUNCE_DEBRIS',fontsize=14)
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
    ax.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper left', ncol=1, borderaxespad=0.19)
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

def compare_hydrographs_differences(units,title,avg_years,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt_ref, snowmelt_runoff_ref, rain_runoff_ref, superimp_icemelt_ref,gl_icemelt_std_ref, snowmelt_runoff_std_ref, rain_runoff_std_ref, superimp_icemelt_std_ref, \
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
    ax.set_ylim(-daily_runoff_upperlim,daily_runoff_upperlim)

    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)

    #ax0 = ax.twinx()
    #ax0.set_ylim(-cumu_runoff_upperlim,cumu_runoff_upperlim)
    #ax0.tick_params(axis='y',labelsize=14)
    #ax0.margins(x=0)
    
    #if units=='mwe':
    #    ax.set_ylabel('Difference (m w.e. day$^{-1}$)',fontsize=14)
    #    ax0.set_ylabel('Cumulative difference (m w.e. a$^{-1}$)',fontsize=14)
    #elif units=='m3':
    #    ax.set_ylabel('Difference (m$^3$ s$^{-1}$)',fontsize=14)
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
    
