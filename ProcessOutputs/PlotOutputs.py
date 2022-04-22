# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 10:54:07 2021
PLOTTING SCRIPT
run this script after you have extracted all the MB, Melt, and Accumulation stuff you want using the "ExtractOutputs.py" script
@author: katierobinson
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\')
import Model_functions_ver4
from Model_functions_ver4 import regridXY_something
from tabulate import tabulate
import matplotlib.cm
from netCDF4 import Dataset
import matplotlib.ticker
# intialize the script#------------------------------------------------------

#where are the files you want to plot (directory not individual file)
Path2files = 'F:\Mass Balance Model\Catchment' #Baseline_NoDebris_Allparams
#Path2files = 'F:\\Mass Balance Model\\OUTPUTS\\Diagnostic\\Tshift=-1.7' #Baseline_NoDebris_Allparams
#Path2files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Kaskonly_R2S=1'
#Path2files = 'F:\\Mass Balance Model\\OUTPUTS\\NetZero_DebrisFree\\Tshift=-1.3'
Case = 'Tshift=-1.7'
debristitle = 'Debris Case'
r2stitle = 'R2S = 1.0$^\circ$C'
#File_glacier_in = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\kaskonly.txt'
File_glacier_in = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\kask_catchment.txt'
debris = False
kaskonly = False

# what do you want the script to actually plot
# MASS BALANCE
Plot_Mean_MB = False
Plot_MB_byYear = False
Plot_Acc_minus_Melt = False
Plot_FullTimeseries_MB = False
Mean_MB_Timeseries = False
Yearly_MB_HydrologicYear = False
Yearly_MB_SummerMelt = False

#ABLATION
Plot_NetMelt = False
Plot_NetMelt_byYear = False
Plot_FullTimeseries_Melt = False
Plot_NetAblation_IcevsSnow = False
Plot_TotalAblation = False
Plot_IceMelt = False
Plot_SnowMelt = False
Plot_IceMelt_byYear = False

#ACCUMULATION
Plot_Accumulation = False
Plot_Accumulation_byYear = False
Plot_Monthly_SnowvZ = False

#RAIN
Plot_Rain = False
Plot_Rain_byYear = False
Plot_FullTimeseries_Precip = False
Avg_Monthly_Rain = False
Plot_Monthly_Distributed_Rain = False
Plot_Monthly_RainvZ = False
Compare_BC_vs_NonBC_Rain = False
Master_BCvNoBC_Rain = False

#REFREEZING
Plot_Refreezing = False
Plot_Refreezing_byYear = False

#MISC.
Plots_vs_elevation = False
Calculate_Pmax = False

#RUNOFF
Plot_TotalRunoff_Yearly = False
Plot_FullTimeseries_Runoff = False
Mean_Runoff_PieChart = False
Yearly_Runoff_PieChart = False
Runoff_Table = False
Avg_Monthly_Runoff = False
Annual_Hydrograph = False
Mean_Runoff_Timeseries = False
Stacked_barchart_monthly_runoff = False

Mean_Runoff_Timeseries_withBCRain = False
Mean_Runoff_PieChart_withBCRain = False
Plot_TotalRunoff_Yearly_withBCRain = False
Yearly_Runoff_PieChart_withBCRain = False

#DOMAIN
Plot_Model_outline = False

#NARR
Plot_NARR_distributed_prcp = False
Plot_Monthly_NARRPrecip_vsZ = False
Plot_Monthly_Distributed_Temp = False
Plot_Monthly_TempvZ = False
Plot_NARR_distributed_acc_and_accVSz = True
Plot_mean_3hourly_temp = False


##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
#if debris is TRUE then we are working from the kaskonly_deb.txt file,
#if debris is FALSE then we are working from the kaskonly.txt file:
if kaskonly == True:
    if debris == True:
        print('kaskonly debris')
        Ix = glacier[:,3] 
        Iy = glacier[:,4] 
        Ih = glacier[:,2] 
        debris_array = glacier[:,6]
        #sfc_type = [:,]
    else:
        print('kaskonly no debris')
        Ix = glacier[:,2]
        Iy = glacier[:,3] 
        Ih = glacier[:,1]
else:
    print('full catchment')
    Ix = glacier[:,4]
    Iy = glacier[:,5] 
    Ih = glacier[:,6]
    sfc_type = glacier[:,8]
        
Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)

# ----------------------------------------------------------------------------
start_year = 2006
end_year = 2018
years = []
year = start_year
while year <= end_year:
  years.append(year)
  year = year + 1
#-----------------------------------------------------------------------------
File_temp_in = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\Tempkaskonly2007.nc'
inT = Dataset(File_temp_in, "r")
T_var = 'Temperature'
T_array = inT.variables[T_var][:]
nanlocs = np.isnan(T_array[0,:,:])

nanlocs2 = np.isnan(Zgrid)
#---------------------------------------------------------------------

if Plot_Mean_MB == True:
    # load the 'allruns_MB' file:
    out_run_means = 'allrunMBs_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_MB = np.load(out_allruns)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    #plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = np.round(np.linspace(-12,7,12),0))
    plt.contourf(np.flipud(np.nanmean(overall_MB, axis = 0)), cmap = 'RdYlBu', levels = np.linspace(-12,7, 20))
    legend = plt.colorbar()
    legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=20)
    plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = 0)
    plt.title('Kaskawulsh Mass Balance (2007-2018) \n' + str(debristitle) + '\n Temperature Shift = -1.7$\degree$C',fontsize=14)
    plt.savefig(os.path.join(Path2files,'Glacier_MassBalance.png'),bbox_inches = 'tight')
    
    ###extract bf
    mean_mb = np.nanmean(overall_MB, axis = 0)
    np.save(os.path.join(Path2files,'Distributed_MB.npy'),mean_mb)
    mb_std = np.std(overall_MB, axis = 0)
    
    #CALCULATE AAR
    Accumulation_Area = np.where(mean_mb > 0)
    Ablation_Area = np.where(mean_mb <= 0)
    off_glacier = np.where(np.isnan(mean_mb))
    
    cell_A = 0.2*0.2 #area per grid cell km^2 
    glacier_area = (len(Accumulation_Area[0]) + len(Ablation_Area[0])) * cell_A
    acc_area = (len(Accumulation_Area[0]))*cell_A
    AAR = acc_area/glacier_area
    
    n = []
    for i in overall_MB:
        n.append(np.nanmean(i))
        
    out_run0 = 'Net_MB0.npy'
    out_allruns = os.path.join(Path2files,out_run0)
    MB0 = np.load(out_allruns)
    
    out_run1 = 'Net_MB1.npy'
    out_allruns = os.path.join(Path2files,out_run1)
    MB1 = np.load(out_allruns)
    
    out_run2 = 'Net_MB2.npy'
    out_allruns = os.path.join(Path2files,out_run2)
    MB2 = np.load(out_allruns)
    
    for i in range(0,12):
        MB0[i,:,:][nanlocs] = np.nan
        MB1[i,:,:][nanlocs] = np.nan
        MB2[i,:,:][nanlocs] = np.nan
        
    MB0_y = np.nanmean(MB0,axis=1)
    MB0_x = np.nanmean(MB0_y,axis=1)
    
    MB1_y = np.nanmean(MB1,axis=1)
    MB1_x = np.nanmean(MB1_y,axis=1)
    
    MB2_y = np.nanmean(MB2,axis=1)
    MB2_x = np.nanmean(MB2_y,axis=1)
        
    
    print('Mean mass balance plotted from the file: ' + str(out_run_means))
else:
    pass


#------------------------------------------------------------------------------
    
if Plot_NetMelt == True:
    # load the 'allruns_MB' file:
    out_run_means = 'allrunNetMelts_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_NetMelt = np.load(out_allruns)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    #P_array = T_array*100
    plt.contourf(np.flipud(np.nanmean(overall_NetMelt, axis = 0)), cmap = 'YlOrRd', levels = np.linspace(0,12, 20))
    #plt.contourf(np.flipud(np.sum(T_array, axis = 0)), cmap = 'PuRd', levels = np.linspace(0,1, 20))
    legend = plt.colorbar()
    legend.ax.set_ylabel('Melt (m w.e. $a^{-1}$)', rotation=270)
    #legend.ax.set_ylabel('Precipitation (m w.e./yr)', rotation=270,x=1.1)
    plt.title('Kaskawulsh Catchment Ablation (2007-2008) \n' + str(debristitle))
    #plt.title('NARR Downscaled Precipitation (2018)')
    plt.savefig(os.path.join(Path2files,'Glacier_NetMelt.png'))
    #Path2files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Catchment'
    #plt.savefig(os.path.join(Path2files,'Catchment_Precip_2018.png'))
    distributed_melt = np.nanmean(overall_NetMelt,axis=0)
    np.save(os.path.join(Path2files,'Distributed_Melt.npy'),distributed_melt)
    
    print('Mean Melt plotted from the file: ' + str(out_run_means))
else:
    pass

#------------------------------------------------------------------------------
    
if Plot_Accumulation == True:
    # load the 'allruns_MB' file:
    out_run_means = 'allrunAccumulations_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_Accumulation = np.load(out_allruns)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(np.nanmean(overall_Accumulation, axis = 0)), cmap = 'PuRd', levels = np.linspace(0,4, 20))
    #plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
    legend = plt.colorbar()
    legend.ax.set_ylabel('(m w.e. $a^{-1}$)', rotation=270)
    plt.title('Kaskawulsh Catchment Accumulation (2007-2008)')
    plt.savefig(os.path.join(Path2files,'Glacier_Accumulation.png'))
    
    distributed_acc = np.nanmean(overall_Accumulation,axis=0)
    np.save(os.path.join(Path2files,'Distributed_Accumulation.npy'),distributed_acc)
    
    
    print('Mean Melt plotted from the file: ' + str(out_run_means))
else:
    pass

#------------------------------------------------------------------------------
    
if Plot_Rain == True:
    # load the 'allruns_MB' file:
    out_run_means = 'allrunRains_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_Rain = np.load(out_allruns)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(np.nanmean(overall_Rain, axis = 0)), cmap = 'Blues', levels = np.linspace(0,0.4, 20))
    legend = plt.colorbar()
    legend.ax.set_ylabel('Rain (m w.e. $a^{-1}$)', rotation=270)
    plt.title('Rain \n' + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'Glacier_Rain.pdf'))
    
    distributed_rain = np.nanmean(overall_Rain,axis=0)
    np.save(os.path.join(Path2files,'Distributed_Rain.npy'),distributed_rain)
    
    
    print('Mean Melt plotted from the file: ' + str(out_run_means))
else:
    pass

#------------------------------------------------------------------------------
    
if Plot_Refreezing == True:
    # load the 'allruns_MB' file:
    out_run_means = 'allrunRefreezings_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_Refreezing = np.load(out_allruns)
    #R2S
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(np.nanmean(overall_Refreezing, axis = 0)), cmap = 'GnBu', levels = np.linspace(0,0.2, 20), rotation = 45)
    legend = plt.colorbar()
    legend.ax.set_ylabel('Refreezing (m w.e. $a^{-1}$)', rotation=270)
    plt.title('Refreezing')
    plt.savefig(os.path.join(Path2files,'Glacier_Refreezing.pdf'))
    plt.tight_layout()
    
    distributed_refreeze = np.nanmean(overall_Refreezing,axis=0)
    np.save(os.path.join(Path2files,'Distributed_Refreezing.npy'),distributed_refreeze)
    
    print('Mean Refreezing plotted from the file: ' + str(out_run_means))
else:
    pass

        
#------------------------------------------------------------------------------

if Plots_vs_elevation == True:
    print('plot MB, melt,accumulation as a function of elevation')
    
    mean_mb = np.nanmean(overall_MB, axis = 0) #has shape (218,328) --> same as zgrid
    mean_melt = np.nanmean(overall_NetMelt, axis = 0)
    mean_acc = np.nanmean(overall_Accumulation, axis = 0)
    mean_refreeze = np.nanmean(overall_Refreezing, axis = 0)
    mean_rain = np.nanmean(overall_Rain, axis = 0)
    # start two lists to keep track of elevation and mass balance
    massbalance = []
    netmelt = []
    accumulation = []
    elevation = []
    refreezing = []
    rain = []
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            z = Zgrid[x][y]
            rainn = mean_rain[x][y]
            rain.append(rainn)
            elevation.append(z)
            mb = mean_mb[x][y]
            nm = mean_melt[x][y]
            acc = mean_acc[x][y]

            refreeze = mean_refreeze[x][y]
            massbalance.append(mb)
            elevation.append(z)
            netmelt.append(nm)
            accumulation.append(acc)
            refreezing.append(refreeze)
    
    
    plt.figure(figsize=(10,8))
    plt.scatter(massbalance,elevation,c=massbalance, cmap="RdYlBu")
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Mass Balance (m w.e. $a^{-1}$)')
    plt.title('Mass Balance vs Elevation \n' + str(debristitle) + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'MBvELEVATION.pdf'))
    
    plt.figure(figsize=(10,8))
    plt.scatter(netmelt,elevation,c=netmelt, cmap="YlOrRd")
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Melt (m w.e. $a^{-1}$)')
    plt.title('Net Ablation vs Elevation \n' + str(debristitle) + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'NetMeltvELEVATION.pdf'))
    
    plt.figure(figsize=(10,8))
    plt.scatter(accumulation,elevation)
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Accumulation (m w.e. $a^{-1}$)')
    #plt.xlim(0,0.4)
    plt.title('Accumulation vs Elevation \n' + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'AccvELEVATION.pdf'))
    
    plt.figure(figsize=(10,8))
    plt.scatter(refreezing,elevation,c=refreezing, cmap="GnBu")
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Refreezing (m w.e. $a^{-1}$)')
    plt.title('Refreezing \n' + str(debristitle) + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'REFREEZEvELEVATION.pdf'))
    
    plt.figure(figsize=(10,8))
    plt.scatter(rain,elevation)
    plt.xlim(0,0.4)
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $a^{-1}$)')
    plt.title('Rain \n' + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'RAINvELEVATION.pdf'))
    
    
else:
    pass       
    
#----------------------------------------------------------------------------

if Plot_NetMelt_byYear == True:

    def PlotDailyMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        cumulative_melt = np.cumsum(mean_daily_melt)
        
        return [mean_daily_melt,cumulative_melt]
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    
    year = 2007
    ax0.set_title('Net Melt: ' + str(year))
    ax0.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax1.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax1.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax1.legend(loc='upper left', frameon=False)
    
    year = 2008
    ax2.set_title('Net Melt: ' + str(year))
    ax2.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax3.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax3.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax3.legend(loc='upper left', frameon=False)
    
    
    year = 2009
    ax4.set_title('Net Melt: ' + str(year))
    ax4.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax5.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax5.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax5.legend(loc='upper left', frameon=False)
    
    
    year = 2010
    ax6.set_title('Net Melt: ' + str(year))
    ax6.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax7.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax7.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax7.legend(loc='upper left', frameon=False)
    
    
    year = 2011
    ax8.set_title('Net Melt: ' + str(year))
    ax8.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax9.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax9.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax9.legend(loc='upper left', frameon=False)
    
    
    year = 2012
    ax10.set_title('Net Melt: ' + str(year))
    ax10.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax11.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax11.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax11.legend(loc='upper left', frameon=False)
    
    
    year = 2013
    ax12.set_title('Net Melt: ' + str(year))
    ax12.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax13.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax13.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax13.legend(loc='upper left', frameon=False)
    
    
    year = 2014
    ax14.set_title('Net Melt: ' + str(year))
    ax14.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax15.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax15.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax15.legend(loc='upper left', frameon=False)
    
    
    year = 2015
    ax16.set_title('Net Melt: ' + str(year))
    ax16.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax17.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax17.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax17.legend(loc='upper left', frameon=False)
    
    
    year = 2016
    ax18.set_title('Net Melt: ' + str(year))
    ax18.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax19.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax19.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax19.legend(loc='upper left', frameon=False)
    
    
    year = 2017
    ax20.set_title('Net Melt: ' + str(year))
    ax20.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax21.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax21.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax21.legend(loc='upper left', frameon=False)
    
    
    year = 2018
    ax22.set_title('Net Melt: ' + str(year))
    ax22.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax23.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax23.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax23.legend(loc='upper left', frameon=False)
    
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'YearlyNetMelt.pdf'))
else:
    pass
  
#-------------------------------------------------------------------------------------    

if Plot_Refreezing_byYear == True:

    def PlotDailyMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Refreezing_' + str(Case) + '_' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        cumulative_melt = np.cumsum(mean_daily_melt)
        
        return [mean_daily_melt,cumulative_melt]
    
    plt.figure(figsize=(12,12))
    plt.title('Baseline Run \n Debris-free case', y=1.05)
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()
    
    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    
    year = 2007
    ax0.set_title('Total Refreezing: ' + str(year))
    ax0.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax1.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax1.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2008
    ax2.set_title('Total Refreezing: ' + str(year))
    ax2.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax3.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax3.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2009
    ax4.set_title('Total Refreezing: ' + str(year))
    ax4.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax5.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax5.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2010
    ax6.set_title('Total Refreezing: ' + str(year))
    ax6.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax7.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax7.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2011
    ax8.set_title('Total Refreezing: ' + str(year))
    ax8.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax9.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax9.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2012
    ax10.set_title('Total Refreezing: ' + str(year))
    ax10.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax11.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax11.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2013
    ax12.set_title('Total Refreezing: ' + str(year))
    ax12.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax13.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax13.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2014
    ax14.set_title('Total Refreezing: ' + str(year))
    ax14.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax15.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax15.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2015
    ax16.set_title('Total Refreezing: ' + str(year))
    ax16.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax17.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax17.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2016
    ax18.set_title('Total Refreezing: ' + str(year))
    ax18.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax19.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax19.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2017
    ax20.set_title('Total Refreezing: ' + str(year))
    ax20.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax21.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax21.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    year = 2018
    ax22.set_title('Total Refreezing: ' + str(year))
    ax22.plot(PlotDailyMeltbyYear(year)[0],'r')
    ax23.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Total Refreezing (m w.e. $day^{-1}$)', color='r')
    ax23.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'YearlyRefreezing.pdf'))
else:
    pass
        
#-------------------------------------------------------------------------------------    
    
if Plot_Accumulation_byYear == True:

    def PlotDailyMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Accumulation' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        cumulative_melt = np.cumsum(mean_daily_melt)
        
        return [mean_daily_melt,cumulative_melt]
    
    for year in years[1:]:
        print(year)
        array = PlotDailyMeltbyYear(year)[0]
        np.save(os.path.join(Path2files,'DailyAccumulation_GWA_' + str(year) + '.npy'),array)
        
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()
    
    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    
    year = 2007
    ax0.set_title('Accumulation: ' + str(year))
    ax0.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax1.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax1.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax1.legend(loc='upper left', frameon=False)
    
    
    year = 2008
    ax2.set_title('Accumulation: ' + str(year))
    ax2.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax3.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax3.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax3.legend(loc='upper left', frameon=False)
    
    year = 2009
    ax4.set_title('Accumulation: ' + str(year))
    ax4.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax5.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax5.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax5.legend(loc='upper left', frameon=False)
    
    year = 2010
    ax6.set_title('Accumulation: ' + str(year))
    ax6.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax7.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax7.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax7.legend(loc='upper left', frameon=False)
    
    year = 2011
    ax8.set_title('Accumulation: ' + str(year))
    ax8.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax9.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax9.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax9.legend(loc='upper left', frameon=False)
    
    year = 2012
    ax10.set_title('Accumulation: ' + str(year))
    ax10.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax11.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax11.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax11.legend(loc='upper left', frameon=False)
    
    year = 2013
    ax12.set_title('Accumulation: ' + str(year))
    ax12.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax13.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax13.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax13.legend(loc='upper left', frameon=False)
    
    year = 2014
    ax14.set_title('Accumulation: ' + str(year))
    ax14.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax15.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax15.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax15.legend(loc='upper left', frameon=False)
    
    year = 2015
    ax16.set_title('Accumulation: ' + str(year))
    ax16.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax17.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax17.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax17.legend(loc='upper left', frameon=False)
    
    year = 2016
    ax18.set_title('Accumulation: ' + str(year))
    ax18.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax19.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax19.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax19.legend(loc='upper left', frameon=False)
    
    year = 2017
    ax20.set_title('Accumulation: ' + str(year))
    ax20.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax21.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax21.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax21.legend(loc='upper left', frameon=False)
    
    year = 2018
    ax22.set_title('Accumulation: ' + str(year))
    ax22.plot(PlotDailyMeltbyYear(year)[0],'m')
    ax23.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Accumulation \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Accumulation (m w.e. $day^{-1}$)', color='m')
    ax23.set_ylabel('Cumulative Accumulation (m w.e.)', color='k')
    ax23.legend(loc='upper left', frameon=False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'CatchmentAccumulation.png'))
    
    acc_add = 0
    for year in years[1:]:
        #print(year)
        total_acc = PlotDailyMeltbyYear(year)[1][-1]
        #print(total_acc)
        acc_add = acc_add + total_acc
        print(acc_add)
else:
    pass

#-------------------------------------------------------------------------------------    
    
if Plot_Rain_byYear == True:

    def PlotDailyMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        cumulative_melt = np.cumsum(mean_daily_melt)
        
        return [mean_daily_melt,cumulative_melt]
    
    for year in years[1:]:
        print(year)
        array = PlotDailyMeltbyYear(year)[0]
        np.save(os.path.join(Path2files,'DailyRain_GWA_' + str(year) + '.npy'),array)
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()
    
    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    
    plt.suptitle('Non Bias Corrected Rain \n' + str(r2stitle), y=1.05)
    
    year = 2007
    ax0.set_title('Rain: ' + str(year))
    ax0.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax1.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax1.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax1.legend(loc='upper left', frameon=False)
    
    year = 2008
    ax2.set_title('Rain: ' + str(year))
    ax2.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax3.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax3.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax3.legend(loc='upper left', frameon=False)
    
    year = 2009
    ax4.set_title('Rain: ' + str(year))
    ax4.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax5.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax5.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax5.legend(loc='upper left', frameon=False)
    
    year = 2010
    ax6.set_title('Rain: ' + str(year))
    ax6.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax7.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax7.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax7.legend(loc='upper left', frameon=False)
    
    year = 2011
    ax8.set_title('Rain: ' + str(year))
    ax8.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax9.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax9.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax9.legend(loc='upper left', frameon=False)
    
    year = 2012
    ax10.set_title('Rain: ' + str(year))
    ax10.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax11.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax11.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax11.legend(loc='upper left', frameon=False)
    
    year = 2013
    ax12.set_title('Rain: ' + str(year))
    ax12.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax13.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax13.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax13.legend(loc='upper left', frameon=False)
    
    year = 2014
    ax14.set_title('Rain: ' + str(year))
    ax14.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax15.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax15.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax15.legend(loc='upper left', frameon=False)
    
    year = 2015
    ax16.set_title('Rain: ' + str(year))
    ax16.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax17.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax17.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax17.legend(loc='upper left', frameon=False)
    
    year = 2016
    ax18.set_title('Rain: ' + str(year))
    ax18.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax19.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax19.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax19.legend(loc='upper left', frameon=False)
    
    year = 2017
    ax20.set_title('Rain: ' + str(year))
    ax20.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax21.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax21.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax21.legend(loc='upper left', frameon=False)
    
    year = 2018
    ax22.set_title('Rain: ' + str(year))
    ax22.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax23.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Cumulative \n Rain \n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='b')
    ax23.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax23.legend(loc='upper left', frameon=False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'YearlyNonBC_Rain.png'))
else:
    pass
        
#------------------------------------------------------------------------------------- 
    
if Plot_MB_byYear == True:

    def PlotDailyMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_MB' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        cumulative_melt = np.cumsum(mean_daily_melt)
        
        return [mean_daily_melt,cumulative_melt]
    
    m = []
    for year in years[1:]:
        mb = PlotDailyMeltbyYear(year)[1][-1]
        m.append(mb)
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()
    
    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax1.set_ylim([-1.6, 1.85])
    
    year = 2007
    ax0.set_title('Mass Balance: ' + str(year))
    ax0.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax1.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax1.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2008
    ax2.set_title('Mass Balance: ' + str(year))
    ax2.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax3.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax3.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2009
    ax4.set_title('Mass Balance: ' + str(year))
    ax4.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax5.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax5.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2010
    ax6.set_title('Mass Balance: ' + str(year))
    ax6.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax7.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax7.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2011
    ax8.set_title('Mass Balance: ' + str(year))
    ax8.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax9.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax9.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2012
    ax10.set_title('Mass Balance: ' + str(year))
    ax10.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax11.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax11.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2013
    ax12.set_title('Mass Balance: ' + str(year))
    ax12.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax13.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax13.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2014
    ax14.set_title('Mass Balance: ' + str(year))
    ax14.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax15.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax15.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2015
    ax16.set_title('Mass Balance: ' + str(year))
    ax16.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax17.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax17.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2016
    ax18.set_title('Mass Balance: ' + str(year))
    ax18.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax19.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax19.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2017
    ax20.set_title('Mass Balance: ' + str(year))
    ax20.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax21.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax21.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    year = 2018
    ax22.set_title('Mass Balance: ' + str(year))
    ax22.plot(PlotDailyMeltbyYear(year)[0],'purple')
    ax23.plot(PlotDailyMeltbyYear(year)[1],'k')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Mass Balance (m w.e. $day^{-1}$)', color='purple')
    ax23.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'YearlyMassBalanceNEW.pdf'))
else:
    pass
            
#------------------------------------------------------------------------------
if Plot_TotalAblation == True:
    # we assume that: refreezing + net melt = total melt
    print('total ablation = net melt + refreezing')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        filename2 = 'Refreezing_' + str(Case) + '_' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_refreezing_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        mean_daily_refreezing = np.mean(yearly_refreezing_array, axis=0)
        
        daily_total_melt = np.add(mean_daily_melt,mean_daily_refreezing)
        total_melt = np.cumsum(daily_total_melt)

        return [mean_daily_melt,mean_daily_refreezing,total_melt]
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    
    plt.suptitle('Total Ablation \n Debris-free case \n R2S = 1.0$^\circ$C', y=1.05)
    
    year = 2007
    ax0.set_title('Total Melt: ' + str(year))
    ax0.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax0.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax1.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax0.legend(loc='upper left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    
    year = 2008
    ax2.set_title('Total Melt: ' + str(year))
    ax2.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax2.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax3.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax3.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax2.legend(loc='upper left', frameon=False)
    #ax3.legend(loc='upper right', frameon=False)
    
    year = 2009
    ax4.set_title('Total Melt: ' + str(year))
    ax4.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax4.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax5.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax5.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax4.legend(loc='upper left', frameon=False)
    #ax5.legend(loc='center left', frameon=False)
    
    year = 2010
    ax6.set_title('Total Melt: ' + str(year))
    ax6.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax6.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax7.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax7.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax6.legend(loc='upper left', frameon=False)
    #ax7.legend(loc='center left', frameon=False)
    
    year = 2011
    ax8.set_title('Total Melt: ' + str(year))
    ax8.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax8.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax9.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax9.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax8.legend(loc='upper left', frameon=False)
    #ax9.legend(loc='upper right', frameon=False)
    
    year = 2012
    ax10.set_title('Total Melt: ' + str(year))
    ax10.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax10.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax11.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax11.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax10.legend(loc='upper left', frameon=False)
    #ax11.legend(loc='upper right', frameon=False)
    
    year = 2013
    ax12.set_title('Total Melt: ' + str(year))
    ax12.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax12.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax13.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax13.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax12.legend(loc='upper left', frameon=False)
    #ax13.legend(loc='center left', frameon=False)
    
    year = 2014
    ax14.set_title('Total Melt: ' + str(year))
    ax14.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax14.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax15.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax15.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax14.legend(loc='upper left', frameon=False)
    #ax15.legend(loc='upper right', frameon=False)
    
    year = 2015
    ax16.set_title('Total Melt: ' + str(year))
    ax16.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax16.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax17.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax17.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax16.legend(loc='upper left', frameon=False)
    #ax17.legend(loc='upper right', frameon=False)
    
    year = 2016
    ax18.set_title('Total Melt: ' + str(year))
    ax18.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax18.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax19.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax19.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax18.legend(loc='upper left', frameon=False)
    #ax19.legend(loc='center left', frameon=False)
    
    year = 2017
    ax20.set_title('Total Melt: ' + str(year))
    ax20.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax20.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax21.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax21.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax20.legend(loc='upper left', frameon=False)
    #ax21.legend(loc='upper right', frameon=False)
    
    year = 2018
    ax22.set_title('Total Melt: ' + str(year))
    ax22.plot(PlotTotalMeltbyYear(year)[0],'c', label ='Net Melt')
    ax22.plot(PlotTotalMeltbyYear(year)[1], 'r', label = 'Refreezing')
    ax23.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='c')
    ax23.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax22.legend(loc='upper left', frameon=False)
    #ax23.legend(loc='upper right', frameon=False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'YearlyTotalMelt.pdf'))

else:
    pass

#------------------------------------------------------------------------
    
if Plot_Acc_minus_Melt == True:
    # we assume that: refreezing + net melt = total melt
    print('total ablation = net melt + refreezing')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        filename2 = 'Daily_Accumulation' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_refreezing_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0) * (-1)
        mean_daily_refreezing = np.mean(yearly_refreezing_array, axis=0)
        cumulative_abl = np.cumsum(mean_daily_melt)
        cumulative_acc = np.cumsum(mean_daily_refreezing)
        
        daily_total_melt = np.add(mean_daily_refreezing,mean_daily_melt)
        total_melt = np.cumsum(daily_total_melt)
        total_melt[0] = 0

        return [mean_daily_melt,mean_daily_refreezing,total_melt,cumulative_abl,cumulative_acc]
    
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([-0.06,0.06])
    ax1.set_ylim([-1.8, 1.8])
    
    plt.suptitle('Mass Balance \n' + str(debristitle) + '\n Uncorrected Accumulation', y=1.05)
    
    year = 2007
    ax0.set_title('Mass Balance: ' + str(year))
    ax0.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax0.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax0.margins(x=0)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax0.grid(which='both')
    #ax1.legend(loc='upper right', frameon=False)
    ax0.legend(loc='lower left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    
    year = 2008
    ax2.set_title('Mass Balance: ' + str(year))
    ax2.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax2.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax3.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax3.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax2.margins(x=0)
    ax2.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax3.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax2.grid(which='both')
    ax2.legend(loc='lower left', frameon=False)
    ax3.legend(loc='upper right', frameon=False)
    
    year = 2009
    ax4.set_title('Mass Balance: ' + str(year))
    ax4.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax4.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax5.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax5.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax4.margins(x=0)
    ax4.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax5.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax4.grid(which='both')
    ax4.legend(loc='lower left', frameon=False)
    ax5.legend(loc='upper right', frameon=False)
    
    year = 2010
    ax6.set_title('Mass Balance: ' + str(year))
    ax6.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax6.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax7.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax7.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax6.margins(x=0)
    ax6.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax7.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax6.grid(which='both')
    ax6.legend(loc='lower left', frameon=False)
    ax7.legend(loc='upper right', frameon=False)
    
    year = 2011
    ax8.set_title('Mass Balance: ' + str(year))
    ax8.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax8.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax9.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax9.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax8.margins(x=0)
    ax8.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax9.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax8.grid(which='both')
    ax8.legend(loc='lower left', frameon=False)
    ax9.legend(loc='upper right', frameon=False)
    
    year = 2012
    ax10.set_title('Mass Balance: ' + str(year))
    ax10.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax10.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax11.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax11.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax10.margins(x=0)
    ax10.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax11.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax10.grid(which='both')
    ax10.legend(loc='lower left', frameon=False)
    ax11.legend(loc='upper right', frameon=False)
    
    year = 2013
    ax12.set_title('Mass Balance: ' + str(year))
    ax12.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax12.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax13.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax13.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax12.margins(x=0)
    ax12.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax13.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax12.grid(which='both')
    ax12.legend(loc='lower left', frameon=False)
    ax13.legend(loc='upper right', frameon=False)
    
    year = 2014
    ax14.set_title('Mass Balance: ' + str(year))
    ax14.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax14.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax15.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax15.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax14.margins(x=0)
    ax14.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax15.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax14.grid(which='both')
    ax14.legend(loc='lower left', frameon=False)
    ax15.legend(loc='upper right', frameon=False)
    
    year = 2015
    ax16.set_title('Mass Balance: ' + str(year))
    ax16.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax16.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax17.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax17.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax16.margins(x=0)
    ax16.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax17.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax16.grid(which='both')
    ax16.legend(loc='lower left', frameon=False)
    ax17.legend(loc='upper right', frameon=False)
    
    year = 2016
    ax18.set_title('Mass Balance: ' + str(year))
    ax18.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax18.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax19.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax19.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax18.margins(x=0)
    ax18.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax19.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax18.grid(which='both')
    ax18.legend(loc='lower left', frameon=False)
    ax19.legend(loc='upper right', frameon=False)
    
    year = 2017
    ax20.set_title('Mass Balance: ' + str(year))
    ax20.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax20.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax21.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax21.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax20.margins(x=0)
    ax20.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax21.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax20.grid(which='both')
    ax20.legend(loc='lower left', frameon=False)
    ax21.legend(loc='upper right', frameon=False)
    
    year = 2018
    ax22.set_title('Mass Balance: ' + str(year))
    ax22.plot(PlotTotalMeltbyYear(year)[0],'r', label ='Ablation')
    ax22.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Accumulation')
    ax23.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Mass Balance' )
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax23.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax22.margins(x=0)
    ax22.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax23.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax22.grid(which='both')
    ax22.legend(loc='lower left', frameon=False)
    ax23.legend(loc='upper right', frameon=False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'YearlyAcc-Melt.pdf'),bbox_inches='tight')

else:
    pass

#-------------------------------------------------------------------------------
if Calculate_Pmax == True:
    print('Calculate Pmax for: ' + str(Case))
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Accumulation' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_acc = np.load(filepath)
        
        filename2 = 'Refreezing_' + str(Case) + '_' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_refreezing = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_acc = np.mean(yearly_acc, axis = 0)
        mean_daily_refreezing = np.mean(yearly_refreezing, axis=0)
        
        #calculate sum:
        cumulative_acc = np.cumsum(mean_daily_acc)
        cumulative_refreeze = np.cumsum(mean_daily_refreezing)

        return [yearly_acc,yearly_refreezing,cumulative_acc,cumulative_refreeze]
    
    Pmax = []
    for year in years[1:]:
        print(year)
        print('Refreezing = ' + str(PlotTotalMeltbyYear(year)[3][-1]))
        print('Accumulation = ' + str(PlotTotalMeltbyYear(year)[2][-1]))
        # Pmax = Pr / C
        pmax = PlotTotalMeltbyYear(year)[3][-1]/(PlotTotalMeltbyYear(year)[2][-1]) * 100
        print('Pmax = ' + str(pmax))
        print('\n')
        Pmax.append(pmax)
        
    # Plot Pmax as a function of elevation
    mean_acc = np.nanmean(overall_Accumulation, axis = 0)
    mean_refreeze = np.nanmean(overall_Refreezing, axis = 0)
    # start two lists to keep track of elevation and mass balance
    accumulation = []
    refreezing = []
    elevation = []
    Pmaxs = []
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            z = Zgrid[x][y]
            acc = mean_acc[x][y]
            refreeze = mean_refreeze[x][y]
            pmax = (refreeze/acc) * 100

            elevation.append(z)
            accumulation.append(acc)
            refreezing.append(refreeze)
            Pmaxs.append(pmax)
    
    plt.figure(figsize=(10,8))
    plt.scatter(Pmaxs,elevation)
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270)
    plt.ylabel('Elevation (m)')
    plt.xlabel('$P_{max}$ (%)')
    plt.title('Pmax vs Elevation \n Non-Debris Case \n R2S = 1.0$^\circ$C')
    plt.savefig(os.path.join(Path2files,'PMAXvELEVATION.pdf'))
        

else:
    pass    
    
#-------------------------------------------------------------------------------
    
if Plot_FullTimeseries_MB == True:
    # we assume that: refreezing + net melt = total melt
    print('total ablation = net melt + refreezing')
    def PlotMB_fulltimeseries(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        filename2 = 'Daily_Accumulation' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_accumulation_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0) * (-1)
        dailymelt_list = mean_daily_melt.tolist()
        
        mean_daily_refreezing = np.mean(yearly_accumulation_array, axis=0)
        dailyacc_list = mean_daily_refreezing.tolist()
        
        #daily_total_melt = np.add(mean_daily_refreezing,mean_daily_melt)
        #total_melt = np.cumsum(daily_total_melt)
        #total_melt[0] = 0

        return [dailymelt_list,dailyacc_list]
    
    allyears_melt = []
    allyears_acc = []
    for year in years[1:]:
        allyears_melt = allyears_melt + PlotMB_fulltimeseries(year)[0]
        allyears_acc = allyears_acc + PlotMB_fulltimeseries(year)[1]

    daily_mb = np.add(allyears_acc,allyears_melt)
    cumulative_mb = np.cumsum(daily_mb)
    totalmb = round(cumulative_mb[-1],2)
    
    plt.figure(figsize=(13,5))
    plt.title('Kaskawulsh Mass Balance (2007-2018) \n' + str(debristitle) + '\n' + str(r2stitle))
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    ax1.set_ylim([-7.7,10])
    ax0.plot(allyears_melt,'r', label ='Net Melt')
    ax0.plot(allyears_acc, 'b', label = 'Accumulation')
    ax1.plot(cumulative_mb,'k', label = 'Cumulative Mass Balance: ' + str(round(totalmb/12,2)) + ' (m w.e. $yr^{-1}$)')
    ax0.set_xlabel('Year')
    ax0.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    year_starts = [1,365,731,1096,1461,1826,2192,2557,2922,3287,3653,4018]
    year_names = ['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018'] 

    ax0.set_xticks(year_starts)
    ax0.set_xticklabels(year_names)
    ax0.margins(x=0.001)
    ax0.legend(loc='upper left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    plt.savefig(os.path.join(Path2files,'MB_timeseries.pdf'))
    
else:
    pass

#-------------------------------------------------------------------------------
    
if Plot_FullTimeseries_Precip == True:
    # we assume that: refreezing + net melt = total melt
    print('total precip = rain + accumulation (snow)')
    def PlotMB_fulltimeseries(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_rain_array = np.load(filepath)
        
        filename2 = 'Daily_Accumulation' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_accumulation_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        dailyrain_list = mean_daily_rain.tolist()
        
        mean_daily_refreezing = np.mean(yearly_accumulation_array, axis=0)
        dailyacc_list = mean_daily_refreezing.tolist()
        
        #daily_total_melt = np.add(mean_daily_refreezing,mean_daily_melt)
        #total_melt = np.cumsum(daily_total_melt)
        #total_melt[0] = 0

        return [dailyrain_list,dailyacc_list]
    
    allyears_rain = []
    allyears_acc = []
    for year in years[1:]:
        allyears_rain = allyears_rain + PlotMB_fulltimeseries(year)[0]
        allyears_acc = allyears_acc + PlotMB_fulltimeseries(year)[1]

    daily_precip = np.add(allyears_acc,allyears_rain)
    cumulative_precip = np.cumsum(daily_precip)
    cumulative_rain = np.cumsum(allyears_rain)
    cumulative_acc = np.cumsum(allyears_acc)
    
    total_p = round(cumulative_precip[-1],2)
    total_r = round(cumulative_rain[-1],2)
    total_s = round(cumulative_acc[-1],2)
    
    plt.figure(figsize=(13,5))
    plt.title('Kaskawulsh Precipitation (2007-2018)')
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    ax0.set_ylim([0,0.1])
    ax1.set_ylim([0,20])
    ax0.plot(allyears_acc, 'blue', label = 'Daily Snow')
    ax0.plot(allyears_rain,'cyan', label ='Daily Rain')
    ax1.plot(cumulative_precip,'k', label = 'Cumulative Total Precipitation: ' + str(total_p) + ' (m w.e.)')
    ax1.plot(cumulative_rain,'cyan', label = 'Cumulative Rain: ' + str(total_r) + ' (m w.e.)')
    ax1.plot(cumulative_acc,'blue', label = 'Cumulative Snow: ' + str(total_s) + ' (m w.e.)')
    ax0.set_xlabel('Year')
    ax0.set_ylabel('Daily Precipitation (m w.e. $day^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Precipitation (m w.e.)', color='k')
    year_starts = [1,365,731,1096,1461,1826,2192,2557,2922,3287,3653,4018]
    year_names = ['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018'] 

    ax0.set_xticks(year_starts)
    ax0.set_xticklabels(year_names)
    ax0.margins(x=0.001)
    ax0.legend(loc='upper left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    plt.savefig(os.path.join(Path2files,'Precip_timeseries.png'))
    
else:
    pass

#-------------------------------------------------------------------------------
    
if Plot_FullTimeseries_Melt == True:
    # we assume that: refreezing + net melt = total melt
    print('total ablation = net melt + refreezing')
    def PlotMB_fulltimeseries(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        filename2 = 'Refreezing_' + str(Case) + '_' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_accumulation_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        dailymelt_list = mean_daily_melt.tolist()
        
        mean_daily_refreezing = np.mean(yearly_accumulation_array, axis=0)
        dailyrefreeze_list = mean_daily_refreezing.tolist()
        
        #daily_total_melt = np.add(mean_daily_refreezing,mean_daily_melt)
        #total_melt = np.cumsum(daily_total_melt)
        #total_melt[0] = 0

        return [dailymelt_list,dailyrefreeze_list]
    
    allyears_melt = []
    allyears_refreeze = []
    for year in years[1:]:
        allyears_melt = allyears_melt + PlotMB_fulltimeseries(year)[0]
        allyears_refreeze = allyears_refreeze + PlotMB_fulltimeseries(year)[1]

    daily_totalmelt = np.add(allyears_refreeze,allyears_melt)
    cumulative_totalmelt = np.cumsum(daily_totalmelt)
    cumulative_melt = np.cumsum(allyears_melt)
    cumulative_refreeze = np.cumsum(allyears_refreeze)
    
    total_m = round(cumulative_totalmelt[-1],2)
    net_m = round(cumulative_melt[-1],2)
    total_r = round(cumulative_refreeze[-1],2)
    
    plt.figure(figsize=(13,5))
    plt.title('Kaskawulsh Melt (2007-2018)')
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    ax0.set_ylim([0,0.04])
    ax1.set_ylim([0,10])
    ax0.plot(allyears_melt,'r', label ='Net Melt')
    ax0.plot(allyears_refreeze, 'magenta', label = 'Refreezing')
    ax1.plot(cumulative_melt,'r', label ='Net Melt: ' + str(net_m) + ' (m w.e.)')
    ax1.plot(cumulative_refreeze, 'magenta', label = 'Refreezing: ' + str(total_r) + ' (m w.e.)')
    ax1.plot(cumulative_totalmelt,'k', label = 'Cumulative Total Melt: ' + str(total_m) + ' (m w.e.)' )
    ax0.set_xlabel('Year')
    ax0.set_ylabel('Daily Melt (m w.e. $day^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    year_starts = [1,365,731,1096,1461,1826,2192,2557,2922,3287,3653,4018]
    year_names = ['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018'] 

    ax0.set_xticks(year_starts)
    ax0.set_xticklabels(year_names)
    ax0.margins(x=0.001)
    ax0.legend(loc='upper left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    plt.savefig(os.path.join(Path2files,'Melt_timeseries.pdf'))
    
else:
    pass

#------------------------------------------------------------------------

if Plot_NetAblation_IcevsSnow == True:
    # we assume that: refreezing + net melt = total melt
    print('ice melt = net melt - snow melt')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        total_melt = np.cumsum(daily_total_melt)
        total_melt[0] = 0

        return [mean_daily_icemelt,mean_daily_snowmelt,total_melt]
    
    icemelt_cumu = []
    for year in years[1:]:
        cummelt = np.cumsum(PlotTotalMeltbyYear(year)[0])
        melt = cummelt[-1]
        icemelt_cumu.append(melt)
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([0,0.05])
    ax1.set_ylim([0,2.5])
    
    plt.suptitle('Ablation by Melt source \n' + str(debristitle) + '\n' + str(r2stitle), y=1.05)
    
    year = 2007
    PlotTotalMeltbyYear(year)[1][0] = 0
    ax0.set_title('Net Ablation: ' + str(year))
    ax0.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax0.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Melt'  )
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax1.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax0.legend(loc='upper left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    
    year = 2008
    ax2.set_title('Net Ablation: ' + str(year))
    ax2.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax2.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax3.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax3.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax2.legend(loc='upper left', frameon=False)
    #ax3.legend(loc='upper right', frameon=False)
    
    year = 2009
    ax4.set_title('Net Ablation: ' + str(year))
    ax4.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax4.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax5.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax5.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax4.legend(loc='upper left', frameon=False)
    #ax5.legend(loc='center left', frameon=False)
    
    year = 2010
    ax6.set_title('Net Ablation: ' + str(year))
    ax6.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax6.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax7.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax7.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax6.legend(loc='upper left', frameon=False)
    #ax7.legend(loc='center left', frameon=False)
    
    year = 2011
    ax8.set_title('Net Ablation: ' + str(year))
    ax8.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax8.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax9.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax9.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax8.legend(loc='upper left', frameon=False)
    #ax9.legend(loc='upper right', frameon=False)
    
    year = 2012
    ax10.set_title('Net Ablation: ' + str(year))
    ax10.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax10.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax11.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax11.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax10.legend(loc='upper left', frameon=False)
    #ax11.legend(loc='upper right', frameon=False)
    
    year = 2013
    ax12.set_title('Net Ablation: ' + str(year))
    ax12.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax12.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax13.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax13.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax12.legend(loc='upper left', frameon=False)
    #ax13.legend(loc='center left', frameon=False)
    
    year = 2014
    ax14.set_title('Net Ablation: ' + str(year))
    ax14.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax14.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax15.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax15.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax14.legend(loc='upper left', frameon=False)
    #ax15.legend(loc='upper right', frameon=False)
    
    year = 2015
    ax16.set_title('Net Ablation: ' + str(year))
    ax16.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax16.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax17.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax17.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax16.legend(loc='upper left', frameon=False)
    #ax17.legend(loc='upper right', frameon=False)
    
    year = 2016
    ax18.set_title('Net Ablation: ' + str(year))
    ax18.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax18.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax19.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax19.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax18.legend(loc='upper left', frameon=False)
    #ax19.legend(loc='center left', frameon=False)
    
    year = 2017
    ax20.set_title('Net Ablation: ' + str(year))
    ax20.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax20.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax21.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax21.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax20.legend(loc='upper left', frameon=False)
    #ax21.legend(loc='upper right', frameon=False)
    
    year = 2018
    ax22.set_title('Net Ablation: ' + str(year))
    ax22.plot(PlotTotalMeltbyYear(year)[1], 'b', label = 'Snow Melt')
    ax22.plot(PlotTotalMeltbyYear(year)[0],'deeppink', label ='Ice Melt')
    ax23.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n total melt' )
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Ice/Snow Melt (m w.e. $day^{-1}$)', color='b')
    ax23.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    #ax22.legend(loc='upper left', frameon=False)
    #ax23.legend(loc='upper right', frameon=False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'IcevsSnowmelt2.pdf'),bbox_inches='tight')

else:
    pass

#------------------------------------------------------------------------

if Plot_TotalRunoff_Yearly == True:
    # we assume that: refreezing + net melt = total melt
    print('ice melt = net melt - snow melt')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        #filename = 'Daily_IceMelt' + str(year) + '.npy'
        #filepath = os.path.join(Path2files,filename)
        #yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        for i in range(0,len(mean_daily_snowmelt)):
            if mean_daily_snowmelt[i] < 0:
                mean_daily_snowmelt[i] = mean_daily_snowmelt[i-1]
            else:
                pass
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        neglocs = np.where(mean_daily_icemelt < 0) #NEW
        mean_daily_icemelt[neglocs] = 0 #NEW
        mean_daily_snowmelt = mean_daily_netmelt - mean_daily_icemelt
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0
        
        total_snowmelt = np.cumsum(mean_daily_snowmelt)
        total_snowmelt[0] = 0
        total_icemelt = np.cumsum(mean_daily_icemelt)
        total_icemelt[0] = 0
        total_rain = np.cumsum(mean_daily_rain)
        total_rain[0] = 0

        return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain,daily_runoff]
    
    KWArea = 1095554179
    
    plt.figure(figsize=(14,14))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([0,340])
    ax1.set_ylim([0,2.5])
    
    plt.suptitle('Total Runoff by Source \n' + str(debristitle) + '\n Uncorrected Accumulation', y=1.05)
    
    year = 2007
    PlotTotalMeltbyYear(year)[1][0] = 0
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax0.set_title('Net Runoff: ' + str(year))
    ax0.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax0.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax0.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax0.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    #ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Runoff \n =' + str(cumulative_runoff) + ' m w.e.')
    ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax1.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax1.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax1.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax0.margins(x=0)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax0.grid(which='both')
    ax0.set_zorder(ax1.get_zorder()+1)
    ax0.set_frame_on(False)
    ax0.legend(loc='upper left', frameon=False)
    #ax1.legend(loc='upper right', frameon=False)
    
    year = 2008
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax2.set_title('Net Runoff: ' + str(year))
    ax2.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax2.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax2.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax2.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax3.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax3.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax3.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax3.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax3.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax2.margins(x=0)
    ax2.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax3.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax2.grid(which='both')
    ax2.set_zorder(ax3.get_zorder()+1)
    ax2.set_frame_on(False)
    ax2.legend(loc='upper left', frameon=False)
    #ax3.legend(loc='upper right', frameon=False)
    
    year = 2009
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax4.set_title('Net Runoff: ' + str(year))
    ax4.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax4.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax4.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax4.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax5.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax5.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax5.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax5.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax5.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax4.margins(x=0)
    ax4.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax5.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax4.grid(which='both')
    ax4.set_zorder(ax5.get_zorder()+1)
    ax4.set_frame_on(False)
    ax4.legend(loc='upper left', frameon=False)
    #ax5.legend(loc='center left', frameon=False)
    
    year = 2010
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax6.set_title('Net Runoff: ' + str(year))
    ax6.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax6.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax6.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax6.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax7.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax7.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax7.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax7.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax7.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax6.margins(x=0)
    ax6.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax7.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax6.grid(which='both')
    ax6.set_zorder(ax7.get_zorder()+1)
    ax6.set_frame_on(False)
    ax6.legend(loc='upper left', frameon=False)
    #ax7.legend(loc='center left', frameon=False)
    
    year = 2011
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax8.set_title('Net Runoff: ' + str(year))
    ax8.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax8.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax8.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax8.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax9.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax9.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax9.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax9.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax9.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax8.margins(x=0)
    ax8.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax9.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax8.grid(which='both')
    ax8.set_zorder(ax9.get_zorder()+1)
    ax8.set_frame_on(False)
    ax8.legend(loc='upper left', frameon=False)
    #ax9.legend(loc='upper right', frameon=False)
    
    year = 2012
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax10.set_title('Net Runoff: ' + str(year))
    ax10.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax10.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax10.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax10.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax11.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax11.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax11.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax11.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax11.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax10.margins(x=0)
    ax10.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax11.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax10.grid(which='both')
    ax10.set_zorder(ax11.get_zorder()+1)
    ax10.set_frame_on(False)
    ax10.legend(loc='upper left', frameon=False)
    #ax11.legend(loc='upper right', frameon=False)
    
    year = 2013
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax12.set_title('Net Runoff: ' + str(year))
    ax12.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax12.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax12.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax12.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax13.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax13.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax13.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax13.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax13.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax12.margins(x=0)
    ax12.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax13.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax12.grid(which='both')
    ax12.set_zorder(ax13.get_zorder()+1)
    ax12.set_frame_on(False)
    ax12.legend(loc='upper left', frameon=False)
    #ax13.legend(loc='center left', frameon=False)
    
    year = 2014
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax14.set_title('Net Runoff: ' + str(year))
    ax14.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax14.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax14.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax14.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax15.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax15.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax15.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax15.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax15.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax14.margins(x=0)
    ax14.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax15.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax14.grid(which='both')
    ax14.set_zorder(ax15.get_zorder()+1)
    ax14.set_frame_on(False)
    ax14.legend(loc='upper left', frameon=False)
    #ax15.legend(loc='upper right', frameon=False)
    
    year = 2015
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax16.set_title('Net Runoff: ' + str(year))
    ax16.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax16.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax16.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax16.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax17.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax17.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax17.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax17.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax17.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax16.margins(x=0)
    ax16.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax17.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax16.grid(which='both')
    ax16.set_zorder(ax17.get_zorder()+1)
    ax16.set_frame_on(False)
    ax16.legend(loc='upper left', frameon=False)
    #ax17.legend(loc='upper right', frameon=False)
    
    year = 2016
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax18.set_title('Net Runoff: ' + str(year))
    ax18.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax18.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax18.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax18.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax19.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax19.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax19.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax19.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax19.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax18.margins(x=0)
    ax18.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax19.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax18.grid(which='both')
    ax18.set_zorder(ax19.get_zorder()+1)
    ax18.set_frame_on(False)
    ax18.legend(loc='upper left', frameon=False)
    #ax19.legend(loc='center left', frameon=False)
    
    year = 2017
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax20.set_title('Net Runoff: ' + str(year))
    ax20.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax20.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax20.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax20.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax21.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax21.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax21.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax21.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax21.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax20.margins(x=0)
    ax20.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax21.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax20.grid(which='both')
    ax20.set_zorder(ax21.get_zorder()+1)
    ax20.set_frame_on(False)
    ax20.legend(loc='upper left', frameon=False)
    #ax21.legend(loc='upper right', frameon=False)
    
    year = 2018
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax22.set_title('Net Runoff: ' + str(year))
    ax22.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax22.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax22.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax22.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'mediumaquamarine', label ='Rain')
    ax23.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax23.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax23.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax23.plot(PlotTotalMeltbyYear(year)[6],'mediumaquamarine')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax23.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax22.margins(x=0)
    ax22.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax23.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax22.grid(which='both')
    ax22.set_zorder(ax23.get_zorder()+1)
    ax22.set_frame_on(False)
    ax22.legend(loc='upper left', frameon=False)
    #ax23.legend(loc='upper right', frameon=False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'Yearly_Runoff.pdf'),bbox_inches='tight')

else:
    pass

#-------------------------------------------------------------------------------
    
if Plot_FullTimeseries_Runoff == True:
    # we assume that: refreezing + net melt = total melt
    print('total runoff = ice melt + snow melt + rain')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        #Daily_IceMelt
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        rain_list = mean_daily_rain.tolist()
        
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        snowmelt_list = mean_daily_snowmelt.tolist()
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        icemelt_list = mean_daily_icemelt.tolist()
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0

        return [icemelt_list,snowmelt_list,total_melt,rain_list]
    
        
    allyears_ice = []
    allyears_snow = []
    allyears_rain = []
    for year in years[1:]:
        allyears_ice = allyears_ice + PlotTotalMeltbyYear(year)[0]
        allyears_snow = allyears_snow + PlotTotalMeltbyYear(year)[1]
        allyears_rain = allyears_rain + PlotTotalMeltbyYear(year)[3]

    daily_icesnowmelt = np.add(allyears_ice,allyears_snow)
    daily_runoff = np.add(daily_icesnowmelt,allyears_rain)
    cumulative_runoff = np.cumsum(daily_runoff)
    cumul_icesnowmelt = np.cumsum(daily_icesnowmelt)
    totalmb = round(cumulative_runoff[-1],2)
    
    cumulative_icemelt = np.cumsum(allyears_ice)
    total_icemelt = round(cumulative_icemelt[-1],2)
    cumulative_snowmelt = np.cumsum(allyears_snow)
    total_snowmelt = round(cumulative_snowmelt[-1],2)
    cumulative_rain = np.cumsum(allyears_rain)
    total_rain = round(cumulative_rain[-1],2)
    
    plt.figure(figsize=(13,5))
    plt.title('Catchment Total Runoff (2007-2018) \n' + str(debristitle) + '\n' + str(r2stitle))
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    ax0.set_ylim([0,0.05])
    ax1.set_ylim([0,30])
    ax0.plot(allyears_snow, 'b', label = 'Snow Melt')
    ax0.plot(allyears_ice,'deeppink', label ='Ice Melt')
    ax0.plot(allyears_rain, 'gold', label = 'Rain')
    ax1.plot(cumulative_runoff,'k', label = 'Cumulative Runoff: ' + str(totalmb) + ' (m w. e.)')
    ax1.plot(cumulative_icemelt,'deeppink', label = 'Cumulative Ice Melt: ' + str(total_icemelt) + ' (m w. e.)')
    ax1.plot(cumulative_snowmelt,'b', label = 'Cumulative Snow Melt: ' + str(total_snowmelt) + ' (m w. e.)')
    ax1.plot(cumulative_rain,'gold', label = 'Cumulative Rain: ' + str(total_rain) + ' (m w. e.)')
    ax0.set_xlabel('Year')
    ax0.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax1.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    year_starts = [1,365,731,1096,1461,1826,2192,2557,2922,3287,3653,4018]
    year_names = ['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018'] 

    ax0.set_xticks(year_starts)
    ax0.set_xticklabels(year_names)
    ax0.margins(x=0.001)
    ax0.legend(loc='upper left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    plt.savefig(os.path.join(Path2files,'Runoff_timeseries.png'))
    
else:
    pass

if Mean_Runoff_PieChart == True:
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        
        rain_list = mean_daily_rain.tolist()
        
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        for i in range(0,len(mean_daily_snowmelt)):
            if mean_daily_snowmelt[i] < 0:
                mean_daily_snowmelt[i] = mean_daily_snowmelt[i-1]
            else:
                pass
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        neglocs = np.where(mean_daily_icemelt < 0) #NEW
        mean_daily_icemelt[neglocs] = 0 #NEW
        
        icemelt_list = mean_daily_icemelt.tolist()
        
        #ADDED TO FIX PIE CHARTS FOR NO BC ACCUMULATION CASE
        mean_daily_snowmelt = mean_daily_netmelt - mean_daily_icemelt
        mean_daily_snowmelt[0] = 0
        snowmelt_list = mean_daily_snowmelt.tolist()
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0

        return [icemelt_list,snowmelt_list,total_melt,rain_list]
    
        
    allyears_ice = []
    allyears_snow = []
    allyears_rain = []
    for year in years[1:]:
        allyears_ice = allyears_ice + PlotTotalMeltbyYear(year)[0]
        allyears_snow = allyears_snow + PlotTotalMeltbyYear(year)[1]
        allyears_rain = allyears_rain + PlotTotalMeltbyYear(year)[3]

    daily_icesnowmelt = np.add(allyears_ice,allyears_snow)
    daily_runoff = np.add(daily_icesnowmelt,allyears_rain)
    cumulative_runoff = np.cumsum(daily_runoff)
    cumul_icesnowmelt = np.cumsum(daily_icesnowmelt)
    totalmb = round(cumulative_runoff[-1],2)
    
    cumulative_icemelt = np.cumsum(allyears_ice)
    total_icemelt = round(cumulative_icemelt[-1],2)
    cumulative_snowmelt = np.cumsum(allyears_snow)
    total_snowmelt = round(cumulative_snowmelt[-1],2)
    cumulative_rain = np.cumsum(allyears_rain)
    total_rain = round(cumulative_rain[-1],2)
    
    ice_percent = round((total_icemelt/totalmb)*100,2)
    snow_percent = round((total_snowmelt/totalmb)*100,2)
    rain_percent = round((total_rain/totalmb)*100,2)
    
    y = np.array([snow_percent,ice_percent,rain_percent])
    mylabels = ['Snow Melt \n' + str(snow_percent) + '%','Ice Melt \n' + str(ice_percent) + '%','Rain \n' + str(rain_percent) + '%']
    mycolors = ['b','deeppink','mediumaquamarine']

    plt.figure(figsize=(5,5))
    plt.pie(y, labels = mylabels, colors = mycolors,textprops={'fontsize': 15})
    plt.title('Kaskawulsh Runoff (2007-2018) \n' + str(debristitle) + '\n Temperature Shift = -1.7$\degree$C',fontsize=15,y=1.05)
    plt.savefig(os.path.join(Path2files,'Runoff_piechart.png'),bbox_inches='tight')
    
else:
    pass

if Yearly_Runoff_PieChart == True:
    
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        rain_list = mean_daily_rain.tolist()
        cumrain = np.cumsum(rain_list)
        
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        for i in range(0,len(mean_daily_snowmelt)):
            if mean_daily_snowmelt[i] < 0:
                mean_daily_snowmelt[i] = mean_daily_snowmelt[i-1]
            else:
                pass
        
        #snowmelt_list = mean_daily_snowmelt.tolist()
        #cumsnow = np.cumsum(snowmelt_list)
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        neglocs = np.where(mean_daily_icemelt < 0) #NEW
        mean_daily_icemelt[neglocs] = 0 #NEW
        icemelt_list = mean_daily_icemelt.tolist()
        cumice = np.cumsum(icemelt_list)
        
        mean_daily_snowmelt = mean_daily_netmelt - mean_daily_icemelt
        mean_daily_snowmelt[0] = 0
        snowmelt_list = mean_daily_snowmelt.tolist()
        cumsnow = np.cumsum(snowmelt_list)
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0
        
        yr_total = total_melt[-1]
        total_ice = cumice[-1]
        total_snow = cumsnow[-1]
        total_rain = cumrain[-1]
        
        ice_percent = round((total_ice/yr_total)*100,1)
        snow_percent = round((total_snow/yr_total)*100,1)
        rain_percent = round((total_rain/yr_total)*100,1)
        
        y = np.array([snow_percent,ice_percent,rain_percent])
        mylabels = ['Snow Melt \n' + str(snow_percent) + '%','Ice Melt \n' + str(ice_percent) + '%','Rain \n' + str(rain_percent) + '%']

        return [y,mylabels]
    
    mycolors = ['b','deeppink','mediumaquamarine']
    
    plt.figure(figsize=(10,14))
    plt.suptitle('Kaskawulsh Runoff (2007-2018) \n' + str(debristitle) + '\n Uncorrected Accumulation',fontsize=15)
    plt.subplot(4,3,1)
    plt.title('2007',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2007)[0], labels = PlotTotalMeltbyYear(2007)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,2)
    plt.title('2008',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2008)[0], labels = PlotTotalMeltbyYear(2008)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,3)
    plt.title('2009',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2009)[0], labels = PlotTotalMeltbyYear(2009)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,4)
    plt.title('2010',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2010)[0], labels = PlotTotalMeltbyYear(2010)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,5)
    plt.title('2011',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2011)[0], labels = PlotTotalMeltbyYear(2011)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,6)
    plt.title('2012',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2012)[0], labels = PlotTotalMeltbyYear(2012)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,7)
    plt.title('2013',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2013)[0], labels = PlotTotalMeltbyYear(2013)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,8)
    plt.title('2014',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2014)[0], labels = PlotTotalMeltbyYear(2014)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,9)
    plt.title('2015',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2015)[0], labels = PlotTotalMeltbyYear(2015)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,10)
    plt.title('2016',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2016)[0], labels = PlotTotalMeltbyYear(2016)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,11)
    plt.title('2017',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2017)[0], labels = PlotTotalMeltbyYear(2017)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,12)
    plt.title('2018',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2018)[0], labels = PlotTotalMeltbyYear(2018)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.savefig(os.path.join(Path2files,'Yearly_Runoff_piechartV2.pdf'),bbox_inches='tight')
    
    
else:
    pass

if Runoff_Table == True:
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0
        
        total_snowmelt = np.cumsum(mean_daily_snowmelt)
        total_snowmelt[0] = 0
        total_icemelt = np.cumsum(mean_daily_icemelt)
        total_icemelt[0] = 0
        total_rain = np.cumsum(mean_daily_rain)
        total_rain[0] = 0

        return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain]
    
    total_yearly_runoff = ['Total Runoff (m w.e.)']
    total_icemelt_runoff = ['Total Ice Melt (m w.e.)']
    total_snowmelt_runoff = ['Total Snow Melt (m w.e.)']
    total_rain_runoff = ['Total Rain (m w.e.)']
    yrs = []
    #titles = ['Year','Total Runoff','Total Ice Melt','Total Snow Melt','Total Rain']
    for year in years[1:]:
        yrs.append(year)
        total_yearly_runoff.append(round((PlotTotalMeltbyYear(year)[2])[-1],2))
        total_snowmelt_runoff.append(round((PlotTotalMeltbyYear(year)[4])[-1],2))
        total_icemelt_runoff.append(round((PlotTotalMeltbyYear(year)[5])[-1],2))
        total_rain_runoff.append(round((PlotTotalMeltbyYear(year)[6])[-1],2))
    
    
    runoff_data = []
    runoff_data.append(yrs)
    runoff_data.append(total_yearly_runoff)
    runoff_data.append(total_snowmelt_runoff)
    runoff_data.append(total_icemelt_runoff)
    runoff_data.append(total_rain_runoff)
    
    runoff_table = tabulate(runoff_data,headers='firstrow',tablefmt='fancy_grid')
    print(runoff_table)
    
    plt.figure(figsize=(10,8))
    plt.title('Annual Runoff from the Kaskawulsh Glacier (2007-2018) \n' + str(debristitle) + '\n' + str(r2stitle))
    plt.plot(yrs,np.array(total_yearly_runoff[1:]),'k')
    plt.plot(yrs,np.array(total_snowmelt_runoff[1:]),'b')
    plt.plot(yrs,np.array(total_icemelt_runoff[1:]),'deeppink')
    plt.plot(yrs,np.array(total_rain_runoff[1:]),'gold')
    plt.legend(['Total Runoff','Snow Melt','Ice Melt','Rain'])
    plt.ylabel('Annual Runoff (m w.e.)')
    plt.margins(x=0)
    plt.savefig(os.path.join(Path2files,'Runoff_annual_linegraph.pdf'))
    #ax.set_xticklabels(['2007','2008','2009','2010','2011])

else:
    pass

if Avg_Monthly_Runoff == True:
    print('average monthly runoff over entire simulation period')
    def PlotTotalRunoffbyMonth(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        #mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        #mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        #mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        #mean_daily_snowmelt[0] = 0
        
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        
       # mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        #mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        daily_runoff = np.add(mean_daily_netmelt,mean_daily_rain)
        #daily_runoff = mean_daily_icemelt
        jan = daily_runoff[0:31]
        feb = daily_runoff[31:59]
        march = daily_runoff[59:90]
        april = daily_runoff[90:120]
        may = daily_runoff[120:151]
        june = daily_runoff[151:181]
        july = daily_runoff[181:212]
        aug = daily_runoff[212:243]
        sept = daily_runoff[243:273]
        octo = daily_runoff[273:304]
        nov = daily_runoff[304:334]
        dec = daily_runoff[334:365] 
        
        #total_snowmelt = np.cumsum(mean_daily_snowmelt)
        #total_snowmelt[0] = 0
        #total_icemelt = np.cumsum(mean_daily_icemelt)
        #total_icemelt[0] = 0
        #total_rain = np.cumsum(mean_daily_rain)
        #total_rain[0] = 0

        #return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain]
        return [daily_runoff,jan,feb,march,april,may,june,july,aug,sept,octo,nov,dec]
    
    jan = []
    feb = []
    march = []
    april = []
    may = []
    june = []
    july = []
    aug = []
    sept = []
    octo = []
    nov = []
    dec = []
    
    for year in years[1:]:
        jan.append(PlotTotalRunoffbyMonth(year)[1])
        feb.append(PlotTotalRunoffbyMonth(year)[2])
        march.append(PlotTotalRunoffbyMonth(year)[3])
        april.append(PlotTotalRunoffbyMonth(year)[4])
        may.append(PlotTotalRunoffbyMonth(year)[5])
        june.append(PlotTotalRunoffbyMonth(year)[6])
        july.append(PlotTotalRunoffbyMonth(year)[7])
        aug.append(PlotTotalRunoffbyMonth(year)[8])
        sept.append(PlotTotalRunoffbyMonth(year)[9])
        octo.append(PlotTotalRunoffbyMonth(year)[10])
        nov.append(PlotTotalRunoffbyMonth(year)[11])
        dec.append(PlotTotalRunoffbyMonth(year)[12])
    
    jan_array = jan[0]
    feb_array = feb[0]
    march_array = march[0]
    april_array = april[0]
    may_array = may[0]
    june_array = june[0]
    july_array = july[0]
    aug_array = aug[0]
    sept_array = sept[0]
    octo_array = octo[0]
    nov_array = nov[0]
    dec_array = dec[0]
    
    for i in range(1,12):
        jan_array += jan[i]
        feb_array += feb[i]
        march_array += march[i]
        april_array += april[i]
        may_array += may[i]
        june_array += june[i]
        july_array += july[i]
        aug_array += aug[i]
        sept_array += sept[i]
        octo_array += octo[i]
        nov_array += nov[i]
        dec_array += dec[i]
    
    jan_array = jan_array/len(years[1:])
    feb_array = feb_array/len(years[1:])
    march_array = march_array/len(years[1:])
    april_array = april_array/len(years[1:])
    may_array = may_array/len(years[1:])
    june_array = june_array/len(years[1:])
    july_array = july_array/len(years[1:])
    aug_array = aug_array/len(years[1:])
    sept_array = sept_array/len(years[1:])
    octo_array = octo_array/len(years[1:])
    nov_array = nov_array/len(years[1:])
    dec_array = dec_array/len(years[1:])
        
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    #ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    #ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    #ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    #ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    #ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    #ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    #ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    #ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    #ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    #ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    #ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    #ax23 = ax22.twinx()

    #ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([0,0.04])
    #ax1.set_ylim([0,2.7])
    
    plt.suptitle('2007-2018 Average Monthly Runoff \n' + str(debristitle) + '\n' + str(r2stitle), y=1.05)
    
    ax0.set_title('Jan')
    jancumu = np.cumsum(jan_array)
    ax0.plot(jan_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(jancumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax0.margins(x=0.001)
    ax0.set_xlabel('Day of Month')
    ax0.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax0.legend(loc='upper left', frameon=False)
    
    ax2.set_title('Feb')
    febcumu = np.cumsum(feb_array)
    ax2.plot(feb_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(febcumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax2.margins(x=0.001)
    ax2.set_xlabel('Day of Month')
    ax2.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax2.legend(loc='upper left', frameon=False)
    
    ax4.set_title('March')
    marchcumu = np.cumsum(march_array)
    ax4.plot(march_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(marchcumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax4.margins(x=0.001)
    ax4.set_xlabel('Day of Month')
    ax4.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax4.legend(loc='upper left', frameon=False)
    
    ax6.set_title('April')
    aprilcumu = np.cumsum(april_array)
    ax6.plot(april_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(aprilcumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax6.margins(x=0.001)
    ax6.set_xlabel('Day of Month')
    ax6.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax6.legend(loc='upper left', frameon=False)
    
    ax8.set_title('May')
    maycumu = np.cumsum(may_array)
    ax8.plot(may_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(maycumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax8.margins(x=0.001)
    ax8.set_xlabel('Day of Month')
    ax8.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax8.legend(loc='upper left', frameon=False)
    
    ax10.set_title('June')
    junecumu = np.cumsum(june_array)
    ax10.plot(june_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(junecumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax10.margins(x=0.001)
    ax10.set_xlabel('Day of Month')
    ax10.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax10.legend(loc='upper left', frameon=False)
    
    ax12.set_title('July')
    julycumu = np.cumsum(july_array)
    ax12.plot(july_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(julycumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax12.margins(x=0.001)
    ax12.set_xlabel('Day of Month')
    ax12.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax12.legend(loc='upper left', frameon=False)
    
    ax14.set_title('August')
    augcumu = np.cumsum(aug_array)
    ax14.plot(aug_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(augcumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax14.margins(x=0.001)
    ax14.set_xlabel('Day of Month')
    ax14.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax14.legend(loc='upper left', frameon=False)
    
    ax16.set_title('September')
    septcumu = np.cumsum(sept_array)
    ax16.plot(sept_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(septcumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax16.margins(x=0.001)
    ax16.set_xlabel('Day of Month')
    ax16.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax16.legend(loc='upper left', frameon=False)
    
    ax18.set_title('October')
    octcumu = np.cumsum(octo_array)
    ax18.plot(octo_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(octcumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax18.margins(x=0.001)
    ax18.set_xlabel('Day of Month')
    ax18.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax18.legend(loc='upper left', frameon=False)
    
    ax20.set_title('November')
    novcumu = np.cumsum(nov_array)
    ax20.plot(nov_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(novcumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax20.margins(x=0.001)
    ax20.set_xlabel('Day of Month')
    ax20.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax20.legend(loc='upper left', frameon=False)
    
    ax22.set_title('December')
    deccumu = np.cumsum(dec_array)
    ax22.plot(dec_array,'b',label = 'total runoff (2007-2018 avg) \n = '+ str(round(deccumu[-1],2)) + ' m w.e. $month^{-1}$')
    ax22.margins(x=0.001)
    ax22.set_xlabel('Day of Month')
    ax22.set_ylabel('Daily Runoff (m w.e. $day^{-1}$)', color='b')
    ax22.legend(loc='upper left', frameon=False)
    
    #ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Runoff \n =' + str(cumulative_runoff) + ' m w.e.')
    
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'MonthlyRunoff.png'),bbox_inches='tight')

else:
    pass

if Avg_Monthly_Rain == True:
    print('average monthly runoff over entire simulation period')
    def PlotTotalRunoffbyMonth(year):
        # load data file for this year from the extracted outputs
        #filename = 'Daily_Net_Melt' + str(year) + '.npy'
        #filepath = os.path.join(Path2files,filename)
        #yearly_netmelt_array = np.load(filepath)
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        #mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        #mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        #mean_daily_snowmelt[0] = 0
        
        #mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        #daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        #daily_runoff = np.add(mean_daily_netmelt,mean_daily_rain)
        #total_runoff = np.cumsum(daily_runoff)
        #total_melt[0] = 0
        jan = mean_daily_rain[0:31]
        feb = mean_daily_rain[32:59]
        march = mean_daily_rain[60:90]
        april = mean_daily_rain[91:120]
        may = mean_daily_rain[121:151]
        june = mean_daily_rain[152:181]
        july = mean_daily_rain[182:212]
        aug = mean_daily_rain[213:243]
        sept = mean_daily_rain[244:273]
        octo = mean_daily_rain[274:304]
        nov = mean_daily_rain[305:334]
        dec = mean_daily_rain[335:365]
        #total_snowmelt = np.cumsum(mean_daily_snowmelt)
        #total_snowmelt[0] = 0
        #total_icemelt = np.cumsum(mean_daily_icemelt)
        #total_icemelt[0] = 0
        #total_rain = np.cumsum(mean_daily_rain)
        #total_rain[0] = 0

        #return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain]
        return [mean_daily_rain,jan,feb,march,april,may,june,july,aug,sept,octo,nov,dec]
    
    jan = []
    feb = []
    march = []
    april = []
    may = []
    june = []
    july = []
    aug = []
    sept = []
    octo = []
    nov = []
    dec = []
    
    for year in years[1:]:
        jan.append(PlotTotalRunoffbyMonth(year)[1])
        feb.append(PlotTotalRunoffbyMonth(year)[2])
        march.append(PlotTotalRunoffbyMonth(year)[3])
        april.append(PlotTotalRunoffbyMonth(year)[4])
        may.append(PlotTotalRunoffbyMonth(year)[5])
        june.append(PlotTotalRunoffbyMonth(year)[6])
        july.append(PlotTotalRunoffbyMonth(year)[7])
        aug.append(PlotTotalRunoffbyMonth(year)[8])
        sept.append(PlotTotalRunoffbyMonth(year)[9])
        octo.append(PlotTotalRunoffbyMonth(year)[10])
        nov.append(PlotTotalRunoffbyMonth(year)[11])
        dec.append(PlotTotalRunoffbyMonth(year)[12])
    
    jan_array = jan[0]
    feb_array = feb[0]
    march_array = march[0]
    april_array = april[0]
    may_array = may[0]
    june_array = june[0]
    july_array = july[0]
    aug_array = aug[0]
    sept_array = sept[0]
    octo_array = octo[0]
    nov_array = nov[0]
    dec_array = dec[0]
    
    for i in range(1,12):
        jan_array += jan[i]
        feb_array += feb[i]
        march_array += march[i]
        april_array += april[i]
        may_array += may[i]
        june_array += june[i]
        july_array += july[i]
        aug_array += aug[i]
        sept_array += sept[i]
        octo_array += octo[i]
        nov_array += nov[i]
        dec_array += dec[i]
    
    jan_array = jan_array/len(years[1:])
    feb_array = feb_array/len(years[1:])
    march_array = march_array/len(years[1:])
    april_array = april_array/len(years[1:])
    may_array = may_array/len(years[1:])
    june_array = june_array/len(years[1:])
    july_array = july_array/len(years[1:])
    aug_array = aug_array/len(years[1:])
    sept_array = sept_array/len(years[1:])
    octo_array = octo_array/len(years[1:])
    nov_array = nov_array/len(years[1:])
    dec_array = dec_array/len(years[1:])
        
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    #ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    #ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    #ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    #ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    #ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    #ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    #ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    #ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    #ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    #ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    #ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    #ax23 = ax22.twinx()

    #ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([0,0.0025])
    #ax1.set_ylim([0,2.7])
    
    plt.suptitle('Average Monthly Rainfall (2007-2018) \n' + str(debristitle) + '\n' + str(r2stitle), y=1.05)
    
    ax0.set_title('Average Rain: Jan')
    jancumu = np.cumsum(jan_array)
    ax0.plot(jan_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(jancumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax0.margins(x=0.001)
    ax0.set_xlabel('Day of Month')
    ax0.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax0.legend(loc='upper left', frameon=False)
    
    ax2.set_title('Average Rain: Feb')
    febcumu = np.cumsum(feb_array)
    ax2.plot(feb_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(febcumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax2.margins(x=0.001)
    ax2.set_xlabel('Day of Month')
    ax2.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax2.legend(loc='upper left', frameon=False)
    
    ax4.set_title('Average Rain: March')
    marchcumu = np.cumsum(march_array)
    ax4.plot(march_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(marchcumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax4.margins(x=0.001)
    ax4.set_xlabel('Day of Month')
    ax4.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax4.legend(loc='upper left', frameon=False)
    
    ax6.set_title('Average Rain: April')
    aprilcumu = np.cumsum(april_array)
    ax6.plot(april_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(aprilcumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax6.margins(x=0.001)
    ax6.set_xlabel('Day of Month')
    ax6.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax6.legend(loc='upper left', frameon=False)
    
    ax8.set_title('Average Rain: May')
    maycumu = np.cumsum(may_array)
    ax8.plot(may_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(maycumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax8.margins(x=0.001)
    ax8.set_xlabel('Day of Month')
    ax8.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax8.legend(loc='upper left', frameon=False)
    
    ax10.set_title('Average Rain: June')
    junecumu = np.cumsum(june_array)
    ax10.plot(june_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(junecumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax10.margins(x=0.001)
    ax10.set_xlabel('Day of Month')
    ax10.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax10.legend(loc='upper left', frameon=False)
    
    ax12.set_title('Average Rain: July')
    julycumu = np.cumsum(july_array)
    ax12.plot(july_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(julycumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax12.margins(x=0.001)
    ax12.set_xlabel('Day of Month')
    ax12.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax12.legend(loc='upper left', frameon=False)
    
    ax14.set_title('Average Rain: August')
    augcumu = np.cumsum(aug_array)
    ax14.plot(aug_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(augcumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax14.margins(x=0.001)
    ax14.set_xlabel('Day of Month')
    ax14.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax14.legend(loc='upper left', frameon=False)
    
    ax16.set_title('Average Rain: September')
    septcumu = np.cumsum(sept_array)
    ax16.plot(sept_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(septcumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax16.margins(x=0.001)
    ax16.set_xlabel('Day of Month')
    ax16.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax16.legend(loc='upper left', frameon=False)
    
    ax18.set_title('Average Rain: October')
    octcumu = np.cumsum(octo_array)
    ax18.plot(octo_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(octcumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax18.margins(x=0.001)
    ax18.set_xlabel('Day of Month')
    ax18.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax18.legend(loc='upper left', frameon=False)
    
    ax20.set_title('Average Rain: November')
    novcumu = np.cumsum(nov_array)
    ax20.plot(nov_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(novcumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax20.margins(x=0.001)
    ax20.set_xlabel('Day of Month')
    ax20.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax20.legend(loc='upper left', frameon=False)
    
    ax22.set_title('Average Rain: December')
    deccumu = np.cumsum(dec_array)
    ax22.plot(dec_array,'lightseagreen',label = 'average total monthly Rain \n = '+ str(round(deccumu[-1],3)) + ' m w.e. $month^{-1}$')
    ax22.margins(x=0.001)
    ax22.set_xlabel('Day of Month')
    ax22.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='lightseagreen')
    ax22.legend(loc='upper left', frameon=False)
    
    #ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Runoff \n =' + str(cumulative_runoff) + ' m w.e.')
    
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'MonthlyRain.pdf'),bbox_inches='tight')

else:
    pass

#-----------------------------------------------------------------------------

if Plot_IceMelt == True:
    # load the 'allruns_MB' file:
    out_run_means = 'allrunIceMelts_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_MB = np.load(out_allruns)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(np.nanmean(overall_MB, axis = 0)), cmap = 'cool', levels = np.linspace(-1,5, 20))
    legend = plt.colorbar()
    legend.ax.set_ylabel('Ice Melt (m w.e. $a^{-1}$)', rotation=270)
    #plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = np.nan)
    plt.title('Kaskawulsh Ice Melt (2007-2018) \n' + str(debristitle) + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'Glacier_IceMelt.png'),bbox_inches = 'tight')
    
    ###extract bf
    mean_mb = np.nanmean(overall_MB, axis = 0)
    mb_std = np.std(overall_MB, axis = 0)
 
    n = []
    for i in overall_MB:
        n.append(np.nanmean(i))
    
    print('Mean icemelt plotted from the file: ' + str(out_run_means))
else:
    pass

#-----------------------------------------------------------------------------

if Plot_SnowMelt == True:
    # load the 'allruns_MB' file:
    out_run_means = 'allrunSnowMelts_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_MB = np.load(out_allruns)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(np.nanmean(overall_MB, axis = 0)), cmap = 'winter', levels = np.linspace(0,6, 20))
    legend = plt.colorbar()
    legend.ax.set_ylabel('Snow Melt (m w.e. $a^{-1}$)', rotation=270)
    #plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = 0)
    plt.title('Kaskawulsh Snow Melt (2007-2018) \n' + str(debristitle) + '\n' + str(r2stitle))
    plt.savefig(os.path.join(Path2files,'Glacier_SnowMelt.png'),bbox_inches = 'tight')
    
    ###extract bf
    mean_mb = np.nanmean(overall_MB, axis = 0)
    mb_std = np.std(overall_MB, axis = 0)
    
    n = []
    for i in overall_MB:
        n.append(np.nanmean(i))
    
    print('Mean snowmelt plotted from the file: ' + str(out_run_means))
else:
    pass

if Annual_Hydrograph == True:
    print('plotting annual hydrograph')
    def PlotTotalRunoffbyMonth(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        #mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        #mean_daily_snowmelt[0] = 0
        
        #mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        #daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(mean_daily_netmelt,mean_daily_rain)
        #total_runoff = np.cumsum(daily_runoff)
        #total_melt[0] = 0
        jan = np.sum(daily_runoff[0:31])
        feb = np.sum(daily_runoff[32:59])
        march = np.sum(daily_runoff[60:90])
        april = np.sum(daily_runoff[91:120])
        may = np.sum(daily_runoff[121:151])
        june = np.sum(daily_runoff[152:181])
        july = np.sum(daily_runoff[182:212])
        aug = np.sum(daily_runoff[213:243])
        sept = np.sum(daily_runoff[244:273])
        octo = np.sum(daily_runoff[274:304])
        nov = np.sum(daily_runoff[305:334])
        dec = np.sum(daily_runoff[335:365])
        
        year_curve = [jan,feb,march,april,may,june,july,aug,sept,octo,nov,dec]
        #total_snowmelt = np.cumsum(mean_daily_snowmelt)
        #total_snowmelt[0] = 0
        #total_icemelt = np.cumsum(mean_daily_icemelt)
        #total_icemelt[0] = 0
        #total_rain = np.cumsum(mean_daily_rain)
        #total_rain[0] = 0

        #return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain]
        return [year_curve,jan,feb,march,april,may,june,july,aug,sept,octo,nov,dec]
    
    #all_hydrographs = []
    #for year in years[1:]:
     #   all_hydrographs.append(PlotTotalRunoffbyMonth(year)[0])
    
    sum_list = []
    for (a,b,c,d,e,f,g,h,i,j,k,l) in zip(PlotTotalRunoffbyMonth(2007)[0], PlotTotalRunoffbyMonth(2008)[0],PlotTotalRunoffbyMonth(2009)[0],PlotTotalRunoffbyMonth(2010)[0], PlotTotalRunoffbyMonth(2011)[0],PlotTotalRunoffbyMonth(2012)[0],PlotTotalRunoffbyMonth(2013)[0], PlotTotalRunoffbyMonth(2014)[0],PlotTotalRunoffbyMonth(2015)[0],PlotTotalRunoffbyMonth(2016)[0], PlotTotalRunoffbyMonth(2017)[0],PlotTotalRunoffbyMonth(2018)[0]):
        sum_list.append(a+b+c+d+e+f+g+h+i+j+k+l)
    
    mean_list = []
    for i in range(0,len(years[1:])):
        mean_list.append(sum_list[i]/12)
    

    plt.figure(figsize=(10,6))
    plt.title('Kaskawulsh Hydrograph (2007-2018)',fontsize=15)
    plt.plot(PlotTotalRunoffbyMonth(2007)[0],'silver')
    plt.plot(mean_list,'black',linewidth=4)
    plt.plot(PlotTotalRunoffbyMonth(2008)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2009)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2010)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2011)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2012)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2013)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2014)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2015)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2016)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2017)[0],'silver')
    plt.plot(PlotTotalRunoffbyMonth(2018)[0],'silver')
    plt.legend(['single year','mean (2007-2018)'])
    #plt.legend(['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','Mean'])
    #plt.set_xticklabels('jan','feb')
    month_starts = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    month_names = ['Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'] 
    plt.xticks(month_starts,month_names)
    plt.ylabel('Runoff (m w.e. $month^{-1}$)',fontsize=12)
    plt.margins(x=0)
    plt.savefig(os.path.join(Path2files,'Annual_Hydrographs.png'),bbox_inches = 'tight')

else:
    pass

#------------------------------------------------------------------------------

if Plot_Model_outline == True:
    print('plotting model outline from .txt file')
    Sfc_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, sfc_type)
    Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD: 
    #Sfc_kask = np.zeros(Zgrid.shape)
    #Sfc_kask[nanlocs] = np.nan
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(Sfc_grid), cmap = 'bwr', levels = np.linspace(0,1, 5))
    #plt.contourf(np.flipud(Sfc_grid), cmap = 'RdYlBu', levels = np.linspace(0,1,5))
    #plt.contourf(np.flipud(Sfc_kask), cmap = 'bwr', levels = np.linspace(0,1,5))
    legend = plt.colorbar()
    legend.ax.set_ylabel('OFF ice   /   ON ice', rotation=270,labelpad=20)
    #plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = 0)
    plt.title('Kaskawulsh Catchment')
    #plt.savefig(os.path.join(Path2files,'Glacier_Outline.png'),bbox_inches = 'tight')
    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(Zgrid), cmap = 'bone', levels = np.linspace(700,3700,16))
    #plt.contourf(np.flipud(Sfc_grid), cmap = 'RdYlBu', levels = np.linspace(0,1,5))
    #plt.contourf(np.flipud(Sfc_kask), cmap = 'bwr', levels = np.linspace(0,1,5))
    legend = plt.colorbar()
    legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
    #plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = 0)
    #plt.title('Kaskawulsh Catchment')
    plt.savefig(os.path.join(Path2files,'Catchment_Elevation.pdf'),bbox_inches = 'tight')
    
    #icelocs = np.where(Sfc_grid == 0)
    #officelocs = np.where(Sfc_grid == 1)
    #nanlocs = np.isnan(Sfc_grid)
    onice = []
    office = []
    offsite = []
    oniceZ = []
    officeZ = []
    for x in range(len(Xgrid)):
        for y in range(len(Xgrid[0])):
            surface_type = Sfc_grid[x][y]
            elev = Zgrid[x][y]
            if surface_type == 1:
                office.append(surface_type)
                officeZ.append(elev)
            elif surface_type == 0:
                onice.append(surface_type)
                oniceZ.append(elev)
            else:
                offsite.append(surface_type)
                
    
    
else:
    pass

#------------------------------------------------------------------------------

if Mean_Runoff_Timeseries == True:
    print('plotting mean annual runoff')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        
        #filename3 = 'Daily_IceMelt' + str(year) + '.npy'
        #filepath3 = os.path.join(Path2files,filename3)
        #yearly_icemelt_array = np.load(filepath3)
        #Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
        #filename = 'Daily_Rain' + str(year) + '.npy'
        #filepath = os.path.join(Path2_BC_files,filename)
        #yearly_rain_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        #negsnow = np.where(mean_daily_snowmelt < 0)
        for i in range(0,len(mean_daily_snowmelt)):
            if mean_daily_snowmelt[i] < 0:
                mean_daily_snowmelt[i] = mean_daily_snowmelt[i-1]
            else:
                pass
            
        #for i in range(0,len(mean_daily_snowmelt)):
            #if (mean_daily_snowmelt[i] - mean_daily_snowmelt[i-1]) > 200:
                #mean_daily_snowmelt[i] = mean_daily_snowmelt[i-1]
            #else:
                #pass
        
        #mean_daily_snowmelt[0] = 0
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        neglocs = np.where(mean_daily_icemelt < 0) #NEW
        mean_daily_icemelt[neglocs] = 0 #NEW
        
        mean_daily_snowmelt = mean_daily_netmelt - mean_daily_icemelt
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0
        
        total_snowmelt = np.cumsum(mean_daily_snowmelt)
        total_snowmelt[0] = 0
        total_icemelt = np.cumsum(mean_daily_icemelt)
        total_icemelt[0] = 0
        total_rain = np.cumsum(mean_daily_rain)
        total_rain[0] = 0

        return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain]
    
    allyrs_ice = []
    allyrs_snow = []
    allyrs_rain = []
    for year in years[1:]:
        allyrs_ice.append(PlotTotalMeltbyYear(year)[0][0:365])
        allyrs_snow.append(PlotTotalMeltbyYear(year)[1][0:365])
        allyrs_rain.append(PlotTotalMeltbyYear(year)[3][0:365])
        
    ice_sum = allyrs_ice[0]
    snow_sum = allyrs_snow[0]
    rain_sum = allyrs_rain[0]
    for i in range(1,len(years[1:])):
        ice_sum += allyrs_ice[i]
        snow_sum += allyrs_snow[i]
        rain_sum += allyrs_rain[i]
    
    KWArea = 1095554179
    
    ice_avg = ice_sum/len(years[1:])
    ice_avg[0] = 0
    snow_avg = snow_sum/len(years[1:])
    snow_avg[0]= 0
    rain_avg = rain_sum/len(years[1:])
    
    ice_discharge = ice_avg*(KWArea/86400) #units = m3/s
    snow_discharge = snow_avg*(KWArea/86400) #units = m3/s
    rain_discharge = rain_avg*(KWArea/86400) #units = m3/s
    
    
    total_runoff = (np.array(ice_avg)) + (np.array(snow_avg)) + (np.array(rain_avg))
    cumulative_runoff = np.nancumsum(total_runoff)
    cumulative_snow = np.cumsum(snow_avg)
    cumulative_ice = np.cumsum(ice_avg)
    cumulative_rain = np.cumsum(rain_avg)
    
    #total_discharge = total_runoff*(KWArea/86400) #units = m3/s
    total_discharge = ice_discharge + snow_discharge + rain_discharge
    
    plt.figure(figsize=(8,8))
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    plt.title('Runoff (2007-2018) \n Debris-Free Case \n Temperature Shift = -1.7$\degree$C',fontsize=15)
    ax1.plot(cumulative_runoff,'black',zorder=5)
    ax1.plot(cumulative_snow,'b',zorder=3)
    ax1.plot(cumulative_ice,'deeppink',zorder=2)
    ax1.plot(cumulative_rain,'mediumaquamarine',zorder=1)
    #ax0.plot(snow_avg,'b')
    #ax0.fill_between(dayofyear,snow_avg+0.0045,snow_avg-0.0045)
    #ax0.plot(ice_avg,'deeppink',zorder=5)
    #ax0.plot(rain_avg,'mediumaquamarine',zorder=6)
    #ax0.plot(total_runoff,'black',linewidth=2,zorder=7)
    ax0.plot(total_discharge,'black',linewidth=2)
    ax0.plot(snow_discharge,'b')
    ax0.plot(ice_discharge,'deeppink')
    ax0.plot(rain_discharge,'mediumaquamarine')
    plt.legend(['Total Runoff: ' + str(np.round(cumulative_runoff[-1],2)) + ' m w.e.','Snow Melt: ' + str(np.round(cumulative_snow[-1],2)) + ' m w.e.','Ice Melt: ' + str(np.round(cumulative_ice[-1],2)) + ' m w.e.','Rain: ' + str(np.round(cumulative_rain[-1],2)) + ' m w.e.'],fontsize=15)
    ax1.set_ylabel('Cumulative Runoff (m w.e. year$^{-1}$)', color='k',fontsize=15)
    ax0.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k',fontsize=15)
    ax0.set_xlabel('Day Of Year',fontsize=15)
    #plt.legend(['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','Mean'])
    #plt.set_xticklabels('jan','feb')
    #month_starts = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    #month_names = ['Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'] 
    plt.margins(x=0)
    ax0.margins(y=0)
    ax1.margins(y=0)
    ax0.set_ylim(0,400)
    ax1.set_ylim(0,2.5)
    plt.setp(ax0.get_xticklabels(), fontsize=15)
    plt.setp(ax0.get_yticklabels(), fontsize=15)
    plt.setp(ax1.get_yticklabels(), fontsize=15)
    #ax0.set_ylim(0,0.03) works for runoff in mw.e./day
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax0.grid(which='both')
    #ax1.set_yticks(0,1.9)
    ax0.set_zorder(ax1.get_zorder()+1)
    ax0.set_frame_on(False)
    plt.savefig(os.path.join(Path2files,'Mean_Runoff_timeseries.png'),bbox_inches = 'tight')
    
else:
    pass
#--------------------------------------------------------------------------------------
if Plot_NARR_distributed_prcp == True:
    print('NARR downscaled Precip')
    out_run_means = 'allrunAccumulations_' + Case + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_Accumulation = np.load(out_allruns)
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(overall_Accumulation), cmap = 'BuPu', levels = np.linspace(0.4,0.8, 21))
    #plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
    legend = plt.colorbar()
    legend.ax.tick_params(labelsize=15) 
    legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,labelpad=20,fontsize=15)
    #legend.ax.set_ylabel('Bed Elevation (m a.s.l)',labelpad=15, rotation=270,fontsize=15)
    plt.title('Uncorrected Precipitation (2007-2008)',fontsize=15)
    plt.savefig(os.path.join(Path2files,'Glacier_Mean_Precipitation.png'))
    
    print('Mean Melt plotted from the file: ' + str(out_run_means))
else:
    pass

#------------------------------------------------------------------------------
if Plot_Monthly_Distributed_Rain == True:
    print('monthly distributed rain')
    janpath = os.path.join(Path2files,'jan_distributedRain.npy')
    febpath = os.path.join(Path2files,'feb_distributedRain.npy')
    marchpath = os.path.join(Path2files,'march_distributedRain.npy')
    aprilpath = os.path.join(Path2files,'april_distributedRain.npy')
    maypath = os.path.join(Path2files,'may_distributedRain.npy')
    junepath = os.path.join(Path2files,'june_distributedRain.npy')
    julypath = os.path.join(Path2files,'july_distributedRain.npy')
    augpath = os.path.join(Path2files,'aug_distributedRain.npy')
    septpath = os.path.join(Path2files,'sept_distributedRain.npy')
    octpath = os.path.join(Path2files,'oct_distributedRain.npy')
    novpath = os.path.join(Path2files,'nov_distributedRain.npy')
    decpath = os.path.join(Path2files,'dec_distributedRain.npy')

    janfile = np.load(janpath)
    febfile = np.load(febpath)
    marchfile = np.load(marchpath)
    aprilfile = np.load(aprilpath)
    mayfile = np.load(maypath)
    junefile = np.load(junepath)
    julyfile = np.load(julypath)
    augfile = np.load(augpath)
    septfile = np.load(septpath)
    octfile = np.load(octpath)
    novfile = np.load(novpath)
    decfile = np.load(decpath)
    
    fig = plt.figure(figsize=(18,15))
    plt.suptitle('Distributed Monthly Rain (Downscaled NARR) \n' + str(r2stitle),y=0.94)
    #fig.colorbar(m, label='colorbar with length <1', ticks=0.1, loc='r', length=0.7)
    plt.subplot(4,3,1)
    CS = plt.contourf(np.flipud(janfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30)) 
    plt.title('January')
    plt.subplot(4,3,2)
    CS = plt.contourf(np.flipud(febfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('February')
    plt.subplot(4,3,3)
    CS = plt.contourf(np.flipud(marchfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('March')
    legend = plt.colorbar()
    legend.ax.set_ylabel('Precipitation (m w.e. $day^{-1}$)', rotation=270)
    plt.subplot(4,3,4)
    CS = plt.contourf(np.flipud(aprilfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('April')
    plt.subplot(4,3,5)
    CS = plt.contourf(np.flipud(mayfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('May')
    plt.subplot(4,3,6)
    CS = plt.contourf(np.flipud(junefile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('June')
    plt.subplot(4,3,7)
    CS = plt.contourf(np.flipud(julyfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('July')
    plt.subplot(4,3,8)
    CS = plt.contourf(np.flipud(augfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('August')
    plt.subplot(4,3,9)
    CS = plt.contourf(np.flipud(septfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('September')
    plt.subplot(4,3,10)
    CS = plt.contourf(np.flipud(octfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('October')
    plt.subplot(4,3,11)
    CS = plt.contourf(np.flipud(novfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('November')
    plt.subplot(4,3,12)
    CS = plt.contourf(np.flipud(decfile), cmap = 'Blues', levels = np.linspace(0,0.0025, 30))   
    plt.title('December')
    #fig.colorbar(cmap = 'Blues')
    #plt.show()
    plt.savefig(os.path.join(Path2files,'Monthly_Distributed_Rain.pdf'))

else:
    pass


#------------------------------------------------------------------------------
if Plot_Monthly_Distributed_Temp == True:
    print('monthly distributed temp')
    janpath = os.path.join(Path2files,'jan_distributedTemp.npy')
    febpath = os.path.join(Path2files,'feb_distributedTemp.npy')
    marchpath = os.path.join(Path2files,'march_distributedTemp.npy')
    aprilpath = os.path.join(Path2files,'april_distributedTemp.npy')
    maypath = os.path.join(Path2files,'may_distributedTemp.npy')
    junepath = os.path.join(Path2files,'june_distributedTemp.npy')
    julypath = os.path.join(Path2files,'july_distributedTemp.npy')
    augpath = os.path.join(Path2files,'aug_distributedTemp.npy')
    septpath = os.path.join(Path2files,'sept_distributedTemp.npy')
    octpath = os.path.join(Path2files,'oct_distributedTemp.npy')
    novpath = os.path.join(Path2files,'nov_distributedTemp.npy')
    decpath = os.path.join(Path2files,'dec_distributedTemp.npy')

    janfile = np.load(janpath)
    febfile = np.load(febpath)
    marchfile = np.load(marchpath)
    aprilfile = np.load(aprilpath)
    mayfile = np.load(maypath)
    junefile = np.load(junepath)
    julyfile = np.load(julypath)
    augfile = np.load(augpath)
    septfile = np.load(septpath)
    octfile = np.load(octpath)
    novfile = np.load(novpath)
    decfile = np.load(decpath)
    
    fig = plt.figure(figsize=(18,15))
    plt.suptitle('Distributed Monthly Temperature',y=0.93)
    #fig.colorbar(m, label='colorbar with length <1', ticks=0.1, loc='r', length=0.7)
    plt.subplot(4,3,1)
    CS = plt.contourf(np.flipud(janfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30)) 
    plt.title('January')
    plt.subplot(4,3,2)
    CS = plt.contourf(np.flipud(febfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))     
    plt.title('February')
    plt.subplot(4,3,3)
    CS = plt.contourf(np.flipud(marchfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('March')
    legend = plt.colorbar()
    legend.ax.set_ylabel('Temperature (C)', rotation=270)
    plt.subplot(4,3,4)
    CS = plt.contourf(np.flipud(aprilfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('April')
    plt.subplot(4,3,5)
    CS = plt.contourf(np.flipud(mayfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('May')
    plt.subplot(4,3,6)
    CS = plt.contourf(np.flipud(junefile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('June')
    legend = plt.colorbar()
    legend.ax.set_ylabel('Temperature (C)', rotation=270)
    plt.subplot(4,3,7)
    CS = plt.contourf(np.flipud(julyfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('July')
    plt.subplot(4,3,8)
    CS = plt.contourf(np.flipud(augfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('August')
    plt.subplot(4,3,9)
    CS = plt.contourf(np.flipud(septfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('September')
    legend = plt.colorbar()
    legend.ax.set_ylabel('Temperature (C)', rotation=270)
    plt.subplot(4,3,10)
    CS = plt.contourf(np.flipud(octfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('October')
    plt.subplot(4,3,11)
    CS = plt.contourf(np.flipud(novfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('November')
    plt.subplot(4,3,12)
    CS = plt.contourf(np.flipud(decfile), cmap = 'YlOrRd', levels = np.linspace(-20,15, 30))   
    plt.title('December')
    legend = plt.colorbar()
    legend.ax.set_ylabel('Temperature (C)', rotation=270)
    #fig.colorbar(cmap = 'Blues')
    #plt.show()
    plt.savefig(os.path.join(Path2files,'Monthly_Distributed_Temp.pdf'))

else:
    pass

#---------------------------------------------------------------------------------

if Plot_Monthly_RainvZ == True:
    print('Rain vs Z')
    janpath = os.path.join(Path2files,'jan_distributedRain.npy')
    febpath = os.path.join(Path2files,'feb_distributedRain.npy')
    marchpath = os.path.join(Path2files,'march_distributedRain.npy')
    aprilpath = os.path.join(Path2files,'april_distributedRain.npy')
    maypath = os.path.join(Path2files,'may_distributedRain.npy')
    junepath = os.path.join(Path2files,'june_distributedRain.npy')
    julypath = os.path.join(Path2files,'july_distributedRain.npy')
    augpath = os.path.join(Path2files,'aug_distributedRain.npy')
    septpath = os.path.join(Path2files,'sept_distributedRain.npy')
    octpath = os.path.join(Path2files,'oct_distributedRain.npy')
    novpath = os.path.join(Path2files,'nov_distributedRain.npy')
    decpath = os.path.join(Path2files,'dec_distributedRain.npy')

    janfile = np.load(janpath)
    febfile = np.load(febpath)
    marchfile = np.load(marchpath)
    aprilfile = np.load(aprilpath)
    mayfile = np.load(maypath)
    junefile = np.load(junepath)
    julyfile = np.load(julypath)
    augfile = np.load(augpath)
    septfile = np.load(septpath)
    octfile = np.load(octpath)
    novfile = np.load(novpath)
    decfile = np.load(decpath)
    
    janrain = []
    febrain = []
    marchrain = [] 
    aprilrain = [] 
    mayrain = [] 
    junerain = []
    julyrain = [] 
    augrain = [] 
    septrain = [] 
    octrain = [] 
    novrain = [] 
    decrain = []
    elevation= []
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            elevation.append(Zgrid[x][y])
            janrain.append(janfile[x][y])
            febrain.append(febfile[x][y])
            marchrain.append(marchfile[x][y]) 
            aprilrain.append(aprilfile[x][y]) 
            mayrain.append(mayfile[x][y]) 
            junerain.append(junefile[x][y])
            julyrain.append(julyfile[x][y]) 
            augrain.append(augfile[x][y]) 
            septrain.append(septfile[x][y]) 
            octrain.append(octfile[x][y]) 
            novrain.append(novfile[x][y]) 
            decrain.append(decfile[x][y])
    
    fig = plt.figure(figsize=(15,15))
    plt.suptitle('NARR Rain vs Elevation',y=1.02)
    #fig.colorbar(m, label='colorbar with length <1', ticks=0.1, loc='r', length=0.7)
    plt.subplot(4,3,1)
    plt.scatter(janrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('January')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,2)
    plt.scatter(febrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('February')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,3)
    plt.scatter(marchrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('March')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,4)
    plt.scatter(aprilrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('April')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,5)
    plt.scatter(mayrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('May')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,6)
    plt.scatter(junerain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('June')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,7)
    plt.scatter(julyrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('July')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,8)
    plt.scatter(augrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('August')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,9)
    plt.scatter(septrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)') 
    plt.title('September')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,10)
    plt.scatter(octrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('October')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,11)
    plt.scatter(novrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('November')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,12)
    plt.scatter(decrain,elevation)
    plt.xlim(0,0.002)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Rain (m w.e. $day^{-1}$)')
    plt.title('December')
    plt.locator_params(axis="x", nbins=5)
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'Monthly_RainvZ.png'))

else:
    pass

#-------------------------------------------------------------------------------
if Plot_Monthly_TempvZ == True:
    print('Temp vs Z')
    janpath = os.path.join(Path2files,'jan_distributedTemp.npy')
    febpath = os.path.join(Path2files,'feb_distributedTemp.npy')
    marchpath = os.path.join(Path2files,'march_distributedTemp.npy')
    aprilpath = os.path.join(Path2files,'april_distributedTemp.npy')
    maypath = os.path.join(Path2files,'may_distributedTemp.npy')
    junepath = os.path.join(Path2files,'june_distributedTemp.npy')
    julypath = os.path.join(Path2files,'july_distributedTemp.npy')
    augpath = os.path.join(Path2files,'aug_distributedTemp.npy')
    septpath = os.path.join(Path2files,'sept_distributedTemp.npy')
    octpath = os.path.join(Path2files,'oct_distributedTemp.npy')
    novpath = os.path.join(Path2files,'nov_distributedTemp.npy')
    decpath = os.path.join(Path2files,'dec_distributedTemp.npy')

    janfile = np.load(janpath)
    febfile = np.load(febpath)
    marchfile = np.load(marchpath)
    aprilfile = np.load(aprilpath)
    mayfile = np.load(maypath)
    junefile = np.load(junepath)
    julyfile = np.load(julypath)
    augfile = np.load(augpath)
    septfile = np.load(septpath)
    octfile = np.load(octpath)
    novfile = np.load(novpath)
    decfile = np.load(decpath)
    
    jantemp = []
    febtemp = []
    marchtemp = [] 
    apriltemp = [] 
    maytemp = [] 
    junetemp = []
    julytemp = [] 
    augtemp = [] 
    septtemp = [] 
    octtemp = [] 
    novtemp = [] 
    dectemp = []
    elevation= []
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            elevation.append(Zgrid[x][y])
            jantemp.append(janfile[x][y])
            febtemp.append(febfile[x][y])
            marchtemp.append(marchfile[x][y]) 
            apriltemp.append(aprilfile[x][y]) 
            maytemp.append(mayfile[x][y]) 
            junetemp.append(junefile[x][y])
            julytemp.append(julyfile[x][y]) 
            augtemp.append(augfile[x][y]) 
            septtemp.append(septfile[x][y]) 
            octtemp.append(octfile[x][y]) 
            novtemp.append(novfile[x][y]) 
            dectemp.append(decfile[x][y])
    
    fig = plt.figure(figsize=(15,15))
    plt.suptitle('NARR Temperature vs Elevation',y=1.02)
    #fig.colorbar(m, label='colorbar with length <1', ticks=0.1, loc='r', length=0.7)
    plt.subplot(4,3,1)
    plt.scatter(jantemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('January')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,2)
    plt.scatter(febtemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('February')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,3)
    plt.scatter(marchtemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('March')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,4)
    plt.scatter(apriltemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('April')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,5)
    plt.scatter(maytemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('May')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,6)
    plt.scatter(junetemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('June')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,7)
    plt.scatter(julytemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('July')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,8)
    plt.scatter(augtemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('August')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,9)
    plt.scatter(septtemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)') 
    plt.title('September')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,10)
    plt.scatter(octtemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('October')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,11)
    plt.scatter(novtemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('November')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,12)
    plt.scatter(dectemp,elevation)
    plt.xlim(-25,15)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Temperature (C)')
    plt.title('December')
    plt.locator_params(axis="x", nbins=5)
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'Monthly_TempvZ.pdf'))

else:
    pass

#---------------------------------------------------------------------------------

if Plot_Monthly_NARRPrecip_vsZ == True:
    print('Rain vs Z')
    janpath = os.path.join(Path2files,'jan_distributedPrecip.npy')
    febpath = os.path.join(Path2files,'feb_distributedPrecip.npy')
    marchpath = os.path.join(Path2files,'march_distributedPrecip.npy')
    aprilpath = os.path.join(Path2files,'april_distributedPrecip.npy')
    maypath = os.path.join(Path2files,'may_distributedPrecip.npy')
    junepath = os.path.join(Path2files,'june_distributedPrecip.npy')
    julypath = os.path.join(Path2files,'july_distributedPrecip.npy')
    augpath = os.path.join(Path2files,'aug_distributedPrecip.npy')
    septpath = os.path.join(Path2files,'sept_distributedPrecip.npy')
    octpath = os.path.join(Path2files,'oct_distributedPrecip.npy')
    novpath = os.path.join(Path2files,'nov_distributedPrecip.npy')
    decpath = os.path.join(Path2files,'dec_distributedPrecip.npy')

    janfile = np.load(janpath)
    febfile = np.load(febpath)
    marchfile = np.load(marchpath)
    aprilfile = np.load(aprilpath)
    mayfile = np.load(maypath)
    junefile = np.load(junepath)
    julyfile = np.load(julypath)
    augfile = np.load(augpath)
    septfile = np.load(septpath)
    octfile = np.load(octpath)
    novfile = np.load(novpath)
    decfile = np.load(decpath)
    
    janrain = []
    febrain = []
    marchrain = [] 
    aprilrain = [] 
    mayrain = [] 
    junerain = []
    julyrain = [] 
    augrain = [] 
    septrain = [] 
    octrain = [] 
    novrain = [] 
    decrain = []
    elevation= []
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            elevation.append(Zgrid[x][y])
            janrain.append(janfile[x][y])
            febrain.append(febfile[x][y])
            marchrain.append(marchfile[x][y]) 
            aprilrain.append(aprilfile[x][y]) 
            mayrain.append(mayfile[x][y]) 
            junerain.append(junefile[x][y])
            julyrain.append(julyfile[x][y]) 
            augrain.append(augfile[x][y]) 
            septrain.append(septfile[x][y]) 
            octrain.append(octfile[x][y]) 
            novrain.append(novfile[x][y]) 
            decrain.append(decfile[x][y])
    
    fig = plt.figure(figsize=(15,15))
    plt.suptitle('NARR Precip vs Elevation',y=1.02)
    #fig.colorbar(m, label='colorbar with length <1', ticks=0.1, loc='r', length=0.7)
    plt.subplot(4,3,1)
    plt.scatter(janrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('January')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,2)
    plt.scatter(febrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('February')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,3)
    plt.scatter(marchrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('March')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,4)
    plt.scatter(aprilrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('April')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,5)
    plt.scatter(mayrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('May')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,6)
    plt.scatter(junerain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('June')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,7)
    plt.scatter(julyrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('July')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,8)
    plt.scatter(augrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('August')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,9)
    plt.scatter(septrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)') 
    plt.title('September')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,10)
    plt.scatter(octrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('October')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,11)
    plt.scatter(novrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('November')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,12)
    plt.scatter(decrain,elevation)
    plt.xlim(0.001,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Precipitation (m w.e. $day^{-1}$)')
    plt.title('December')
    plt.locator_params(axis="x", nbins=5)
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'Monthly_PrecipitationvZ.pdf'))

else:
    pass

#---------------------------------------------------------------------------------

if Plot_Monthly_SnowvZ == True:
    print('Rain vs Z')
    janpath = os.path.join(Path2files,'jan_distributedSnow.npy')
    febpath = os.path.join(Path2files,'feb_distributedSnow.npy')
    marchpath = os.path.join(Path2files,'march_distributedSnow.npy')
    aprilpath = os.path.join(Path2files,'april_distributedSnow.npy')
    maypath = os.path.join(Path2files,'may_distributedSnow.npy')
    junepath = os.path.join(Path2files,'june_distributedSnow.npy')
    julypath = os.path.join(Path2files,'july_distributedSnow.npy')
    augpath = os.path.join(Path2files,'aug_distributedSnow.npy')
    septpath = os.path.join(Path2files,'sept_distributedSnow.npy')
    octpath = os.path.join(Path2files,'oct_distributedSnow.npy')
    novpath = os.path.join(Path2files,'nov_distributedSnow.npy')
    decpath = os.path.join(Path2files,'dec_distributedSnow.npy')

    janfile = np.load(janpath)
    febfile = np.load(febpath)
    marchfile = np.load(marchpath)
    aprilfile = np.load(aprilpath)
    mayfile = np.load(maypath)
    junefile = np.load(junepath)
    julyfile = np.load(julypath)
    augfile = np.load(augpath)
    septfile = np.load(septpath)
    octfile = np.load(octpath)
    novfile = np.load(novpath)
    decfile = np.load(decpath)
    
    janrain = []
    febrain = []
    marchrain = [] 
    aprilrain = [] 
    mayrain = [] 
    junerain = []
    julyrain = [] 
    augrain = [] 
    septrain = [] 
    octrain = [] 
    novrain = [] 
    decrain = []
    elevation= []
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            elevation.append(Zgrid[x][y])
            janrain.append(janfile[x][y])
            febrain.append(febfile[x][y])
            marchrain.append(marchfile[x][y]) 
            aprilrain.append(aprilfile[x][y]) 
            mayrain.append(mayfile[x][y]) 
            junerain.append(junefile[x][y])
            julyrain.append(julyfile[x][y]) 
            augrain.append(augfile[x][y]) 
            septrain.append(septfile[x][y]) 
            octrain.append(octfile[x][y]) 
            novrain.append(novfile[x][y]) 
            decrain.append(decfile[x][y])
    
    fig = plt.figure(figsize=(15,15))
    plt.suptitle('NARR Snow vs Elevation',y=1.02)
    #fig.colorbar(m, label='colorbar with length <1', ticks=0.1, loc='r', length=0.7)
    plt.subplot(4,3,1)
    plt.scatter(janrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('January')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,2)
    plt.scatter(febrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('February')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,3)
    plt.scatter(marchrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('March')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,4)
    plt.scatter(aprilrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('April')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,5)
    plt.scatter(mayrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('May')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,6)
    plt.scatter(junerain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('June')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,7)
    plt.scatter(julyrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('July')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,8)
    plt.scatter(augrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('August')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,9)
    plt.scatter(septrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)') 
    plt.title('September')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,10)
    plt.scatter(octrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('October')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,11)
    plt.scatter(novrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('November')
    plt.locator_params(axis="x", nbins=5)
    plt.subplot(4,3,12)
    plt.scatter(decrain,elevation)
    plt.xlim(0,0.02)
    plt.ylabel('Elevation (m)')
    plt.xlabel('Snow (m w.e. $day^{-1}$)')
    plt.title('December')
    plt.locator_params(axis="x", nbins=5)
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'Monthly_SnowvZ.png'))

else:
    pass

#------------------------------------------------------------------------------

if Mean_MB_Timeseries == True:
    print('plotting mean annual runoff')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename2 = 'Daily_Accumulation' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_accumulation_array = np.load(filepath2) #timeseries array
        
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_acc = np.mean(yearly_accumulation_array, axis = 0)

        return [mean_daily_netmelt,mean_daily_acc,yearly_netmelt_array]
    

    allyrs_abl = []
    allyrs_acc = []
    for year in years[1:]:
        allyrs_abl.append(PlotTotalMeltbyYear(year)[0][0:365])
        allyrs_acc.append(PlotTotalMeltbyYear(year)[1][0:365])

    abl_sum = allyrs_abl[0]
    acc_sum = allyrs_acc[0]
    for i in range(1,len(years[1:])):
        abl_sum += allyrs_abl[i]
        acc_sum += allyrs_acc[i]

    
    abl_avg = (abl_sum/len(years[1:]))*(-1)
    acc_avg = acc_sum/len(years[1:])
    
    avg_mb = np.add(abl_avg,acc_avg)
    cumulative_mb = np.cumsum(avg_mb)
    
    plt.figure(figsize=(8,8))
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    plt.title('Mass Balance (2007-2018) \n Debris Case',fontsize=15)
    ax0.plot(abl_avg,'r',linewidth=3,label='Ablation')
    ax0.plot(acc_avg,'b',zorder=5,linewidth=3,label='Accumulation')
    ax1.plot(cumulative_mb,'black',linewidth=4,label='Average Mass Balance: \n' + str(np.round(cumulative_mb[-1],2)) + ' m w.e. year$^{-1}$')
    #plt.legend(['Ablation','Accumulation','Average Mass Balance: ' + str(np.round(cumulative_mb[-1],2)) + ' m w.e. year${^-1}$'],fontsize=15)
    ax1.set_ylabel('Cumulative Mass Balance (m w.e. year$^{-1}$)', color='k',fontsize=15)
    ax0.set_ylabel('Daily Mass Change (m w.e. day$^{-1}$)', color='k',fontsize=15)
    ax0.set_xlabel('Day Of Year',fontsize=15)
    #plt.legend(['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','Mean'])
    #plt.set_xticklabels('jan','feb')
    #month_starts = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    #month_names = ['Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'] 
    plt.margins(x=0)
    ax0.margins(y=0)
    ax1.margins(y=0)
    ax1.set_ylim(-0.9,0.9)
    plt.setp(ax0.get_xticklabels(), fontsize=15)
    plt.setp(ax0.get_yticklabels(), fontsize=15)
    plt.setp(ax1.get_yticklabels(), fontsize=15)
    #ax0.set_ylim(0,0.03) works for runoff in mw.e./day
    ax0.set_ylim(-0.03,0.03)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax0.grid(which='both')
    #ax1.set_yticks(0,1.9)
    ax0.set_zorder(ax1.get_zorder()+1)
    ax0.set_frame_on(False)
    ax0.legend(loc='lower left', frameon=False,fontsize=15)
    ax1.legend(loc='upper right', frameon=False,fontsize=15)
    plt.savefig(os.path.join(Path2files,'Mean_MB_timeseries.pdf'),bbox_inches = 'tight')
    
else:
    pass

#-------------------------------------------------------------------------------------    
    
if Compare_BC_vs_NonBC_Rain == True:

    def PlotDailyMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_noBCrain_array = np.load(filepath)
        
        Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2_BC_files,filename)
        yearly_BCrain_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_noBCrain = np.mean(yearly_noBCrain_array, axis = 0)
        cumulative_noBCrain = np.cumsum(mean_daily_noBCrain)
        
        #mean_daily_BCrain = np.mean(yearly_BCrain_array, axis = 0)
        cumulative_BCrain = np.cumsum(yearly_BCrain_array)
        
        return [mean_daily_noBCrain,cumulative_noBCrain,yearly_BCrain_array,cumulative_BCrain]
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()
    
    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    
    plt.suptitle('Non Bias Corrected vs Bias Corrected Rain \n' + str(debristitle) + '\n' + str(r2stitle), y=1.05)
    
    year = 2007
    ax0.set_title('Rain: ' + str(year))
    ax0.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax0.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax1.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax1.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain') #\n =' + str(round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    #ax1.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain') #\n =' + str(round(PlotDailyMeltbyYear(year)[3][-1],2)) + ' m w.e.')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax1.legend(loc='upper left', frameon=False)
    ax0.margins(x=0)
    
    year = 2008
    ax2.set_title('Rain: ' + str(year))
    ax2.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax3.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax2.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax3.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax3.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax3.legend(loc='upper left', frameon=False)
    ax2.margins(x=0)
    
    year = 2009
    ax4.set_title('Rain: ' + str(year))
    ax4.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax5.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax4.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax5.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax5.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax5.legend(loc='upper left', frameon=False)
    ax4.margins(x=0)
    
    year = 2010
    ax6.set_title('Rain: ' + str(year))
    ax6.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax7.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax6.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax7.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax7.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax7.legend(loc='upper left', frameon=False)
    ax6.margins(x=0)
    
    year = 2011
    ax8.set_title('Rain: ' + str(year))
    ax8.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax9.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax8.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax9.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax9.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax9.legend(loc='upper left', frameon=False)
    ax8.margins(x=0)
    
    year = 2012
    ax10.set_title('Rain: ' + str(year))
    ax10.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax11.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax10.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax11.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax11.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax11.legend(loc='upper left', frameon=False)
    ax10.margins(x=0)
    
    year = 2013
    ax12.set_title('Rain: ' + str(year))
    ax12.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax13.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax12.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax13.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax13.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax13.legend(loc='upper left', frameon=False)
    ax12.margins(x=0)
    
    year = 2014
    ax14.set_title('Rain: ' + str(year))
    ax14.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax15.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax14.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax15.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax15.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax15.legend(loc='upper left', frameon=False)
    ax14.margins(x=0)
    
    year = 2015
    ax16.set_title('Rain: ' + str(year))
    ax16.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax17.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax16.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax17.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax17.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax17.legend(loc='upper left', frameon=False)
    ax16.margins(x=0)
    
    year = 2016
    ax18.set_title('Rain: ' + str(year))
    ax18.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax19.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax18.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax19.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax19.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax19.legend(loc='upper left', frameon=False)
    ax18.margins(x=0)
    
    year = 2017
    ax20.set_title('Rain: ' + str(year))
    ax20.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax21.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax20.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax21.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax21.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax21.legend(loc='upper left', frameon=False)
    ax20.margins(x=0)
    
    year = 2018
    ax22.set_title('Rain: ' + str(year))
    ax22.plot(PlotDailyMeltbyYear(year)[2],'r')
    ax23.plot(PlotDailyMeltbyYear(year)[3],'r',label = 'BC Rain')
    ax22.plot(PlotDailyMeltbyYear(year)[0],'b')
    ax23.plot(PlotDailyMeltbyYear(year)[1],'b',label = 'Non BC Rain')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Rain (m w.e. $day^{-1}$)', color='k')
    ax23.set_ylabel('Cumulative Rain (m w.e.)', color='k')
    ax23.legend(loc='upper left', frameon=False)
    ax22.margins(x=0)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'BCvsNonBC_Rain.pdf'))
else:
    pass
        
#------------------------------------------------------------------------------------- 
    

if Master_BCvNoBC_Rain == True:
    print('plotting mean annual runoff')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs

        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2_BC_files,filename)
        mean_daily_BCrain = np.load(filepath)
        
        mean_daily_NoBCrain = np.mean(yearly_rain_array, axis = 0)

        total_NoBCrain = np.cumsum(mean_daily_NoBCrain)
        total_NoBCrain[0] = 0
        
        total_BCrain = np.cumsum(mean_daily_BCrain)
        total_BCrain[0] = 0

        return [mean_daily_NoBCrain,total_NoBCrain,mean_daily_BCrain,total_BCrain]
    

    allyrs_NoBCrain = []
    allyrs_BCrain = []
    for year in years[1:]:
        allyrs_NoBCrain.append(PlotTotalMeltbyYear(year)[0][0:365])
        allyrs_BCrain.append(PlotTotalMeltbyYear(year)[2][0:365])
        

    NoBCrain_sum = allyrs_NoBCrain[0]
    BCrain_sum = allyrs_BCrain[0]
    for i in range(1,len(years[1:])):
        NoBCrain_sum += allyrs_NoBCrain[i]
        BCrain_sum += allyrs_BCrain[i]
    
    KWArea = 1095554179
    
    NoBCrain_avg = NoBCrain_sum/len(years[1:])
    BCrain_avg = BCrain_sum/len(years[1:])

    NoBCrain_discharge = NoBCrain_avg*(KWArea/86400) #units = m3/s
    BCrain_discharge = BCrain_avg*(KWArea/86400)
    
    cumulative_NoBCrain = np.cumsum(NoBCrain_avg)
    cumulative_BCrain = np.cumsum(BCrain_avg)
    
    plt.figure(figsize=(8,8))
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    plt.title('Bias-Corrected vs Non Bias-Corrected Rain (2007-2018)',fontsize=15)
    ax1.plot(cumulative_BCrain,'darkgreen')
    ax1.plot(cumulative_NoBCrain,'mediumaquamarine')
    ax0.plot(BCrain_discharge,'darkgreen')
    ax0.plot(NoBCrain_discharge,'mediumaquamarine')
    ax0.legend(['BC Rain: ' + str(np.round(cumulative_BCrain[-1],2)) + ' m w.e. year$^{-1}$','Non BC Rain: ' + str(np.round(cumulative_NoBCrain[-1],2)) + ' m w.e. year$^{-1}$'],fontsize=15)
    ax1.set_ylabel('Cumulative Rainfall (m w.e. year$^{-1}$)', color='k',fontsize=15)
    ax0.set_ylabel('Daily Rainfall (m${^3}$ s$^{-1}$)', color='k',fontsize=15)
    ax0.set_xlabel('Day Of Year',fontsize=15)
    #plt.legend(['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','Mean'])
    #plt.set_xticklabels('jan','feb')
    #month_starts = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    #month_names = ['Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'] 
    plt.margins(x=0)
    ax0.margins(y=0)
    ax1.margins(y=0)
    ax1.set_ylim(0,0.25)
    plt.setp(ax0.get_xticklabels(), fontsize=15)
    plt.setp(ax0.get_yticklabels(), fontsize=15)
    plt.setp(ax1.get_yticklabels(), fontsize=15)
    #ax0.set_ylim(0,0.03) works for runoff in mw.e./day
    ax0.set_ylim(0,40)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax0.grid(which='both')
    #ax1.set_yticks(0,1.9)
    ax0.set_zorder(ax1.get_zorder()+1)
    ax0.set_frame_on(False)
    plt.savefig(os.path.join(Path2files,'Master_BCvNoBC_Rain.pdf'),bbox_inches = 'tight')
    
else:
    pass


if Mean_Runoff_Timeseries_withBCRain == True:
    print('plotting mean annual runoff')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        #filename3 = 'Daily_Rain' + str(year) + '.npy'
        #filepath3 = os.path.join(Path2files,filename3)
        #yearly_rain_array = np.load(filepath3) #timeseries array
        
        Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2_BC_files,filename)
        mean_daily_rain = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        #mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        #mean_daily_snowmelt[0] = 0
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0
        
        total_snowmelt = np.cumsum(mean_daily_snowmelt)
        total_snowmelt[0] = 0
        total_icemelt = np.cumsum(mean_daily_icemelt)
        total_icemelt[0] = 0
        total_rain = np.cumsum(mean_daily_rain)
        total_rain[0] = 0

        return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain]
    
    allyrs_ice = []
    allyrs_snow = []
    allyrs_rain = []
    for year in years[1:]:
        allyrs_ice.append(PlotTotalMeltbyYear(year)[0][0:365])
        allyrs_snow.append(PlotTotalMeltbyYear(year)[1][0:365])
        allyrs_rain.append(PlotTotalMeltbyYear(year)[3][0:365])
        
    ice_sum = allyrs_ice[0]
    snow_sum = allyrs_snow[0]
    rain_sum = allyrs_rain[0]
    for i in range(1,len(years[1:])):
        ice_sum += allyrs_ice[i]
        snow_sum += allyrs_snow[i]
        rain_sum += allyrs_rain[i]
    
    KWArea = 1095554179
    
    ice_avg = ice_sum/len(years[1:])
    ice_avg[0] = 0
    snow_avg = snow_sum/len(years[1:])
    snow_avg[0]= 0
    rain_avg = rain_sum/len(years[1:])
    
    ice_discharge = ice_avg*(KWArea/86400) #units = m3/s
    snow_discharge = snow_avg*(KWArea/86400) #units = m3/s
    rain_discharge = rain_avg*(KWArea/86400) #units = m3/s
    
    
    total_runoff = (np.array(ice_avg)) + (np.array(snow_avg)) + (np.array(rain_avg))
    cumulative_runoff = np.nancumsum(total_runoff)
    cumulative_snow = np.cumsum(snow_avg)
    cumulative_ice = np.cumsum(ice_avg)
    cumulative_rain = np.cumsum(rain_avg)
    
    total_discharge = total_runoff*(KWArea/86400) #units = m3/s
    Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
    
    plt.figure(figsize=(8,8))
    ax0 = plt.subplot()
    ax1 = ax0.twinx()
    plt.title('Runoff (2007-2018) \n Debris Case \n Bias Corrected Rain',fontsize=15)
    ax1.plot(cumulative_runoff,'black',zorder=5)
    ax1.plot(cumulative_snow,'b',zorder=3)
    ax1.plot(cumulative_ice,'deeppink',zorder=2)
    ax1.plot(cumulative_rain,'darkgreen',zorder=1)
    #ax0.plot(snow_avg,'b')
    #ax0.fill_between(dayofyear,snow_avg+0.0045,snow_avg-0.0045)
    #ax0.plot(ice_avg,'deeppink',zorder=5)
    #ax0.plot(rain_avg,'mediumaquamarine',zorder=6)
    #ax0.plot(total_runoff,'black',linewidth=2,zorder=7)
    ax0.plot(snow_discharge,'b')
    ax0.plot(ice_discharge,'deeppink',zorder=5)
    ax0.plot(rain_discharge,'darkgreen',zorder=6)
    ax0.plot(total_discharge,'black',linewidth=2,zorder=7)
    plt.legend(['Total Runoff: ' + str(np.round(cumulative_runoff[-1],2)) + ' m w.e.','Snow Melt: ' + str(np.round(cumulative_snow[-1],2)) + ' m w.e.','Ice Melt: ' + str(np.round(cumulative_ice[-1],2)) + ' m w.e.','BC Rain: ' + str(np.round(cumulative_rain[-1],2)) + ' m w.e.'],fontsize=15)
    ax1.set_ylabel('Cumulative Runoff (m w.e. year$^{-1}$)', color='k',fontsize=15)
    ax0.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k',fontsize=15)
    ax0.set_xlabel('Day Of Year',fontsize=15)
    #plt.legend(['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','Mean'])
    #plt.set_xticklabels('jan','feb')
    #month_starts = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    #month_names = ['Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'] 
    plt.margins(x=0)
    ax0.margins(y=0)
    ax1.margins(y=0)
    ax1.set_ylim(0,2.10)
    plt.setp(ax0.get_xticklabels(), fontsize=15)
    plt.setp(ax0.get_yticklabels(), fontsize=15)
    plt.setp(ax1.get_yticklabels(), fontsize=15)
    #ax0.set_ylim(0,0.03) works for runoff in mw.e./day
    ax0.set_ylim(0,360)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax0.grid(which='both')
    #ax1.set_yticks(0,1.9)
    ax0.set_zorder(ax1.get_zorder()+1)
    ax0.set_frame_on(False)
    plt.savefig(os.path.join(Path2_BC_files,'Mean_Runoff_BCRain.pdf'),bbox_inches = 'tight')
    
else:
    pass

if Mean_Runoff_PieChart_withBCRain == True:
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        #filename3 = 'Daily_Rain' + str(year) + '.npy'
        #filepath3 = os.path.join(Path2files,filename3)
        #yearly_rain_array = np.load(filepath3) #timeseries array
        
        Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2_BC_files,filename)
        mean_daily_rain = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        #mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        rain_list = mean_daily_rain.tolist()
        
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        snowmelt_list = mean_daily_snowmelt.tolist()
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        icemelt_list = mean_daily_icemelt.tolist()
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0

        return [icemelt_list,snowmelt_list,total_melt,rain_list]
    
        
    allyears_ice = []
    allyears_snow = []
    allyears_rain = []
    for year in years[1:]:
        allyears_ice = allyears_ice + PlotTotalMeltbyYear(year)[0]
        allyears_snow = allyears_snow + PlotTotalMeltbyYear(year)[1]
        allyears_rain = allyears_rain + PlotTotalMeltbyYear(year)[3]

    daily_icesnowmelt = np.add(allyears_ice,allyears_snow)
    daily_runoff = np.add(daily_icesnowmelt,allyears_rain)
    cumulative_runoff = np.cumsum(daily_runoff)
    cumul_icesnowmelt = np.cumsum(daily_icesnowmelt)
    totalmb = (cumulative_runoff[-1])
    
    cumulative_icemelt = np.cumsum(allyears_ice)
    total_icemelt = cumulative_icemelt[-1]
    cumulative_snowmelt = np.cumsum(allyears_snow)
    total_snowmelt = cumulative_snowmelt[-1]
    cumulative_rain = np.cumsum(allyears_rain)
    total_rain = cumulative_rain[-1]
    
    ice_percent = round((total_icemelt/totalmb)*100,2)
    snow_percent = round((total_snowmelt/totalmb)*100,2)
    rain_percent = round((total_rain/totalmb)*100,2)
    
    Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
    y = np.array([snow_percent,ice_percent,rain_percent])
    mylabels = ['Snow Melt \n' + str(snow_percent) + '%','Ice Melt \n' + str(ice_percent) + '%','BC Rain \n' + str(rain_percent) + '%']
    mycolors = ['b','deeppink','darkgreen']

    plt.figure(figsize=(5,5))
    plt.pie(y, labels = mylabels, colors = mycolors,textprops={'fontsize': 15})
    plt.title('Kaskawulsh Runoff (2007-2018) \n' + str(debristitle) + '\n Bias Corrected Rain',fontsize=15)
    plt.savefig(os.path.join(Path2_BC_files,'Runoff_piechart_BCRain.pdf'),bbox_inches='tight')
    
else:
    pass

if Plot_TotalRunoff_Yearly_withBCRain == True:
    # we assume that: refreezing + net melt = total melt
    print('ice melt = net melt - snow melt')
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        #filename = 'Daily_IceMelt' + str(year) + '.npy'
        #filepath = os.path.join(Path2files,filename)
        #yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        #filename3 = 'Daily_Rain' + str(year) + '.npy'
        #filepath3 = os.path.join(Path2files,filename3)
        #yearly_rain_array = np.load(filepath3) #timeseries array
        
        Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2_BC_files,filename)
        mean_daily_rain = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        #mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0
        
        total_snowmelt = np.cumsum(mean_daily_snowmelt)
        total_snowmelt[0] = 0
        total_icemelt = np.cumsum(mean_daily_icemelt)
        total_icemelt[0] = 0
        total_rain = np.cumsum(mean_daily_rain)
        total_rain[0] = 0

        return [mean_daily_icemelt,mean_daily_snowmelt,total_melt,mean_daily_rain,total_snowmelt,total_icemelt,total_rain,daily_runoff]
    
    KWArea = 1095554179
    Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
    
    plt.figure(figsize=(14,14))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([0,600])
    ax1.set_ylim([0,2.7])
    
    plt.suptitle('Total Runoff by Source \n' + str(debristitle) + '\n Bias Corrected Rain', y=1.05)
    
    year = 2007
    PlotTotalMeltbyYear(year)[1][0] = 0
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax0.set_title('Net Runoff: ' + str(year))
    ax0.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax0.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax0.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax0.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    #ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Cumulative \n Runoff \n =' + str(cumulative_runoff) + ' m w.e.')
    ax1.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax1.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax1.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax1.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax0.margins(x=0)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax0.grid(which='both')
    ax0.set_zorder(ax1.get_zorder()+1)
    ax0.set_frame_on(False)
    ax0.legend(loc='upper left', frameon=False)
    #ax1.legend(loc='upper right', frameon=False)
    
    year = 2008
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax2.set_title('Net Runoff: ' + str(year))
    ax2.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax2.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax2.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax2.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax3.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax3.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax3.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax3.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax3.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax2.margins(x=0)
    ax2.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax3.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax2.grid(which='both')
    ax2.set_zorder(ax3.get_zorder()+1)
    ax2.set_frame_on(False)
    ax2.legend(loc='upper left', frameon=False)
    #ax3.legend(loc='upper right', frameon=False)
    
    year = 2009
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax4.set_title('Net Runoff: ' + str(year))
    ax4.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax4.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax4.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax4.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax5.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax5.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax5.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax5.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax5.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax4.margins(x=0)
    ax4.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax5.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax4.grid(which='both')
    ax4.set_zorder(ax5.get_zorder()+1)
    ax4.set_frame_on(False)
    ax4.legend(loc='upper left', frameon=False)
    #ax5.legend(loc='center left', frameon=False)
    
    year = 2010
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax6.set_title('Net Runoff: ' + str(year))
    ax6.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax6.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax6.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax6.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax7.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax7.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax7.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax7.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax7.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax6.margins(x=0)
    ax6.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax7.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax6.grid(which='both')
    ax6.set_zorder(ax7.get_zorder()+1)
    ax6.set_frame_on(False)
    ax6.legend(loc='upper left', frameon=False)
    #ax7.legend(loc='center left', frameon=False)
    
    year = 2011
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax8.set_title('Net Runoff: ' + str(year))
    ax8.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax8.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax8.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax8.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax9.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax9.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax9.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax9.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax9.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax8.margins(x=0)
    ax8.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax9.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax8.grid(which='both')
    ax8.set_zorder(ax9.get_zorder()+1)
    ax8.set_frame_on(False)
    ax8.legend(loc='upper left', frameon=False)
    #ax9.legend(loc='upper right', frameon=False)
    
    year = 2012
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax10.set_title('Net Runoff: ' + str(year))
    ax10.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax10.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax10.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax10.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax11.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax11.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax11.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax11.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax11.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax10.margins(x=0)
    ax10.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax11.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax10.grid(which='both')
    ax10.set_zorder(ax11.get_zorder()+1)
    ax10.set_frame_on(False)
    ax10.legend(loc='upper left', frameon=False)
    #ax11.legend(loc='upper right', frameon=False)
    
    year = 2013
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax12.set_title('Net Runoff: ' + str(year))
    ax12.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax12.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax12.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax12.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax13.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax13.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax13.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax13.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax13.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax12.margins(x=0)
    ax12.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax13.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax12.grid(which='both')
    ax12.set_zorder(ax13.get_zorder()+1)
    ax12.set_frame_on(False)
    ax12.legend(loc='upper left', frameon=False)
    #ax13.legend(loc='center left', frameon=False)
    
    year = 2014
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax14.set_title('Net Runoff: ' + str(year))
    ax14.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax14.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax14.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax14.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax15.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax15.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax15.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax15.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax15.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax14.margins(x=0)
    ax14.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax15.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax14.grid(which='both')
    ax14.set_zorder(ax15.get_zorder()+1)
    ax14.set_frame_on(False)
    ax14.legend(loc='upper left', frameon=False)
    #ax15.legend(loc='upper right', frameon=False)
    
    year = 2015
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax16.set_title('Net Runoff: ' + str(year))
    ax16.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax16.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax16.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax16.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax17.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax17.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax17.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax17.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax17.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax16.margins(x=0)
    ax16.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax17.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax16.grid(which='both')
    ax16.set_zorder(ax17.get_zorder()+1)
    ax16.set_frame_on(False)
    ax16.legend(loc='upper left', frameon=False)
    #ax17.legend(loc='upper right', frameon=False)
    
    year = 2016
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax18.set_title('Net Runoff: ' + str(year))
    ax18.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax18.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax18.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax18.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax19.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax19.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax19.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax19.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax19.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax18.margins(x=0)
    ax18.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax19.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax18.grid(which='both')
    ax18.set_zorder(ax19.get_zorder()+1)
    ax18.set_frame_on(False)
    ax18.legend(loc='upper left', frameon=False)
    #ax19.legend(loc='center left', frameon=False)
    
    year = 2017
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax20.set_title('Net Runoff: ' + str(year))
    ax20.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax20.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax20.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax20.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax21.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax21.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax21.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax21.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax21.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax20.margins(x=0)
    ax20.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax21.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax20.grid(which='both')
    ax20.set_zorder(ax21.get_zorder()+1)
    ax20.set_frame_on(False)
    ax20.legend(loc='upper left', frameon=False)
    #ax21.legend(loc='upper right', frameon=False)
    
    year = 2018
    cumulative_runoff = round((PlotTotalMeltbyYear(year)[2])[-1],2)
    ax22.set_title('Net Runoff: ' + str(year))
    ax22.plot(PlotTotalMeltbyYear(year)[7]*(KWArea/86400),'k', label = 'Total \n Runoff')
    ax22.plot(PlotTotalMeltbyYear(year)[1]*(KWArea/86400), 'b', label = 'Snow Melt')
    ax22.plot(PlotTotalMeltbyYear(year)[0]*(KWArea/86400),'deeppink', label ='Ice Melt')
    ax22.plot(PlotTotalMeltbyYear(year)[3]*(KWArea/86400),'darkgreen', label = 'BC Rain')
    ax23.plot(PlotTotalMeltbyYear(year)[2],'k', label = 'Total \n Runoff')
    ax23.plot(PlotTotalMeltbyYear(year)[4],'b')
    ax23.plot(PlotTotalMeltbyYear(year)[5],'deeppink')
    ax23.plot(PlotTotalMeltbyYear(year)[6],'darkgreen')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Runoff (m${^3}$ s$^{-1}$)', color='k')
    ax23.set_ylabel('Cumulative Runoff (m w.e.)', color='k')
    ax22.margins(x=0)
    ax22.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax23.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax22.grid(which='both')
    ax22.set_zorder(ax23.get_zorder()+1)
    ax22.set_frame_on(False)
    ax22.legend(loc='upper left', frameon=False)
    #ax23.legend(loc='upper right', frameon=False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2_BC_files,'Yearly_Runoff_BCRain.pdf'),bbox_inches='tight')

else:
    pass

if Yearly_Runoff_PieChart_withBCRain == True:
    def PlotTotalMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        #filename3 = 'Daily_Rain' + str(year) + '.npy'
        #filepath3 = os.path.join(Path2files,filename3)
        #yearly_rain_array = np.load(filepath3) #timeseries array
        
        Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
        filename = 'Daily_Rain' + str(year) + '.npy'
        filepath = os.path.join(Path2_BC_files,filename)
        mean_daily_rain = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        #mean_daily_icemelt = np.mean(yearly_icemelt_array, axis = 0)
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        #mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        rain_list = mean_daily_rain.tolist()
        cumrain = np.cumsum(rain_list)
        
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        snowmelt_list = mean_daily_snowmelt.tolist()
        cumsnow = np.cumsum(snowmelt_list)
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        icemelt_list = mean_daily_icemelt.tolist()
        cumice = np.cumsum(icemelt_list)
        
        daily_total_melt = np.add(mean_daily_snowmelt,mean_daily_icemelt)
        daily_runoff = np.add(daily_total_melt,mean_daily_rain)
        total_melt = np.cumsum(daily_runoff)
        total_melt[0] = 0
        
        yr_total = total_melt[-1]
        total_ice = cumice[-1]
        total_snow = cumsnow[-1]
        total_rain = cumrain[-1]
        
        ice_percent = round((total_ice/yr_total)*100,1)
        snow_percent = round((total_snow/yr_total)*100,1)
        rain_percent = round((total_rain/yr_total)*100,1)
        
        y = np.array([snow_percent,ice_percent,rain_percent])
        mylabels = ['Snow Melt \n' + str(snow_percent) + '%','Ice Melt \n' + str(ice_percent) + '%','BC Rain \n' + str(rain_percent) + '%']

        return [y,mylabels]
    
    mycolors = ['b','deeppink','darkgreen']
    Path2_BC_files = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\BiasCorrected_Rain'
    
    plt.figure(figsize=(10,12))
    plt.suptitle('Kaskawulsh Runoff (2007-2018) \n' + str(debristitle) + '\n Bias Corrected Rain',fontsize=15)
    plt.subplot(4,3,1)
    plt.title('2007',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2007)[0], labels = PlotTotalMeltbyYear(2007)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,2)
    plt.title('2008',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2008)[0], labels = PlotTotalMeltbyYear(2008)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,3)
    plt.title('2009',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2009)[0], labels = PlotTotalMeltbyYear(2009)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,4)
    plt.title('2010',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2010)[0], labels = PlotTotalMeltbyYear(2010)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,5)
    plt.title('2011',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2011)[0], labels = PlotTotalMeltbyYear(2011)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,6)
    plt.title('2012',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2012)[0], labels = PlotTotalMeltbyYear(2012)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,7)
    plt.title('2013',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2013)[0], labels = PlotTotalMeltbyYear(2013)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,8)
    plt.title('2014',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2014)[0], labels = PlotTotalMeltbyYear(2014)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,9)
    plt.title('2015',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2015)[0], labels = PlotTotalMeltbyYear(2015)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,10)
    plt.title('2016',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2016)[0], labels = PlotTotalMeltbyYear(2016)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,11)
    plt.title('2017',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2017)[0], labels = PlotTotalMeltbyYear(2017)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.subplot(4,3,12)
    plt.title('2018',fontsize=12)
    plt.pie(PlotTotalMeltbyYear(2018)[0], labels = PlotTotalMeltbyYear(2018)[1], colors = mycolors,textprops={'fontsize': 12})
    plt.savefig(os.path.join(Path2_BC_files,'Yearly_Runoff_piechart_withBCRain.pdf'),bbox_inches='tight')
    
    
else:
    pass

if Plot_NARR_distributed_acc_and_accVSz == True:
    print('NARR downscaled Precip')
    
    out_run_means = 'allrunPrecipitation_' + 'DownscaledNoBC' + '.npy'
    out_allruns = os.path.join(Path2files,out_run_means)
    overall_Precipitation = np.load(out_allruns)
    
    out_run_rains = 'allrunRain_' + 'DownscaledNoBC' + '.npy'
    out_allruns = os.path.join(Path2files,out_run_rains)
    overall_Rain = np.load(out_allruns)
    
    overall_Accumulation = overall_Precipitation - overall_Rain
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD:    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(overall_Accumulation), cmap = 'BuPu', levels = np.linspace(0,4, 21))
    #plt.contourf(np.flipud(Topo), cmap = 'PuRd', levels = np.linspace(0,4, 20))
    legend = plt.colorbar()
    legend.ax.tick_params(labelsize=15) 
    legend.ax.set_ylabel('Accumulation (m w.e. $a^{-1}$)', rotation=270,labelpad=20,fontsize=15)
    #legend.ax.set_ylabel('Bed Elevation (m a.s.l)',labelpad=15, rotation=270,fontsize=15)
    plt.title('Bias Corrected Accumulation (2007-2008)',fontsize=15)
    plt.savefig(os.path.join(Path2files,'Mean_Accumulation_UnCorrected.pdf'))
    
    accumulation = []
    elevation = []
    for x in range(len(Zgrid)):
        for y in range(len(Zgrid[0])):
            z = Zgrid[x][y]
            elevation.append(z)
            acc = overall_Accumulation[x][y]
            accumulation.append(acc)
    
    
    plt.figure(figsize=(12,8))
    plt.scatter(accumulation,elevation)
    #legend = plt.colorbar()
    #legend.ax.set_ylabel('Mass Balance (m w.e. $a^{-1}$)', rotation=270)
    plt.ylabel('Elevation (m)',fontsize=14)
    plt.xlabel('Accumulation (m w.e. $a^{-1}$)',fontsize=14)
    plt.title('Uncorrected Accumulation vs Elevation' + '\n' + str(r2stitle),fontsize=14)
    plt.savefig(os.path.join(Path2files,'AccumulationNoBC_vs_Z.pdf'))
    
    np.save(os.path.join(Path2files,'UncorrectedAccumulation.npy'),overall_Accumulation)
    
    print('Mean Melt plotted from the file: ' + str(out_run_means))
else:
    pass

if Plot_IceMelt_byYear == True: 
    def PlotDailyMeltbyYear(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0)
        cumulative_melt = np.cumsum(mean_daily_melt)
        
        return [mean_daily_melt,cumulative_melt]
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    
    year = 2007
    ax0.set_title('Net Melt: ' + str(year))
    ax0.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax1.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax1.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax1.legend(loc='upper left', frameon=False)
    
    year = 2008
    ax2.set_title('Net Melt: ' + str(year))
    ax2.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax3.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax3.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax3.legend(loc='upper left', frameon=False)
    
    
    year = 2009
    ax4.set_title('Net Melt: ' + str(year))
    ax4.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax5.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax5.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax5.legend(loc='upper left', frameon=False)
    
    
    year = 2010
    ax6.set_title('Net Melt: ' + str(year))
    ax6.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax7.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax7.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax7.legend(loc='upper left', frameon=False)
    
    
    year = 2011
    ax8.set_title('Net Melt: ' + str(year))
    ax8.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax9.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax9.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax9.legend(loc='upper left', frameon=False)
    
    
    year = 2012
    ax10.set_title('Net Melt: ' + str(year))
    ax10.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax11.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax11.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax11.legend(loc='upper left', frameon=False)
    
    
    year = 2013
    ax12.set_title('Net Melt: ' + str(year))
    ax12.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax13.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax13.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax13.legend(loc='upper left', frameon=False)
    
    
    year = 2014
    ax14.set_title('Net Melt: ' + str(year))
    ax14.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax15.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax15.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax15.legend(loc='upper left', frameon=False)
    
    
    year = 2015
    ax16.set_title('Net Melt: ' + str(year))
    ax16.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax17.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax17.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax17.legend(loc='upper left', frameon=False)
    
    
    year = 2016
    ax18.set_title('Net Melt: ' + str(year))
    ax18.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax19.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax19.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax19.legend(loc='upper left', frameon=False)
    
    
    year = 2017
    ax20.set_title('Net Melt: ' + str(year))
    ax20.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax21.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax21.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax21.legend(loc='upper left', frameon=False)
    
    
    year = 2018
    ax22.set_title('Net Melt: ' + str(year))
    ax22.plot(PlotDailyMeltbyYear(year)[0],'c')
    ax23.plot(PlotDailyMeltbyYear(year)[1],'k',label = 'Total \n Melt \n =' + str(np.round(PlotDailyMeltbyYear(year)[1][-1],2)) + ' m w.e.')
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Net Melt (m w.e. $day^{-1}$)', color='c')
    ax23.set_ylabel('Cumulative Melt (m w.e.)', color='k')
    ax23.legend(loc='upper left', frameon=False)
    
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'YearlyNetMelt.pdf'))
else:
    pass

if Stacked_barchart_monthly_runoff == True:    
    print('average monthly runoff over entire simulation period')
    def PlotTotalRunoffbyMonth(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_netmelt_array = np.load(filepath)
        
        filename2 = 'Daily_SnowMelt' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_snowmelt_array = np.load(filepath2) #timeseries array
        
        filename3 = 'Daily_Rain' + str(year) + '.npy'
        filepath3 = os.path.join(Path2files,filename3)
        yearly_rain_array = np.load(filepath3) #timeseries array
        
        mean_daily_netmelt = np.mean(yearly_netmelt_array, axis = 0)
        mean_daily_rain = np.mean(yearly_rain_array, axis = 0)
        mean_daily_snowmelt = np.mean(yearly_snowmelt_array, axis=0)
        mean_daily_snowmelt[0] = 0
        #negsnow = np.where(mean_daily_snowmelt < 0)
        for i in range(0,len(mean_daily_snowmelt)):
            if mean_daily_snowmelt[i] < 0:
                mean_daily_snowmelt[i] = mean_daily_snowmelt[i-1]
            else:
                pass
        
        mean_daily_icemelt = mean_daily_netmelt - mean_daily_snowmelt
        neglocs = np.where(mean_daily_icemelt < 0) #NEW
        mean_daily_icemelt[neglocs] = 0 #NEW
        mean_daily_icemelt[0] = 0
        
        #mean_daily_snowmelt = mean_daily_netmelt - mean_daily_icemelt
        
        # calculate snowmelt 
        jan = np.cumsum(mean_daily_snowmelt[0:31])[-1]
        feb = np.cumsum(mean_daily_snowmelt[31:59])[-1]
        march = np.cumsum(mean_daily_snowmelt[59:90])[-1]
        april = np.cumsum(mean_daily_snowmelt[90:120])[-1]
        may = np.cumsum(mean_daily_snowmelt[120:151])[-1]
        june = np.cumsum(mean_daily_snowmelt[151:181])[-1]
        july = np.cumsum(mean_daily_snowmelt[181:212])[-1]
        aug = np.cumsum(mean_daily_snowmelt[212:243])[-1]
        sept = np.cumsum(mean_daily_snowmelt[243:273])[-1]
        octo = np.cumsum(mean_daily_snowmelt[273:304])[-1]
        nov = np.cumsum(mean_daily_snowmelt[304:334])[-1]
        dec = np.cumsum(mean_daily_snowmelt[334:365])[-1]
        monthly_snowmelt_l = [jan,feb,march,april,may,june,july,aug,sept,octo,nov,dec]
        monthly_snowmelt = np.array(monthly_snowmelt_l)
        
        jan = np.cumsum(mean_daily_icemelt[0:31])[-1]
        feb = np.cumsum(mean_daily_icemelt[31:59])[-1]
        march = np.cumsum(mean_daily_icemelt[59:90])[-1]
        april = np.cumsum(mean_daily_icemelt[90:120])[-1]
        may = np.cumsum(mean_daily_icemelt[120:151])[-1]
        june = np.cumsum(mean_daily_icemelt[151:181])[-1]
        july = np.cumsum(mean_daily_icemelt[181:212])[-1]
        aug = np.cumsum(mean_daily_icemelt[212:243])[-1]
        sept = np.cumsum(mean_daily_icemelt[243:273])[-1]
        octo = np.cumsum(mean_daily_icemelt[273:304])[-1]
        nov = np.cumsum(mean_daily_icemelt[304:334])[-1]
        dec = np.cumsum(mean_daily_icemelt[334:365])[-1]
        monthly_icemelt_l = [jan,feb,march,april,may,june,july,aug,sept,octo,nov,dec]
        monthly_icemelt = np.array(monthly_icemelt_l)
        
        jan = np.cumsum(mean_daily_rain[0:31])[-1]
        feb = np.cumsum(mean_daily_rain[31:59])[-1]
        march = np.cumsum(mean_daily_rain[59:90])[-1]
        april = np.cumsum(mean_daily_rain[90:120])[-1]
        may = np.cumsum(mean_daily_rain[120:151])[-1]
        june = np.cumsum(mean_daily_rain[151:181])[-1]
        july = np.cumsum(mean_daily_rain[181:212])[-1]
        aug = np.cumsum(mean_daily_rain[212:243])[-1]
        sept = np.cumsum(mean_daily_rain[243:273])[-1]
        octo = np.cumsum(mean_daily_rain[273:304])[-1]
        nov = np.cumsum(mean_daily_rain[304:334])[-1]
        dec = np.cumsum(mean_daily_rain[334:365])[-1] 
        monthly_rain_l = [jan,feb,march,april,may,june,july,aug,sept,octo,nov,dec]
        monthly_rain = np.array(monthly_rain_l)
        
        
        return [monthly_snowmelt,monthly_icemelt,monthly_rain]
    
    allyrs_snowmelt = PlotTotalRunoffbyMonth(2007)[0] + PlotTotalRunoffbyMonth(2008)[0] +PlotTotalRunoffbyMonth(2009)[0] + PlotTotalRunoffbyMonth(2010)[0] + PlotTotalRunoffbyMonth(2011)[0] + PlotTotalRunoffbyMonth(2012)[0] + PlotTotalRunoffbyMonth(2013)[0] + PlotTotalRunoffbyMonth(2014)[0] + PlotTotalRunoffbyMonth(2015)[0] + PlotTotalRunoffbyMonth(2016)[0] + PlotTotalRunoffbyMonth(2017)[0] + PlotTotalRunoffbyMonth(2018)[0]
    avg_monthly_snowmelt = (allyrs_snowmelt/12)*(KWArea/2592000)
    
    allyrs_icemelt = PlotTotalRunoffbyMonth(2007)[1] + PlotTotalRunoffbyMonth(2008)[1] +PlotTotalRunoffbyMonth(2009)[1] + PlotTotalRunoffbyMonth(2010)[1] + PlotTotalRunoffbyMonth(2011)[1] + PlotTotalRunoffbyMonth(2012)[1] + PlotTotalRunoffbyMonth(2013)[1] + PlotTotalRunoffbyMonth(2014)[1] + PlotTotalRunoffbyMonth(2015)[1] + PlotTotalRunoffbyMonth(2016)[1] + PlotTotalRunoffbyMonth(2017)[1] + PlotTotalRunoffbyMonth(2018)[1]
    avg_monthly_icemelt = (allyrs_icemelt/12)*(KWArea/2592000)
    
    allyrs_rain = PlotTotalRunoffbyMonth(2007)[2] + PlotTotalRunoffbyMonth(2008)[2] +PlotTotalRunoffbyMonth(2009)[2] + PlotTotalRunoffbyMonth(2010)[2] + PlotTotalRunoffbyMonth(2011)[2] + PlotTotalRunoffbyMonth(2012)[2] + PlotTotalRunoffbyMonth(2013)[2] + PlotTotalRunoffbyMonth(2014)[2] + PlotTotalRunoffbyMonth(2015)[2] + PlotTotalRunoffbyMonth(2016)[2] + PlotTotalRunoffbyMonth(2017)[2] + PlotTotalRunoffbyMonth(2018)[2]
    avg_monthly_rain = allyrs_rain/12*(KWArea/2592000)
    
    labels = ['Jan', 'Feb', 'March', 'April', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
    width = 0.75       # the width of the bars: can also be len(x) sequence
    
    
    fig, ax = plt.subplots(figsize=(8,6))
    ax.bar(labels, avg_monthly_rain, width, label='Rain',color='mediumaquamarine')
    ax.bar(labels, avg_monthly_icemelt, width, bottom=avg_monthly_rain,label='Ice Melt',color='deeppink')
    ax.bar(labels, avg_monthly_snowmelt, width, bottom=(avg_monthly_icemelt+np.array(avg_monthly_rain)),label='Snow Melt',color='b')
    
    ax.set_ylabel('Discharge (m${^3}$ s${^-1}$)',fontsize=15)
    plt.xticks(fontsize=15,rotation=45)
    plt.yticks(fontsize=15)
    plt.ylim(0,360)
    ax.set_title('Monthly Runoff (2007-2018) \n Debris Case',fontsize = 15)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'MonthlyRunoff_DebrisCase.pdf'))
else:
    pass

if Yearly_MB_HydrologicYear == True:
    # hydrological year = oct 1 - sept 31
    print('total ablation = net melt + refreezing')
    def PlotMB_fulltimeseries(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        filename2 = 'Daily_Accumulation' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_accumulation_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0) * (-1)
        dailymelt_list = mean_daily_melt.tolist()
        
        mean_daily_refreezing = np.mean(yearly_accumulation_array, axis=0)
        dailyacc_list = mean_daily_refreezing.tolist()
        
        #daily_total_melt = np.add(mean_daily_refreezing,mean_daily_melt)
        #total_melt = np.cumsum(daily_total_melt)
        #total_melt[0] = 0

        return [dailymelt_list,dailyacc_list]
    
    allyears_melt = []
    allyears_acc = []
    for year in years[1:]:
        allyears_melt = allyears_melt + PlotMB_fulltimeseries(year)[0]
        allyears_acc = allyears_acc + PlotMB_fulltimeseries(year)[1]
        
    daily_mb = np.add(allyears_acc,allyears_melt)
    cumulative_mb = np.cumsum(daily_mb)
    totalmb = round(cumulative_mb[-1],2)
    
    #index for oct1 = 274
    acc_0708 = allyears_acc[274:640]
    acc_0809 = allyears_acc[641:1006]
    acc_0910 = allyears_acc[1007:1372]
    acc_1011 = allyears_acc[1373:1738]
    acc_1112 = allyears_acc[1739:2105]
    acc_1213 = allyears_acc[2106:2471]
    acc_1314 = allyears_acc[2472:2837]
    acc_1415 = allyears_acc[2838:3203]
    acc_1516 = allyears_acc[3204:3570]
    acc_1617 = allyears_acc[3571:3936]
    acc_1718 = allyears_acc[3937:4302]
    acc_18 = allyears_acc[4302:] #only 81 days
    
    abl_0708 = allyears_melt[274:640]
    abl_0809 = allyears_melt[641:1006]
    abl_0910 = allyears_melt[1007:1372]
    abl_1011 = allyears_melt[1373:1738]
    abl_1112 = allyears_melt[1739:2105]
    abl_1213 = allyears_melt[2106:2471]
    abl_1314 = allyears_melt[2472:2837]
    abl_1415 = allyears_melt[2838:3203]
    abl_1516 = allyears_melt[3204:3570]
    abl_1617 = allyears_melt[3571:3936]
    abl_1718 = allyears_melt[3937:4302]
    abl_18 = allyears_melt[4302:] #only 81 days
    
    mb_0708 = np.cumsum(np.add(abl_0708,acc_0708))
    mb_0809 = np.cumsum(np.add(abl_0809,acc_0809))
    mb_0910 = np.cumsum(np.add(abl_0910,acc_0910))
    mb_1011 = np.cumsum(np.add(abl_1011,acc_1011))
    mb_1112 = np.cumsum(np.add(abl_1112,acc_1112))
    mb_1213 = np.cumsum(np.add(abl_1213,acc_1213))
    mb_1314 = np.cumsum(np.add(abl_1314,acc_1314))
    mb_1415 = np.cumsum(np.add(abl_1415,acc_1415))
    mb_1516 = np.cumsum(np.add(abl_1516,acc_1516))
    mb_1617 = np.cumsum(np.add(abl_1617,acc_1617))
    mb_1718 = np.cumsum(np.add(abl_1718,acc_1718))
    mb_18 = np.cumsum(np.add(abl_18,acc_18)) #only 81 days
    
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([-0.06,0.06])
    ax1.set_ylim([-1.8, 1.8])
    
    plt.suptitle('Mass Balance \n' + str(debristitle), y=1.05)
    
    year = '2007-2008'
    ax0.set_title('Mass Balance: ' + str(year))
    ax0.plot(abl_0708,'r', label ='Ablation')
    ax0.plot(acc_0708, 'b', label = 'Accumulation')
    ax1.plot(mb_0708,'k', label = 'Cumulative \n Mass Balance' )
    #ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax0.margins(x=0)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax0.grid(which='both')
    #ax1.legend(loc='upper right', frameon=False)
    ax0.legend(loc='lower left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    year_starts = [0,31,61,92,122,153,183,214,244,275,305,336]
    year_names = ['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'] 
    ax0.set_xticks(year_starts)
    ax0.set_xticklabels(year_names,rotation=90)
    
    year = '2008-2009'
    ax2.set_title('Mass Balance: ' + str(year))
    ax2.plot(abl_0809,'r', label ='Ablation')
    ax2.plot(acc_0809, 'b', label = 'Accumulation')
    ax3.plot(mb_0809,'k', label = 'Cumulative \n Mass Balance' )
    #ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax3.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax2.margins(x=0)
    ax2.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax3.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax2.grid(which='both')
    ax2.legend(loc='lower left', frameon=False)
    ax3.legend(loc='upper right', frameon=False)
    ax2.set_xticks(year_starts)
    ax2.set_xticklabels(year_names,rotation=90)
    
    year = '2009-2010'
    ax4.set_title('Mass Balance: ' + str(year))
    ax4.plot(abl_0910,'r', label ='Ablation')
    ax4.plot(acc_0910, 'b', label = 'Accumulation')
    ax5.plot(mb_0910,'k', label = 'Cumulative \n Mass Balance' )
    #ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax5.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax4.margins(x=0)
    ax4.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax5.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax4.grid(which='both')
    ax4.legend(loc='lower left', frameon=False)
    ax5.legend(loc='upper right', frameon=False)
    ax4.set_xticks(year_starts)
    ax4.set_xticklabels(year_names,rotation=90)
    
    year = '2010-2011'
    ax6.set_title('Mass Balance: ' + str(year))
    ax6.plot(abl_1011,'r', label ='Ablation')
    ax6.plot(acc_1011, 'b', label = 'Accumulation')
    ax7.plot(mb_1011,'k', label = 'Cumulative \n Mass Balance' )
    #ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax7.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax6.margins(x=0)
    ax6.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax7.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax6.grid(which='both')
    ax6.legend(loc='lower left', frameon=False)
    ax7.legend(loc='upper right', frameon=False)
    ax6.set_xticks(year_starts)
    ax6.set_xticklabels(year_names,rotation=90)
    
    year = '2011-2012'
    ax8.set_title('Mass Balance: ' + str(year))
    ax8.plot(abl_1112,'r', label ='Ablation')
    ax8.plot(acc_1112, 'b', label = 'Accumulation')
    ax9.plot(mb_1112,'k', label = 'Cumulative \n Mass Balance' )
    #ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax9.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax8.margins(x=0)
    ax8.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax9.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax8.grid(which='both')
    ax8.legend(loc='lower left', frameon=False)
    ax9.legend(loc='upper right', frameon=False)
    ax8.set_xticks(year_starts)
    ax8.set_xticklabels(year_names,rotation=90)
    
    year = '2012-2013'
    ax10.set_title('Mass Balance: ' + str(year))
    ax10.plot(abl_1213,'r', label ='Ablation')
    ax10.plot(acc_1213, 'b', label = 'Accumulation')
    ax11.plot(mb_1213,'k', label = 'Cumulative \n Mass Balance' )
    #ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax11.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax10.margins(x=0)
    ax10.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax11.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax10.grid(which='both')
    ax10.legend(loc='lower left', frameon=False)
    ax11.legend(loc='upper right', frameon=False)
    ax10.set_xticks(year_starts)
    ax10.set_xticklabels(year_names,rotation=90)
    
    year = '2013-2014'
    ax12.set_title('Mass Balance: ' + str(year))
    ax12.plot(abl_1314,'r', label ='Ablation')
    ax12.plot(acc_1314, 'b', label = 'Accumulation')
    ax13.plot(mb_1314,'k', label = 'Cumulative \n Mass Balance' )
    #ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax13.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax12.margins(x=0)
    ax12.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax13.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax12.grid(which='both')
    ax12.legend(loc='lower left', frameon=False)
    ax13.legend(loc='upper right', frameon=False)
    ax12.set_xticks(year_starts)
    ax12.set_xticklabels(year_names,rotation=90)
    
    year = '2014-2015'
    ax14.set_title('Mass Balance: ' + str(year))
    ax14.plot(abl_1415,'r', label ='Ablation')
    ax14.plot(acc_1415, 'b', label = 'Accumulation')
    ax15.plot(mb_1415,'k', label = 'Cumulative \n Mass Balance' )
    #ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax15.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax14.margins(x=0)
    ax14.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax15.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax14.grid(which='both')
    ax14.legend(loc='lower left', frameon=False)
    ax15.legend(loc='upper right', frameon=False)
    ax14.set_xticks(year_starts)
    ax14.set_xticklabels(year_names,rotation=90)
    
    year = '2015-2016'
    ax16.set_title('Mass Balance: ' + str(year))
    ax16.plot(abl_1516,'r', label ='Ablation')
    ax16.plot(acc_1516, 'b', label = 'Accumulation')
    ax17.plot(mb_1516,'k', label = 'Cumulative \n Mass Balance' )
    #ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax17.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax16.margins(x=0)
    ax16.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax17.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax16.grid(which='both')
    ax16.legend(loc='lower left', frameon=False)
    ax17.legend(loc='upper right', frameon=False)
    ax16.set_xticks(year_starts)
    ax16.set_xticklabels(year_names,rotation=90)
    
    year = '2016-2017'
    ax18.set_title('Mass Balance: ' + str(year))
    ax18.plot(abl_1617,'r', label ='Ablation')
    ax18.plot(acc_1617, 'b', label = 'Accumulation')
    ax19.plot(mb_1617,'k', label = 'Cumulative \n Mass Balance' )
    #ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax19.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax18.margins(x=0)
    ax18.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax19.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax18.grid(which='both')
    ax18.legend(loc='lower left', frameon=False)
    ax19.legend(loc='upper right', frameon=False)
    ax18.set_xticks(year_starts)
    ax18.set_xticklabels(year_names,rotation=90)

    year = '2017-2018'
    ax20.set_title('Mass Balance: ' + str(year))
    ax20.plot(abl_1718,'r', label ='Ablation')
    ax20.plot(acc_1718, 'b', label = 'Accumulation')
    ax21.plot(mb_1718,'k', label = 'Cumulative \n Mass Balance' )
    #ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax21.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax20.margins(x=0)
    ax20.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax21.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax20.grid(which='both')
    ax20.legend(loc='lower left', frameon=False)
    ax21.legend(loc='upper right', frameon=False)
    ax20.set_xticks(year_starts)
    ax20.set_xticklabels(year_names,rotation=90)
    
    year = 2018
    ax22.set_title('Mass Balance: ' + str(year))
    ax22.plot(abl_18,'r', label ='Ablation')
    ax22.plot(acc_18, 'b', label = 'Accumulation')
    ax23.plot(mb_18,'k', label = 'Cumulative \n Mass Balance' )
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax23.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax22.margins(x=0)
    ax22.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax23.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax22.grid(which='both')
    ax22.legend(loc='lower left', frameon=False)
    ax23.legend(loc='upper right', frameon=False)
    #ax0.set_xticks(year_starts)
    #ax0.set_xticklabels(year_names,rotation=90)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'HydrologicYear_MB_DebrisFreeCase.pdf'),bbox_inches='tight')

else:
    pass

if Yearly_MB_SummerMelt == True:
    # year = May1st-April 30th
    print('total ablation = net melt + refreezing')
    def PlotMB_fulltimeseries(year):
        # load data file for this year from the extracted outputs
        filename = 'Daily_Net_Melt' + str(year) + '.npy'
        filepath = os.path.join(Path2files,filename)
        yearly_melt_array = np.load(filepath)
        
        filename2 = 'Daily_Accumulation' + str(year) + '.npy'
        filepath2 = os.path.join(Path2files,filename2)
        yearly_accumulation_array = np.load(filepath2) #timeseries array
        
        
        # calculate mean daily net melt for this year across all 10 simulations
        mean_daily_melt = np.mean(yearly_melt_array, axis = 0) * (-1)
        dailymelt_list = mean_daily_melt.tolist()
        
        mean_daily_refreezing = np.mean(yearly_accumulation_array, axis=0)
        dailyacc_list = mean_daily_refreezing.tolist()
        
        #daily_total_melt = np.add(mean_daily_refreezing,mean_daily_melt)
        #total_melt = np.cumsum(daily_total_melt)
        #total_melt[0] = 0

        return [dailymelt_list,dailyacc_list]
    
    allyears_melt = []
    allyears_acc = []
    for year in years[1:]:
        allyears_melt = allyears_melt + PlotMB_fulltimeseries(year)[0]
        allyears_acc = allyears_acc + PlotMB_fulltimeseries(year)[1]
        
    daily_mb = np.add(allyears_acc,allyears_melt)
    cumulative_mb = np.cumsum(daily_mb)
    totalmb = round(cumulative_mb[-1],2)
    
    #index for may 1: 
    acc_0708 = allyears_acc[121:487]
    acc_0809 = allyears_acc[488:853]
    acc_0910 = allyears_acc[854:1219]
    acc_1011 = allyears_acc[1220:1585]
    acc_1112 = allyears_acc[1586:1952]
    acc_1213 = allyears_acc[1953:2318]
    acc_1314 = allyears_acc[2319:2684]
    acc_1415 = allyears_acc[2685:3050]
    acc_1516 = allyears_acc[3051:3417]
    acc_1617 = allyears_acc[3418:3783]
    acc_1718 = allyears_acc[3784:4149]
    acc_18 = allyears_acc[4150:] #only 81 days
    
    abl_0708 = allyears_melt[121:487]
    abl_0809 = allyears_melt[488:853]
    abl_0910 = allyears_melt[854:1219]
    abl_1011 = allyears_melt[1220:1585]
    abl_1112 = allyears_melt[1586:1952]
    abl_1213 = allyears_melt[1953:2318]
    abl_1314 = allyears_melt[2319:2684]
    abl_1415 = allyears_melt[2685:3050]
    abl_1516 = allyears_melt[3051:3417]
    abl_1617 = allyears_melt[3418:3783]
    abl_1718 = allyears_melt[3784:4149]
    abl_18 = allyears_melt[4150:] #only 81 days
    
    mb_0708 = np.cumsum(np.add(abl_0708,acc_0708))
    mb_0809 = np.cumsum(np.add(abl_0809,acc_0809))
    mb_0910 = np.cumsum(np.add(abl_0910,acc_0910))
    mb_1011 = np.cumsum(np.add(abl_1011,acc_1011))
    mb_1112 = np.cumsum(np.add(abl_1112,acc_1112))
    mb_1213 = np.cumsum(np.add(abl_1213,acc_1213))
    mb_1314 = np.cumsum(np.add(abl_1314,acc_1314))
    mb_1415 = np.cumsum(np.add(abl_1415,acc_1415))
    mb_1516 = np.cumsum(np.add(abl_1516,acc_1516))
    mb_1617 = np.cumsum(np.add(abl_1617,acc_1617))
    mb_1718 = np.cumsum(np.add(abl_1718,acc_1718))
    mb_18 = np.cumsum(np.add(abl_18,acc_18)) #only 81 days
    
    
    plt.figure(figsize=(12,12))
    
    ax0 = plt.subplot(4,3,1)
    ax1 = ax0.twinx()
    
    ax2 = plt.subplot(4,3,2)
    ax3 = ax2.twinx()
    
    ax4 = plt.subplot(4,3,3)
    ax5 = ax4.twinx()
    
    ax6 = plt.subplot(4,3,4)
    ax7 = ax6.twinx()
    
    ax8 = plt.subplot(4,3,5)
    ax9 = ax8.twinx()
    
    ax10 = plt.subplot(4,3,6)
    ax11 = ax10.twinx()
    
    ax12 = plt.subplot(4,3,7)
    ax13 = ax12.twinx()
    
    ax14 = plt.subplot(4,3,8)
    ax15 = ax14.twinx()
    
    ax16 = plt.subplot(4,3,9)
    ax17 = ax16.twinx()
    
    ax18 = plt.subplot(4,3,10)
    ax19 = ax18.twinx()
    
    ax20 = plt.subplot(4,3,11)
    ax21 = ax20.twinx()
    
    ax22 = plt.subplot(4,3,12)
    ax23 = ax22.twinx()

    ax1.get_shared_y_axes().join(ax1, ax3, ax5, ax7, ax9, ax11, ax13, ax15, ax17, ax19, ax21, ax23)
    ax0.get_shared_y_axes().join(ax0, ax2, ax4, ax6, ax8, ax10, ax12, ax14, ax16, ax18, ax20, ax22)
    ax0.set_ylim([-0.06,0.06])
    ax1.set_ylim([-2.4, 2.4])
    
    plt.suptitle('Mass Balance \n' + str(debristitle), y=1.05)
    
    year = '2007-2008'
    ax0.set_title('Mass Balance: ' + str(year))
    ax0.plot(abl_0708,'r', label ='Ablation')
    ax0.plot(acc_0708, 'b', label = 'Accumulation')
    ax1.plot(mb_0708,'k', label = 'Cumulative \n Mass Balance' )
    #ax0.set_xlabel('Day of Year')
    ax0.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax1.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax0.margins(x=0)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax0.grid(which='both')
    #ax1.legend(loc='upper right', frameon=False)
    ax0.legend(loc='lower left', frameon=False)
    ax1.legend(loc='upper right', frameon=False)
    year_starts = [0,31,61,92,122,153,183,214,244,275,305,336]
    year_names = ['May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr'] 
    ax0.set_xticks(year_starts)
    ax0.set_xticklabels(year_names,rotation=90)
    
    year = '2008-2009'
    ax2.set_title('Mass Balance: ' + str(year))
    ax2.plot(abl_0809,'r', label ='Ablation')
    ax2.plot(acc_0809, 'b', label = 'Accumulation')
    ax3.plot(mb_0809,'k', label = 'Cumulative \n Mass Balance' )
    #ax2.set_xlabel('Day of Year')
    ax2.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax3.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax2.margins(x=0)
    ax2.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax3.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax2.grid(which='both')
    ax2.legend(loc='lower left', frameon=False)
    ax3.legend(loc='upper right', frameon=False)
    ax2.set_xticks(year_starts)
    ax2.set_xticklabels(year_names,rotation=90)
    
    year = '2009-2010'
    ax4.set_title('Mass Balance: ' + str(year))
    ax4.plot(abl_0910,'r', label ='Ablation')
    ax4.plot(acc_0910, 'b', label = 'Accumulation')
    ax5.plot(mb_0910,'k', label = 'Cumulative \n Mass Balance' )
    #ax4.set_xlabel('Day of Year')
    ax4.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax5.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax4.margins(x=0)
    ax4.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax5.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax4.grid(which='both')
    ax4.legend(loc='lower left', frameon=False)
    ax5.legend(loc='upper right', frameon=False)
    ax4.set_xticks(year_starts)
    ax4.set_xticklabels(year_names,rotation=90)
    
    year = '2010-2011'
    ax6.set_title('Mass Balance: ' + str(year))
    ax6.plot(abl_1011,'r', label ='Ablation')
    ax6.plot(acc_1011, 'b', label = 'Accumulation')
    ax7.plot(mb_1011,'k', label = 'Cumulative \n Mass Balance' )
    #ax6.set_xlabel('Day of Year')
    ax6.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax7.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax6.margins(x=0)
    ax6.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax7.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax6.grid(which='both')
    ax6.legend(loc='lower left', frameon=False)
    ax7.legend(loc='upper right', frameon=False)
    ax6.set_xticks(year_starts)
    ax6.set_xticklabels(year_names,rotation=90)
    
    year = '2011-2012'
    ax8.set_title('Mass Balance: ' + str(year))
    ax8.plot(abl_1112,'r', label ='Ablation')
    ax8.plot(acc_1112, 'b', label = 'Accumulation')
    ax9.plot(mb_1112,'k', label = 'Cumulative \n Mass Balance' )
    #ax8.set_xlabel('Day of Year')
    ax8.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax9.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax8.margins(x=0)
    ax8.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax9.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax8.grid(which='both')
    ax8.legend(loc='lower left', frameon=False)
    ax9.legend(loc='upper right', frameon=False)
    ax8.set_xticks(year_starts)
    ax8.set_xticklabels(year_names,rotation=90)
    
    year = '2012-2013'
    ax10.set_title('Mass Balance: ' + str(year))
    ax10.plot(abl_1213,'r', label ='Ablation')
    ax10.plot(acc_1213, 'b', label = 'Accumulation')
    ax11.plot(mb_1213,'k', label = 'Cumulative \n Mass Balance' )
    #ax10.set_xlabel('Day of Year')
    ax10.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax11.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax10.margins(x=0)
    ax10.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax11.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax10.grid(which='both')
    ax10.legend(loc='lower left', frameon=False)
    ax11.legend(loc='upper right', frameon=False)
    ax10.set_xticks(year_starts)
    ax10.set_xticklabels(year_names,rotation=90)
    
    year = '2013-2014'
    ax12.set_title('Mass Balance: ' + str(year))
    ax12.plot(abl_1314,'r', label ='Ablation')
    ax12.plot(acc_1314, 'b', label = 'Accumulation')
    ax13.plot(mb_1314,'k', label = 'Cumulative \n Mass Balance' )
    #ax12.set_xlabel('Day of Year')
    ax12.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax13.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax12.margins(x=0)
    ax12.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax13.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax12.grid(which='both')
    ax12.legend(loc='lower left', frameon=False)
    ax13.legend(loc='upper right', frameon=False)
    ax12.set_xticks(year_starts)
    ax12.set_xticklabels(year_names,rotation=90)
    
    year = '2014-2015'
    ax14.set_title('Mass Balance: ' + str(year))
    ax14.plot(abl_1415,'r', label ='Ablation')
    ax14.plot(acc_1415, 'b', label = 'Accumulation')
    ax15.plot(mb_1415,'k', label = 'Cumulative \n Mass Balance' )
    #ax14.set_xlabel('Day of Year')
    ax14.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax15.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax14.margins(x=0)
    ax14.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax15.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax14.grid(which='both')
    ax14.legend(loc='lower left', frameon=False)
    ax15.legend(loc='upper right', frameon=False)
    ax14.set_xticks(year_starts)
    ax14.set_xticklabels(year_names,rotation=90)
    
    year = '2015-2016'
    ax16.set_title('Mass Balance: ' + str(year))
    ax16.plot(abl_1516,'r', label ='Ablation')
    ax16.plot(acc_1516, 'b', label = 'Accumulation')
    ax17.plot(mb_1516,'k', label = 'Cumulative \n Mass Balance' )
    #ax16.set_xlabel('Day of Year')
    ax16.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax17.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax16.margins(x=0)
    ax16.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax17.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax16.grid(which='both')
    ax16.legend(loc='lower left', frameon=False)
    ax17.legend(loc='upper right', frameon=False)
    ax16.set_xticks(year_starts)
    ax16.set_xticklabels(year_names,rotation=90)
    
    year = '2016-2017'
    ax18.set_title('Mass Balance: ' + str(year))
    ax18.plot(abl_1617,'r', label ='Ablation')
    ax18.plot(acc_1617, 'b', label = 'Accumulation')
    ax19.plot(mb_1617,'k', label = 'Cumulative \n Mass Balance' )
    #ax18.set_xlabel('Day of Year')
    ax18.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax19.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax18.margins(x=0)
    ax18.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax19.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax18.grid(which='both')
    ax18.legend(loc='lower left', frameon=False)
    ax19.legend(loc='upper right', frameon=False)
    ax18.set_xticks(year_starts)
    ax18.set_xticklabels(year_names,rotation=90)

    year = '2017-2018'
    ax20.set_title('Mass Balance: ' + str(year))
    ax20.plot(abl_1718,'r', label ='Ablation')
    ax20.plot(acc_1718, 'b', label = 'Accumulation')
    ax21.plot(mb_1718,'k', label = 'Cumulative \n Mass Balance' )
    #ax20.set_xlabel('Day of Year')
    ax20.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax21.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax20.margins(x=0)
    ax20.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax21.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax20.grid(which='both')
    ax20.legend(loc='lower left', frameon=False)
    ax21.legend(loc='upper right', frameon=False)
    ax20.set_xticks(year_starts)
    ax20.set_xticklabels(year_names,rotation=90)
    
    year = 2018
    ax22.set_title('Mass Balance: ' + str(year))
    ax22.plot(abl_18,'r', label ='Ablation')
    ax22.plot(acc_18, 'b', label = 'Accumulation')
    ax23.plot(mb_18,'k', label = 'Cumulative \n Mass Balance' )
    ax22.set_xlabel('Day of Year')
    ax22.set_ylabel('Daily Mass Change (m w.e. $day^{-1}$)', color='k')
    ax23.set_ylabel('Cumulative Mass Balance (m w.e.)', color='k')
    ax22.margins(x=0)
    ax22.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax23.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(7))
    ax22.grid(which='both')
    ax22.legend(loc='lower left', frameon=False)
    ax23.legend(loc='upper right', frameon=False)
    #ax0.set_xticks(year_starts)
    #ax0.set_xticklabels(year_names,rotation=90)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(Path2files,'MaytoApril_MB_DebrisCase.pdf'),bbox_inches='tight')

else:
    pass

if Plot_mean_3hourly_temp == True:
    print('mean 3 hourly temperature')
    # don't need to do any extracting for this one --> just make sure Path2files leads
    # to the downscaled and bias corrected temp files
    # import temp from each year
    allyears_array = np.empty((12,2920))
    for year in years[1:]:
        i = year-2007
        print(year)
        filename = 'Tempkaskonly_BC_' + str(year) + '.nc'
        yearfile = os.path.join(Path2files,filename)
        T_in = Dataset(yearfile, "r")
        T_var = 'Temperature'
        T_array1 = T_in.variables[T_var][:]
        
        T_array = T_array1 - 1.7
        
        mean_x = np.nanmean(T_array,axis=1)
        mean_xy = np.nanmean(mean_x,axis=1)
        
        allyears_array[i] = mean_xy[:2920]
    
    mean_3hourly_T = np.nanmean(allyears_array,axis=0)
    

    plt.figure(figsize=(8,8))
    ax0 = plt.subplot()
    plt.title('Kaskawulsh Mean Temperature (2007-2018)',fontsize=15)
    ax0.plot(mean_3hourly_T,'black',zorder=5)
    ax0.set_ylabel('Temperature (C)', color='k',fontsize=15)
    ax0.set_xlabel('Day Of Year',fontsize=15)
    #plt.legend(['2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','Mean'])
    month_starts = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    month_names = ['Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'] 
    plt.margins(x=0)
    ax0.margins(x=0)
    plt.setp(ax0.get_xticklabels(), fontsize=15)
    #ax0.set_ylim(0,340)
    ax0.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(6))
    ax0.set_frame_on(False)
    #plt.savefig(os.path.join(Path2files,'KaskawulshMeanTemp.pdf'),bbox_inches = 'tight')
    
        
else:
    pass
      

a = np.empty((4,3))
a[0] = [1,2,3]
a[1] = [4,5,6]
a[2] = [7,8,9]
