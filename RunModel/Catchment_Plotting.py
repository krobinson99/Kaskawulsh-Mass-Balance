# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 20:18:09 2022

@author: katierobinson
this script is meant for plotting figures related to the catchment-wide mass balance model. 
outputs from the model can still be plotting using the PlotOutputs.py script but this script 
will be used to generate figs related to the domain, hypsometry, input data, etc. 
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import colors
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

Glacier_outline = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\kaskonly.txt'
Catchment_outline = 'D:\\Katie\\Mass Balance Model\\MassBalanceModel_KatiesVersion\\Final_runs\\kask_catchment.txt'

debris = False

###############################################################################

# what do you want the script to actually plot
Plot_Catchment_sfctype_elevation = False
Catchment_hypsometric_curve = True

###############################################################################

glacier = np.genfromtxt(Glacier_outline, skip_header=1, delimiter=',')
        
#if debris is TRUE then we are working from the kaskonly_deb.txt file,
#if debris is FALSE then we are working from the kaskonly.txt file:
if debris == True:
    print('kaskonly debris outline')
    Ix = glacier[:,3] 
    Iy = glacier[:,4] 
    Ih = glacier[:,2] 
    debris_array = glacier[:,6]
    #sfc_type = [:,]
else:
    print('kaskonly no debris outline')
    Ix = glacier[:,2]
    Iy = glacier[:,3] 
    Ih = glacier[:,1]
    
catchment = np.genfromtxt(Catchment_outline, skip_header=1, delimiter=',')

print('full catchment outline')
Ix_catch = catchment[:,4]
Iy_catch = catchment[:,5] 
Ih_catch = catchment[:,6]
sfc_type = catchment[:,8]
        
Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix_catch, Iy_catch, Ih_catch)
Sfc_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix_catch, Iy_catch, sfc_type)



if Plot_Catchment_outline == True:
    print('plotting model outline from .txt file')
    
    #PLOT MEAN MASS BALANCE OVER ENTIRE SIMULATION PERIOD: 
    #Sfc_kask = np.zeros(Zgrid.shape)
    #Sfc_kask[nanlocs] = np.nan
    # make a color map of fixed colors
    cmap = colors.ListedColormap(['lightsteelblue', 'saddlebrown'])
    bounds=[0,2]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    
    #plt.figure(figsize=(12,7))
    plt.figure(figsize=(10,7))
    #plt.contourf(np.flipud(Sfc_grid), cmap = 'bwr', levels = np.linspace(0,1, 5))
    plt.contourf(np.flipud(Sfc_grid), cmap = cmap, levels = np.linspace(0,1, 5))
    #plt.contourf(np.flipud(Sfc_grid), cmap = 'RdYlBu', levels = np.linspace(0,1,5))
    #plt.contourf(np.flipud(Sfc_kask), cmap = 'bwr', levels = np.linspace(0,1,5))
    #legend = plt.colorbar(ticks=[0,1])
    #legend.ax.set_ylabel('OFF ice   /   ON ice', rotation=270,labelpad=20)
    #plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = 0)
    plt.title('Kaskawulsh Catchment',fontsize=14)
    plt.xlabel('Easting',fontsize=14)
    plt.ylabel('Northing',fontsize=14)
    plt.xticks(ticks=[0,50,100,150,200,250,300],labels=[int(Xgrid[0][0]),int(Xgrid[0][50]),int(Xgrid[0][100]),int(Xgrid[0][150]),int(Xgrid[0][200]),int(Xgrid[0][250]),int(Xgrid[0][300])],fontsize=12,rotation=20)
    plt.yticks(ticks=[0,50,100,150,200],labels=[int(Ygrid[:,0][0]),int(Ygrid[:,0][50]),int(Ygrid[:,0][100]),int(Ygrid[:,0][150]),int(Ygrid[:,0][200])],fontsize=12,rotation=20)
    plt.savefig(os.path.join(Path2files,'Catchment_Outline.pdf'),bbox_inches = 'tight')
    
    plt.figure(figsize=(12,7))
    plt.contourf(np.flipud(Zgrid), cmap = 'bone', levels = np.linspace(700,3700,16))
    #plt.contourf(np.flipud(Sfc_grid), cmap = 'RdYlBu', levels = np.linspace(0,1,5))
    #plt.contourf(np.flipud(Sfc_kask), cmap = 'bwr', levels = np.linspace(0,1,5))
    legend = plt.colorbar()
    legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
    #plt.contour(np.flipud(np.nanmean(overall_MB, axis = 0)), colors = 'black', levels = 0)
    #plt.title('Kaskawulsh Catchment')
    #plt.savefig(os.path.join(Path2files,'Catchment_Elevation.pdf'),bbox_inches = 'tight')
    
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

#############################################################################3
    
if Catchment_hypsometric_curve == True:
    print('hypsometric curves')
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
                    
    #np.nanmin(Zgrid) is 753 m
    #np.nanmax(Zgrid) is 3923 m
                
    #calculate frequency of each elevation in each 100m bin 
    #make a histogram
    office_histogram = np.histogram(officeZ,33,(700,4000))
    onice_histogram = np.histogram(oniceZ,33,(700,4000))
    
    labels = []
    bin_num = [0]
    #for i in range(0,len(office_histogram[1])-1):
    i = 0
    while i < len(office_histogram[1]):
        #print(office_histogram[1][i])
        print(i)
        mid_bin = int((office_histogram[1][i] + office_histogram[1][i+1])/2)
        i += 2
        #print(mid_bin)
        labels.append(mid_bin)
        bin_num.append(i)
        
        
    #seasonlabels = ['DJF','MAM','JJA','SON'    
    x = np.arange(len(labels))
    width = 0.8
    
    plt.figure(figsize=(10,6))
    plt.subplot(121)
    plt.title('Off-Glacier Hypsometry', fontsize=14)
    plt.barh(x,office_histogram[0],width)
    plt.yticks(bin_num[:-1],labels,fontsize=12)
    plt.margins(y=0)
    plt.ylabel('Elevation (m)',fontsize=14)
    plt.xlabel('Frequency',fontsize=14)
    plt.xlim(0,2800)
    plt.subplot(122)
    plt.title('On-Glacier Hypsometry', fontsize=14)
    plt.barh(x,onice_histogram[0],width)
    plt.yticks(bin_num[:-1],labels,fontsize=12)
    plt.margins(y=0)
    plt.ylabel('Elevation (m)',fontsize=14)
    plt.xlabel('Frequency',fontsize=14)
    plt.xlim(0,2800)
    plt.tight_layout(w_pad=1)
    #plt.savefig(os.path.join(Path2files,'Catchment_Hypsometry.pdf'),bbox_inches = 'tight')
     
    # make another hypsometry figure where boxes are overlayed for easier comparison
    # also scale x-axis from frequency to Area
    # each cell = 200m x 200m = 0.04km^2
    
    plt.figure(figsize=(5,6))
    plt.title('Catchment Hypsometry', fontsize=14)
    plt.barh(x,(0.04*onice_histogram[0]),width,color='lightsteelblue')
    plt.barh(x,(0.04*office_histogram[0]),width,color='saddlebrown')
    plt.yticks(bin_num[:-1],labels,fontsize=12)
    plt.margins(y=0)
    plt.legend(['On-ice cells','Off-ice cells'],fontsize=14)
    plt.ylabel('Elevation (m)',fontsize=14)
    plt.xlabel('Area (km$^2$)',fontsize=14)
    #plt.xlim(0,2800)
    #plt.savefig(os.path.join(Path2files,'Catchment_Hypsometry_Overlay.pdf'),bbox_inches = 'tight')
    
else:
    pass

