# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 17:22:10 2022

comparing the glacier outlines from kaskonly.txt vs kaskonly_deb.txt

@author: katierobinson
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import sys
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something

File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kaskonly.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,5]       
        
Zgrid1, Xgrid1, Ygrid1, xbounds1, ybounds1 = regridXY_something(Ix, Iy, Ih)

File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kaskonly_deb.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,2]       
        
Zgrid2, Xgrid2, Ygrid2, xbounds2, ybounds2 = regridXY_something(Ix, Iy, Ih)


File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kask_catchment.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,4] 
Iy = glacier[:,5] 
Ih = glacier[:,6]  
sfc_type = glacier[:,8]        

Zgrid3, Xgrid3, Ygrid3, xbounds3, ybounds3 = regridXY_something(Ix, Iy, Ih)
Sfc_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, sfc_type)

plt.figure(figsize=(7.2,12.6))
plt.subplot(3,1,1)
plt.title('Outline from kaskonly.txt',fontsize=14)
plt.contourf(np.flipud(Zgrid1), cmap = 'GnBu', levels = np.linspace(700,3700,16))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.subplot(3,1,2)
plt.title('Outline from kaskonly_deb.txt',fontsize=14)
plt.contourf(np.flipud(Zgrid2), cmap = 'GnBu', levels = np.linspace(700,3700,16))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.subplot(3,1,3)
plt.title('Outline from kask_catchment.txt',fontsize=14)
plt.contourf(np.flipud(Zgrid3), cmap = 'GnBu', levels = np.linspace(700,3700,16))
#plt.contourf(np.flipud(Zgrid2), cmap = 'bone', levels = np.linspace(700,3700,16))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.contour(np.flipud(Sfc_grid), colors = 'black', levels = 0, linewidth = 0.2)
plt.tight_layout()
#plt.savefig('ComparingGlacierOutlines.png',bbox_inches = 'tight')

plt.figure(figsize=(9,5))
plt.title('Outline from kaskonly_deb.txt',fontsize=14)
plt.contourf(np.flipud(Zgrid2), cmap = 'GnBu', levels = np.linspace(700,3700,16))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.savefig('KWoutline.png',bbox_inches = 'tight')