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
from Model_functions_ver4 import KWtributaries
from Model_functions_ver4 import model_domain

File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kaskonly.txt'
##-------Turn glacier grid vectors into 3D grids--------------------##
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
        
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,5]   
ELA = glacier[:,9]
TOPO = glacier[:,7]    
        
Zgrid1, Xgrid1, Ygrid1, xbounds1, ybounds1 = regridXY_something(Ix, Iy, Ih)
ELA_grid, Xgrid1, Ygrid1, xbounds1, ybounds1 = regridXY_something(Ix, Iy, ELA)
TOPO_grid, Xgrid1, Ygrid1, xbounds1, ybounds1 = regridXY_something(Ix, Iy, TOPO)

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
ELA = glacier[:,10]     

Zgrid3, Xgrid3, Ygrid3, xbounds3, ybounds3 = regridXY_something(Ix, Iy, Ih)
Sfc_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, sfc_type)
#ELA_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, ELA)

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

plt.figure(figsize=(10,5))
#plt.title('Outline from kaskonly_deb.txt',fontsize=14)
plt.contourf(Xgrid2,Ygrid2,np.flipud(Zgrid2), cmap = 'GnBu', levels = np.linspace(700,3600,30))
plt.axis('equal')
legend = plt.colorbar()
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.text(590000,6707000,'____________')
plt.text(590800,6707000+500,'10 km',fontsize=14)
plt.savefig('KW-elevation.pdf',bbox_inches = 'tight')

plt.figure(figsize=(12,7))
plt.title('Outline from kaskonly.txt',fontsize=14)
plt.contourf(np.flipud(Zgrid2), cmap = 'GnBu', levels = np.linspace(700,3700,16))
legend = plt.colorbar()
legend.ax.set_ylabel('Elevation (m a.s.l.)', rotation=270,fontsize=14,labelpad=20)
plt.contourf(np.flipud(Zgrid1), cmap = 'OrRd', levels = np.linspace(700,3700,16))
legend = plt.colorbar()

plt.figure(figsize=(9,5))
plt.contourf(np.flipud(Sfc_grid))
legend = plt.colorbar()


# =============================================================================
# GET TRIBUTARY ARRAY FOR THE CATCHMENT GRID:
# =============================================================================
Zgrid_KW, Xgrid_KW, Ygrid_KW, xbounds, ybounds, Sfc_KW = model_domain(catchment=False)
Zgrid, Xgrid, Ygrid, xbounds, ybounds, Sfc = model_domain(catchment=True)
np.savetxt('KRH_SfcType.txt',Sfc)

tribarray = KWtributaries()[0]
nanlocs_KW = np.where(np.isnan(Zgrid_KW))
tribarray[nanlocs_KW] = np.nan

plt.figure(figsize=(8,6))
plt.contourf(Xgrid,np.flipud(Ygrid),Sfc)
plt.colorbar()
plt.contourf(Xgrid_KW,np.flipud(Ygrid_KW),tribarray,cmap='Accent',levels=np.arange(0,6))
#plt.colorbar()
plt.axis('equal')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)

Xgrid_KW_int = np.array(Xgrid_KW,dtype=int)
Xgrid_int = np.array(Xgrid,dtype=int)
Ygrid_KW_int = np.array(Ygrid_KW,dtype=int)
Ygrid_int = np.array(Ygrid,dtype=int)

np.where(Xgrid_int[0] == Xgrid_KW_int[0,0])
np.where(Xgrid_int[0] == Xgrid_KW_int[0,-1])

np.where(Ygrid_int[:,0] == Ygrid_KW_int[0,0])
np.where(Ygrid_int[:,0] == Ygrid_KW_int[-1,0])

# Test:
newY = np.array(Ygrid_int)
newY[0:217+1,1:328+1] - Ygrid_KW_int # all zeros! that means its a match
newX = np.array(Xgrid_int)
newX[0:217+1,1:328+1] - Xgrid_KW_int # all zeros! that means its a match

# =============================================================================
#  sub the tribarray into the catchment Sfc grid.
# =============================================================================
catchment_tribarray = np.zeros(Zgrid.shape)
catchment_tribarray[:,:] = np.nan
# 0 = small glaciers
# -1 = offglacier
# SA = 1
# SW = 2
# CA = 3
# NA = 4
# Trunk = 5
#catchment_tribarray[0:217+1,1:328+1] = tribarray
catchment_tribarray[12:229+1,1:328+1] = tribarray

plt.figure(figsize=(8,6))
plt.contourf(Xgrid,np.flipud(Ygrid),Zgrid,cmap='Greys_r')
plt.colorbar()
#plt.contourf(Xgrid_KW,np.flipud(Ygrid_KW),tribarray,cmap='Accent',levels=np.arange(0,6))
plt.contourf(Xgrid,np.flipud(Ygrid),catchment_tribarray)
plt.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0,colors='k',linestyles='dashed')
#plt.colorbar()
plt.axis('equal')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)

#np.savetxt('KRH_Tributaries.txt',catchment_tribarray)

# =============================================================================
# Sub the balance flux map onto the catchment Sfc grid
# =============================================================================
KW_fluxgates = KWtributaries()[1]

catchment_fluxarray = np.zeros(Zgrid.shape)
catchment_fluxarray[:,:] = np.nan
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
#catchment_tribarray[0:217+1,1:328+1] = tribarray
catchment_fluxarray[12:229+1,1:328+1] = KW_fluxgates

plt.figure(figsize=(8,6))
plt.contourf(Xgrid,np.flipud(Ygrid),Zgrid,cmap='Greys_r')
plt.colorbar()
#plt.contourf(Xgrid_KW,np.flipud(Ygrid_KW),tribarray,cmap='Accent',levels=np.arange(0,6))
plt.contourf(Xgrid,np.flipud(Ygrid),catchment_fluxarray)
plt.contour(Xgrid,np.flipud(Ygrid),Sfc,levels=0,colors='k',linestyles='dashed')
#plt.colorbar()
plt.axis('equal')
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)

#np.savetxt('KRH_Fluxgates.txt',catchment_fluxarray)
