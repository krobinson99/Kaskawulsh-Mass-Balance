# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 10:34:32 2022

@author: katierobinson
"""
#the purpose of this script is to covert .TIF data to a numpy array
# for converting Rounce et al. (2020) debris thickness data for the Kaskawulsh to a useable data type

import numpy as np
import PIL
from PIL import Image
import matplotlib.pyplot as plt
import os
#from Model_functions_ver4 import regridXY_something

Path2files = 'F:/Mass Balance Model/DebrisThickness'

debristhickness = Image.open('F:/Mass Balance Model/DebrisThickness/HMA_DTE_1.16201_hdts_m.tif')
#debristhickness.show()
debristhickness_array = np.array(debristhickness)
debristhickness_array.shape

nanlocs = np.where(debristhickness_array >= 1e20 )
#zerolocs = np.where(debristhickness_array == 0 )
debristhickness_array[nanlocs] = np.nan
#debristhickness_array[zerolocs] = np.nan

plt.figure(figsize=(8,5))
#plt.contourf(debristhickness_array, cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
plt.contourf(np.flipud(debristhickness_array[:300,200:]), cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Debris Thickness (m)', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
plt.title('Kaskawulsh Debris Thickness Estimate',fontsize=14)
plt.savefig(os.path.join(Path2files,'KW_debristhickness.png'),bbox_inches = 'tight')


plt.figure(figsize=(8,5))
plt.contourf(np.flipud(debristhickness_array[:300,200:]), cmap = 'PuOr', levels = np.round(np.linspace(0,1,10),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Debris Thickness (m)', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
plt.title('Kaskawulsh Debris Thickness Estimate',fontsize=14)
plt.savefig(os.path.join(Path2files,'KW_debristhickness_zoomed.png'),bbox_inches = 'tight')

#max thickness = 3m, min thickness = 0m, mean thickness = 0.18 m

#now do the same thing with the melt factors
meltfactors = Image.open('F:/Mass Balance Model/DebrisThickness/HMA_DTE_1.16201_meltfactor.tif')
#meltfactors.show()
MF_array = np.array(meltfactors)
MF_array.shape

nanlocs2 = np.where(MF_array >= 1e20 )
MF_array[nanlocs2] = np.nan

plt.figure(figsize=(8,5))
#plt.contourf(debristhickness_array, cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
plt.contourf(np.flipud(MF_array[:300,200:]), cmap = 'coolwarm', levels = np.round(np.linspace(0,2,9),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Sub-Debris Melt Factor', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
plt.title('Kaskawulsh Sub-Debris Melt Enhancement Factors',fontsize=14)
plt.savefig(os.path.join(Path2files,'KW_meltfactors.png'),bbox_inches = 'tight')

####Plot debris mask from the original model
File_glacier_in = os.path.join(Path2files,'Kaskonly_deb.txt')
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,2] 
debris_array = glacier[:,6]

#-------Turn vectors into 3D gridded inputs--------------------

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs3 = np.where(np.isnan(Zgrid))

#Setup debris mask for use in radiation parameters
debris_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, debris_array)
debris_m = np.zeros(debris_grid.shape)
debris_m[np.where(debris_grid > 100)] = 1
debris_m[np.where(debris_grid <= 100)] = np.nan
debris_m[nanlocs3] = np.nan

plt.figure(figsize=(8,5))
#plt.contourf(debristhickness_array, cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
plt.contourf(np.flipud(debris_m[:180,100:]), cmap = 'PuOr', levels = np.round(np.linspace(0,1,3),1))
#plt.contourf((debris_m[:180,100:]), cmap = 'PuOr', levels = np.round(np.linspace(0,1,3),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Debris True / False', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
plt.title('Kaskawulsh Boolean Debris Map',fontsize=14)
#plt.savefig(os.path.join(Path2files,'KW_booleandebris.png'),bbox_inches = 'tight')



######## plot melt factor vs debris thickkness
meltfacts = []
thickness = []
for x in range(len(MF_array)):
    for y in range(len(MF_array[0])):
        h =  debristhickness_array[x][y]
        MF = MF_array[x][y]
        meltfacts.append(MF)
        thickness.append(h)
        
plt.figure(figsize=(8,6))
plt.scatter(thickness,meltfacts)
plt.xlabel('Debris Thickness (m)',fontsize=14)
plt.ylabel('Melt Enhancement Factor',fontsize=14)
plt.axhline(y=1,color='k',linestyle='--')
plt.axvline(x=0.13,color='k',linestyle='--')
#plt.axvline(x=0.03,color='k',linestyle='--')
plt.title('Transition thickness = 13 cm',fontsize=14)
plt.xlim(0,0.45)
plt.savefig(os.path.join(Path2files,'debristhickness_vs_MF_Ostrem.pdf'),bbox_inches = 'tight')

transitionMF = np.where(MF_array > 0.99 and MF_array < 0.101)
transitionh = debristhickness_array[transitionMF]

transition_thickness = np.nanmean(transitionh)


