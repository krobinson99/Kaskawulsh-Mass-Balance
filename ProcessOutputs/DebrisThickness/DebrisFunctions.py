# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:33:21 2022

This code is for generating a sub-debris melt enhancement factor map, based on some 
user-defined transition debris-thickness.
The melt-enhancement map will then be used as a model input

@author: katierobinson
"""
import numpy as np
import PIL
from PIL import Image
import os
import matplotlib.pyplot as plt
#from osgeo import gdal

#sys.path.insert(1,'D:\Katie\Mass Balance Model\MassBalanceModel_KatiesVersion\RunModel')
#from Model_functions_ver4 import regridXY_something


# for now this script will use a "fake" debris thickness map, until I can
# figure out how to turn the 100m res. debris map into a 200m one that matches
# the model domain

def DebrisMap(Path2files='D:\Katie\Mass Balance Model\MassBalanceModel_KatiesVersion\ProcessOutputs\DebrisThickness'):
    """
    this function takes the .tif files from Rounce et al. (2021) as input and
    returns a numpy array.
    """
    thicknessmap = os.path.join(Path2files,'HMA_DTE_1.16201_hdts_m.tif')
    debristhickness = Image.open(thicknessmap)
    debristhickness_array = np.array(debristhickness)
    debristhickness_array.shape

    nanlocs = np.where(debristhickness_array >= 1e20 )
    debristhickness_array[nanlocs] = np.nan
    
    return debristhickness_array
    
def MeltEnhancementMap(Path2files='D:\Katie\Mass Balance Model\MassBalanceModel_KatiesVersion\ProcessOutputs\DebrisThickness'):
    """
    this function takes the .tif files from Rounce et al. (2021) as input and
    returns a numpy array.
    """
    mf_map = os.path.join(Path2files,'HMA_DTE_1.16201_meltfactor.tif')
    meltfactors = Image.open(mf_map)
    #meltfactors.show()
    MF_array = np.array(meltfactors)
    MF_array.shape
    
    nanlocs2 = np.where(MF_array >= 1e20 )
    MF_array[nanlocs2] = np.nan

    return MF_array    

def MeltFactors(debris_array,tt,Mcleanice=2.9200183038347,M2cm=6.62306028549649,bo=11.0349260206858,k=1.98717418666925):
    """
    this function calculates the melt enhacement field using a user-specified transition thickness
    (tt) in meters. The tt is where the melt changes from enhanced to reduced
    """
    SubdebrisMelt = M2cm/(1+(k*M2cm*debris_array)) # eq(1) from Rounce et al. (2021)
    tt_melt = M2cm/(1+(k*M2cm*tt)) #calculate the melt at the transition thickness using eq(1) from Rounce et al (2021)
    MEF = SubdebrisMelt/tt_melt # calculate the melt factors by dividing the melt array by the melt at the transition thickness
    
    return MEF

deb = DebrisMap()
mef = MeltFactors(deb,0.05)

debthickness = []
MEF = []
for x in range(len(deb)):
    for y in range(len(deb[0])):
        h =  deb[x][y]
        mefs =  mef[x][y]
        debthickness.append(h)
        MEF.append(mef)

#plt.figure(figsize=(8,6))
#plt.scatter(debthickness,MEF)
#plt.xlabel('Debris Thickness (m)',fontsize=14)
#plt.ylabel('Melt Enhancement Factor (m w.e. a$^{-1}$)',fontsize=14)
#plt.axhline(y=1,color='k',linestyle='--')
#plt.axvline(x=0.03,color='k',linestyle='--')
#plt.axvline(x=0.03,color='k',linestyle='--')
#plt.title('Transition Thickness = 3cm \n M$_0$ = 6.62 \n M$_{clean ice}$ = M$_{(h=3cm)}$ = 4.74',fontsize=14)
#plt.xlim(0,0.5)