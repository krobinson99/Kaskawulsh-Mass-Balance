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
import sys
import matplotlib.pyplot as plt
import random
#from osgeo import gdal

sys.path.insert(1,'D:\Katie\Mass Balance Model\MassBalanceModel_KatiesVersion\RunModel')
from Model_functions_ver4 import regridXY_something


# for now this script will use a "fake" debris thickness map, until I can
# figure out how to turn the 100m res. debris map into a 200m one that matches
# the model domain

def LoadRounceDebrisMap(DebrisMapPath):
    """
    this function takes the .tif files from Rounce et al. (2021) as input and
    returns a numpy array.
    """
    debristhickness = Image.open(DebrisMapPath)
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

def MeltFactors(debris_array,tt,Mcleanice=2.9200183038347,M2cm=6.62306028549649,b0=11.0349260206858,k=1.98717418666925):
    """
    this function calculates the melt enhacement field using a user-specified transition thickness
    (tt) in meters. The tt is where the melt changes from enhanced to reduced.
    """
    # don't use clean ice melt at all, use tt_melt
    tt_melt = b0/(1+(k*b0*tt))
    
    Melt_array = np.zeros(debris_array.shape)
    Melt_array = (((M2cm-tt_melt)/0.02)*debris_array) + tt_melt #linear portion of melt curve for debris 0-2cm
    
    thickerthan2cm = np.where(debris_array > 0.02)
    Melt_array[thickerthan2cm] = b0/(1+(k*b0*debris_array[thickerthan2cm]))
    # cap melt at the value of M2cm ??
    meltcap = np.where(Melt_array > M2cm) # comment out these lines if it should NOT be capped
    Melt_array[meltcap] = M2cm 
    
    MeltFacts = Melt_array/tt_melt
    
    #return Melt_array
    return MeltFacts

def MeltFactorsv2(debris_array,tt,Mcleanice=2.9200183038347,M2cm=6.62306028549649,b0=11.0349260206858,k=1.98717418666925):
    """
    THIS FUNCTION IS NOW IN THE DEBRISFUNCTIONS.PY FILE. EDIT THERE INSTEAD OF HERE.
    this function calculates the melt enhacement field using a user-specified transition thickness
    (tt) in meters. The tt is where the melt changes from enhanced to reduced.
    """
    # don't use clean ice melt at all, use tt_melt
    tt_melt = b0/(1+(k*b0*tt))
    
    Melt_array = np.zeros(debris_array.shape)
    Melt_array = (((M2cm-tt_melt)/0.02)*debris_array) + tt_melt #linear portion of melt curve for debris 0-2cm
    
    thickerthan2cm = np.where(debris_array > 0.02)
    Melt_array[thickerthan2cm] = b0/(1+(k*b0*debris_array[thickerthan2cm]))
    # cap melt at the value of M2cm ??
    #meltcap = np.where(Melt_array > M2cm) # comment out these lines if it should NOT be capped
    #Melt_array[meltcap] = M2cm 
    
    MeltFactors = Melt_array/tt_melt
    
    # make a copy of melt factors to get the indices where MF is incorrect (if any)
    MeltFactors_copy = np.empty(MeltFactors.shape)
    MeltFactors_copy[:] = MeltFactors[:]
    MeltFactors_copy[thickerthan2cm] = np.nan
    
    incorrectMFs = np.where(MeltFactors_copy < 1)
    MeltFactors[incorrectMFs] = 1
    
    #return Melt_array
    return MeltFactors

def SyntheticDebrisMap(File_glacier_in):
    """
    this function generates a "fake" debris thickness map that is the same size as the
    mass-balance model domain. This map can then be used to generate the melt factors
    to input into the model.
    
    Will not need to use this function once the debris map from Rounce et al. (2021) has been
    projected onto the coordinates of the model domain
    """
    #generate model domain with 1's where debris should be:
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
    
    for i in range(len(debris_m)):
        for j in range(len(debris_m[0])):
            if debris_m[i,j] == 1:
                debris_m[i,j] = random.uniform(0,0.5)
            else:
                pass
                
    return debris_m
