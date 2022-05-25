# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:26:33 2022

script to process raw NARR files (cut down to KW domain)

@author: katierobinson
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

infile = 'F:\\Mass Balance Model\\Kaskawulsh-Mass-Balance\\Downscaling\\temp_KW.nc'
inF = Dataset(infile,'r')
Fvar = 'var11'
F_array = inF.variables[Fvar][:] # has shape (time, lev, y, x)

xvar = 'x'
yvar = 'y'
x = inF.variables[xvar][:]
y = inF.variables[yvar][:]

F_array2 = np.nanmean(F_array,axis=0)
F_array3 = np.nanmean(F_array2,axis=0)


infile = 'F:\\Mass Balance Model\\Kaskawulsh-Mass-Balance\\Downscaling\\temp_KW.nc'
inF = Dataset(infile,'r')
Fvar = 'var11'
F_array = inF.variables[Fvar][:] # has shape (time, lev, y, x)

# open KW_domain_coords, get x y data for selected region
# shift x y data to the 0 - end scale
# cut that domain out of full NARR file
#then, write a function to do that automatically

# range of x and y vals from raw narr files: 
xrange = np.linspace(-5629.3*1000,5667.8*1000,349) #got these vals from loading raw narr file into PANOPLY and looking at x, y coordinates
yrange = np.linspace(-4609.8*1000,4349.9*1000,277)

#KW domain (delineated in ArcMap)

coordinate_file = 'F:\\Mass Balance Model\\Kaskawulsh-Mass-Balance\\Downscaling\\KW_domain_coords.csv'
domain = np.loadtxt(coordinate_file,delimiter=",",skiprows=1,dtype='float')

x = domain[:,4]
y = domain[:,5]