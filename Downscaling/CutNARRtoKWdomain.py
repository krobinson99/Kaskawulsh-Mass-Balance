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
NARR = np.nanmean(F_array2,axis=0) #has shape (y,x)

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

xKW = domain[:,4]
yKW = domain[:,5]

plt.plot(xKW)
plt.plot(xrange)

plt.plot(yKW)
plt.plot(yrange)

#####################################################
KWdomain_xleft = np.where(xrange > np.min(xKW)+(-15000))[0][0] #add/subtract half a cell (15km) on either side to catch any cells that are inside the bounding box but have centres outside
KWdomain_xright = np.where(xrange < np.max(xKW)+15000)[0][-1]

KWdomain_ybott = np.where(yrange > np.min(yKW)-15000)[0][0]
KWdomain_ytop = np.where(yrange < np.max(yKW)+15000)[0][-1]
# gives 6 x 6 grid

# COORDS TO TAKE FROM FULL RAW NARR DOMAIN: X[118-123], Y[188-193]
KWDOMAIN = NARR[188:194,118:124]
X_KW = x[188:194]
Y_KW = y[118:124]

#NARR[188:194,118:124] = 305


#write KWDOMAIN to a netcdf file so  it can be opened in panoply:
f = Dataset('KWdomain.nc','w',format='NETCDF4')
datagroup = f.createGroup('Data')
datagroup.createDimension('x',6)
datagroup.createDimension('y',6)

xvar = datagroup.createVariable('easting','f4','x')
yvar = datagroup.createVariable('northing','f4','y')
temp = datagroup.createVariable('Temp','f4',('y','x'))

print(f)
print(f.groups['Data'])

xvar[:] = X_KW
yvar[:] = Y_KW
temp[:] = KWDOMAIN
f.close()
