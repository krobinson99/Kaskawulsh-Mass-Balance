# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 15:25:49 2023

1. Convert DEM's to aspect, slope grids
2. Save DEM, aspect, and slope as txt files in the correct format for HOCK model
3. scp files to cedar and run solar model
4. scp solar outputs back here
5. Plot solar results! compare different years

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import os
import pandas as pd

# Where to save the DEM/slope/aspect inputs generated by this script:
path_to_files = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Solar/Hock_model_inputs'

# Note: Only need the next two funcntions if you do not already have a georeferenced DEM (.tif) to start with:
def read_geotiff(filename):
    '''
    Function that reads in a geotiff, and returns an array + a gdal datset containing the geospatial information
    
    filename: path to georeferenced tiff file
    returns ds: a dataset containing the geospatial information from the tiff file. 
    '''
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr, ds

def write_geotiff(filename, arr, in_ds):
    '''
    Convert a numpy array to a georeferenced tif file (given the geospatial info from another tif file)
    
    filename: output path for newly generated files
    arr: numpy array to be converted to .tif
    in_ds: geospatial dataset, generated using the read_geotiff function
    '''
    if arr.dtype == np.float32:
        arr_type = gdal.GDT_Float32
    else:
        arr_type = gdal.GDT_Int32

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(filename, arr.shape[1], arr.shape[0], 1, arr_type)
    # Get geospatial information from input dataset
    out_ds.SetProjection(in_ds.GetProjection()) 
    out_ds.SetGeoTransform(in_ds.GetGeoTransform())
    band = out_ds.GetRasterBand(1)
    # Write the input array to file
    band.WriteArray(arr)
    band.FlushCache()
    band.ComputeStatistics(False)
    
# Get the geospatial information from tif file sent by E.B.
EtienneBerthier_DEM, DEM_georeference = read_geotiff('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/Coregistered_DEMs/gdal/2018DEM_coregistered_regridded_average.tif')

# =============================================================================
# GENERATE INPUT FILES FOR SHADING MODEL (DEM, ASPECT, SLOPE)
# =============================================================================

dem = np.load('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/SurfaceZ/DEM2018_gapfilled.npy')
year = 2018    

# Convert DEM array to georeferenced tif:
# Note: skip this step if you already have a DEM georeferenced tif file!!
fname = 'KRH_DEM_' + str(year) + '.tif'
DEMpath = os.path.join(path_to_files,fname)
write_geotiff(DEMpath,dem,DEM_georeference)
    
# Open geotiff DEM with GDAL
dem_dataset = gdal.Open(DEMpath)

# Calculate aspect on the DEM grid:
asp = gdal.DEMProcessing(os.path.join(path_to_files,'KRH_Aspect_' + str(year) + '.tif'), dem_dataset, "aspect",  computeEdges=True) 
aspect_array = asp.GetRasterBand(1).ReadAsArray()
aspect_array[np.where(aspect_array==-9999)] = 0

# Calculate slope on the DEM grid:
slp = gdal.DEMProcessing(os.path.join(path_to_files,'KRH_Slope_' + str(year) + '.tif'), dem_dataset, "slope",  computeEdges=True)
slope_array = slp.GetRasterBand(1).ReadAsArray()

# Calculate hillshade on the DEM grid:
hill = gdal.DEMProcessing(os.path.join(path_to_files,'KRH_Hillshade_' + str(year) + '.tif'), dem_dataset, "hillshade",  computeEdges=True)
hillshade_array = hill.GetRasterBand(1).ReadAsArray()

# Close datsets
slp = asp = hill = dem_dataset = None

# =============================================================================
#  Convert slope/aspect/dem to text files (formatted for input to Hock solar model)
# =============================================================================
# ncols, nrows, nodata must be integers (see thesis appendix B for details on how to format)
# Note: If you run this code multiple times it will just append lines to the end of the .txt files that were already generated
# Need to delete all .txt files each time to create fresh ones

lowerleft_x = 567282.84286 # x coordinate of lower left corner of grid (UTM)
lowerleft_y = 6702785.56577 # y coordinate of lower left corner of grid (UTM)
gridcell_size = 200 # grid cell size in meters

data = [['ncols',int(dem.shape[1])],['nrows',int(dem.shape[0])],['xllcorner',lowerleft_x],['yllcorner',lowerleft_y],['cellsize',gridcell_size],['nodata_value',-9999]]
df = pd.DataFrame(data,columns=['name','number'])

with open(os.path.join(path_to_files,'KRH_DEM_' + str(year) + '.txt'),'a') as f:
    f.seek(0)
    np.savetxt(f,data,fmt='%s') # adds the array to the bottom of the textfile
    np.savetxt(f,dem,fmt='%1.0f') # adds the array to the bottom of the textfile
    f.close()
with open(os.path.join(path_to_files,'KRH_Slope_' + str(year) + '.txt'),'a') as f:
    f.seek(0)
    np.savetxt(f,data,fmt='%s') # adds the array to the bottom of the textfile
    np.savetxt(f,slope_array,fmt='%1.0f') # adds the array to the bottom of the textfile
    f.close()
with open(os.path.join(path_to_files,'KRH_Aspect_' + str(year) + '.txt'),'a') as f:
    f.seek(0)
    np.savetxt(f,data,fmt='%s') # adds the array to the bottom of the textfile
    np.savetxt(f,aspect_array,fmt='%1.0f') # adds the array to the bottom of the textfile
    f.close()
   
# plotting the aspect and slope arrays:
plt.figure(figsize=(16,5))
plt.subplot(1,2,1)
plt.contourf(aspect_array,cmap='twilight',levels=np.linspace(0,360,37))
legend = plt.colorbar()
legend.ax.set_ylabel('Aspect ($\circ$)', rotation=270,fontsize=14,labelpad=20)
plt.axis('equal') 
plt.title('Aspect (' + str(year) + ')',fontsize=14)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
 
plt.subplot(1,2,2)
plt.contourf(slope_array,cmap='YlGnBu')
legend = plt.colorbar()
legend.ax.set_ylabel('Slope ($\circ$)', rotation=270,fontsize=14,labelpad=20)
plt.axis('equal') 
plt.title('Slope (' + str(year) + ')',fontsize=14)
plt.xlabel('Easting',fontsize=14)
plt.ylabel('Northing',fontsize=14)
 
plt.tight_layout()
plt.close()


