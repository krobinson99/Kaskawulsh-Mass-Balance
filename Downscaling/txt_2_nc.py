####Tool for converting large csv textfiles to netcdf####

import matplotlib.dates as date
import netCDF4
from netCDF4 import Dataset
import numpy as np
import datetime as dt
from Model_functions_ver4 import regridXY_something


#name = 'GL_mid'  
#name_list = ['Meltmid_dd', 'Meltmid_eddarc', 'Meltmid_eddre', 'MBGmid_dd', 'MBGmid_eddarc', 'MBGmid_eddre', 'Snowmid_dd', 'Snowmid_edd', 'Sradmid'] 

name = 'kaskonly'
#name_list =  [ 'Meltkaskice2017_dd']

#name = 'WS' 
#longform for all files 
#name_list = ['Meltkaskice_dd', 'Meltkaskice_eddarc', 'MBGkaskice_dd', 'MBGkaskice_eddarc', 'Tempkaskice', 'Sradkaskice', 'Snowkaskice_dd', 'Snowkaskice_edd', 'precipaskice', 'netSnowkaskice']
#name_list = ['Meltkasktopo_dd', 'Meltkasktopo_eddarc', 'MBGkasktopo_dd', 'MBGkasktopo_eddarc', 'Tempkasktopo', 'Sradkasktopo', 'Snowkasktopo_dd', 'Snowkasktopo_edd', 'precipkasktopo', 'netSnowkasktopo']
#name_list = ['netSnowkask2010','netSnowkask2011','netSnowkask2012','netSnowkask2013','netSnowkaskonly2014','netSnowkaskonly2015', 'netSnowkaskonly2016','netSnowkaskonly2017','netSnowkaskonly2018']
#name_list = ['Sradkask2006']
#shortform for ETIM outpputs
name_list = ['Tempkaskonly2006']

#WS outputs
#name_list = ['TempWS', 'SradWS', 'precipWS']

start_year = 2006

for n in name_list:
    print(n)

    File_glacier_in = name + '.txt'
    File_in = n + '.txt'
    File_out = n + '.nc'
    File_grid_out = name + '_grid.txt'
    file_in = 'ref' + str(start_year) + '.nc'
    print('using ' + file_in)
    
    fh = Dataset(file_in, "r")
    
    
    ###Extract grid geometry###
    
    glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
    
    IH = glacier[:,6]
    inval_loc = np.where(IH == -9999)
    IH[inval_loc] = IH[inval_loc] * np.nan
    Ih = IH
#Ih = np.delete(IH, inval_loc)
    
    IX = glacier[:,4]
    Ix = IX
#    Ix = np.delete(IX, inval_loc)
    IY = glacier[:,5]
    Iy = IY
#    Iy = np.delete(IY, inval_loc)
    
    #note that Zgrid is also used as a nan mask
    Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
    
    np.savetxt(File_grid_out, Zgrid)
    
    ###Create empty .nc file for MB values###
    
    f = Dataset(File_out, 'w', format='NETCDF4') #'w' stands for write
    
    time = f.createDimension('time', None)
    #z = f.createDimension('z', len(Ih))
    y = f.createDimension('y', len(ybounds))
    x = f.createDimension('x', len(xbounds))
    #f.createDimension('date', None)
    
    #load subset into new netcdffile
    times = f.createVariable('time', np.float64, 'time')
    #Zs = f.createVariable('z', np.int32, 'z')
    Ys = f.createVariable('y', np.float32, 'y')
    Xs = f.createVariable('x', np.float32, 'x')  
    mb = f.createVariable('Precipitation', np.float32, ('time', 'y', 'x'))
    
    # Global Attributes 
    f.description = 'Model input total precipitation'   
    f.source = 'Downscaled NARR' 
    # Variable Attributes  
    Ys.units = 'm'  
    Xs.units = 'm'  
    #Zs.units = 'm' 
    mb.units = 'w m2' 
    times.units = 'hours since 1800-01-01 00:00:00' 
    times.calendar = 'gregorian' 
    times.delta_t = '0000-00-00 03:00:00'
    #time.actual_range = [1835688. , 1836405]
    
    #insert metadata into .nc file 
    XX = xbounds
    YY = ybounds
    #Z = Ih
    time_nc = fh.variables['time'][:]
    b = fh.variables['time']
    dtime = netCDF4.num2date(b[:], b.units)
    
    Xs[:] = XX #The "[:]" at the end of the variable instance is necessary
    Ys[:] = YY
    #Zs[:] = Z
    times[:] = netCDF4.date2num(dtime, units = times.units, calendar = times.calendar)
    
    ###Get values to input###
    t = 0
    with open(File_in) as infile:
        
        for line in infile:
            
            if "[" or "]" in line:
                line_clini = line.replace("[", "")
                line_cl = line_clini.replace("]", "")
            else:
                line_cl = line
            
            text = line_cl[1:-2].split(',')
            MaBa = np.loadtxt(text, delimiter=',')
            
            maba, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, MaBa)
    
            #insert values into netcdf file        
            mb[t,:,:] = maba
            t = t + 1 
            print("\r"+ "Day:" + str(t/8))
            
    print("Finished:"+  n)
    
            
    infile.close()
    f.close()
    fh.close()
    
    start_year += 1
        

