####Tool for converting large csv textfiles to netcdf####

#import matplotlib.dates as date
import netCDF4
from netCDF4 import Dataset
import numpy as np
#import datetime as dt
import os
from Model_functions_ver4 import regridXY_something

name = 'kaskonly'

path = 'D:/Downscaled_files/includes_CAtrib/2006-2022' 
variable = 'netSnow'
years = [2006]
file_list = []
for y in years:
    fname = str(variable) + str(name) + str(y)
    #file_list.append(os.path.join(path,fname))
    file_list.append(fname)

ref_path = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Ref_files'

start_year = years[0]

for n in file_list:
    print(n)

    File_glacier_in = 'kaskonly_deb.txt'
    File_in = os.path.join(path, str(n) + '.txt')
    File_out = os.path.join(path, str(n) + '.nc')
    File_grid_out = os.path.join(path, str(name) + '_grid.txt')
    file_in = os.path.join(ref_path, 'ref' + str(start_year) + '.nc')
    
    fh = Dataset(file_in, "r")
    
    
    ###Extract grid geometry###
    
    glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
    Ih = glacier[:,2]
    inval_loc = np.where(Ih == -9999)
    Ih[inval_loc] = Ih[inval_loc] * np.nan
    Ix = glacier[:,3]
    Iy = glacier[:,4]

    
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
        

