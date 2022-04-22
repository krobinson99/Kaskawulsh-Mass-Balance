# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:33:44 2020

@author: Juhyo
"""
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from netCDF4 import Dataset
from Model_functions_ver4 import regridXY_something
from Model_functions_ver4 import Check_out_that_MB


###get nan locs just in case
File_temp_in = 'Tempkaskonly2007.nc'
inT = Dataset(File_temp_in, "r")
T_var = 'Temperature'
T_array = inT.variables[T_var][:]
nanlocs = np.isnan(T_array[0,:,:])

###and get Zgrid
glacier = np.genfromtxt('kaskonly.txt', skip_header=1, delimiter=',')
Ix = glacier[:,2]
Iy = glacier[:,3]
Ih = glacier[:,1]
Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)


#script to load yearly net mass balance for each model run with final values

####1 plot asnow/aice/MF
#load params
params25 = np.loadtxt('final_params_katietest.csv')
aice_pp = params25[0,:]
asnow_pp = params25[1,:]
MF_pp = params25[2,:]
#print(params)

fig = plt.figure(1,figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(aice_p, asnow_p, MF_p)
plt.xlabel('aice_p')
plt.ylabel('asnow_p')
#ax.set_zlabel('MF_p')
#plt.zlabel('MF_p')
plt.xticks(np.arange(min(aice_p), max(aice_p),0.0000004))
plt.yticks(np.arange(min(asnow_p), max(asnow_p),0.0000004))
#ax.plot(aice_p, asnow_p, MF_p)

params10 = np.loadtxt('final_params_katietest2.csv')
aice_p = params10[0,:]
asnow_p = params10[1,:]
MF_p = params10[2,:]
#print(params)

fig = plt.figure(2,figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(aice_p, asnow_p, MF_p)
plt.xlabel('aice_p')
plt.ylabel('asnow_p')
plt.xticks(np.arange(min(aice_p), max(aice_p),0.0000004))
plt.yticks(np.arange(min(asnow_p), max(asnow_p),0.0000004))
#plt.zlabel('MF_p')
#ax.plot(aice_p, asnow_p, MF_p)

# 2d projections to compare param distributions
plt.figure(3,figsize=(14,14))
plt.subplot(3,2,1)
plt.title('aice vs asnow - 10 params')
plt.xlabel('a_ice')
plt.ylabel('a_snow')
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=45)
plt.plot(aice_p,asnow_p,'.')
plt.subplot(3,2,2)
plt.title('aice vs asnow - 25 params')
plt.xlabel('a_ice')
plt.ylabel('a_snow')
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=45)
plt.plot(aice_pp,asnow_pp,'.')
plt.subplot(3,2,3)
plt.title('aice vs MF - 10 params')
plt.xlabel('a_ice')
plt.ylabel('MF')
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=45)
plt.plot(aice_p,MF_p,'.')
plt.subplot(3,2,4)
plt.title('aice vs MF - 25 params')
plt.xlabel('a_ice')
plt.ylabel('MF')
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=45)
plt.plot(aice_pp,MF_pp,'.')
plt.subplot(3,2,5)
plt.title('asnow vs MF - 10 params')
plt.xlabel('a_snow')
plt.ylabel('MF')
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=45)
plt.plot(asnow_p,MF_p,'.')
plt.subplot(3,2,6)
plt.title('asnow vs MF - 25 params')
plt.xlabel('a_snow')
plt.ylabel('MF')
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=45)
plt.plot(asnow_pp,MF_pp,'.')
plt.tight_layout()

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

x =aice_p
y =asnow_p
z =MF_p

ax.scatter(x, y, z, c='r', marker='o')
plt.xticks(np.arange(min(aice_p), max(aice_p),0.0000004))
plt.yticks(np.arange(min(asnow_p), max(asnow_p),0.0000004))
plt.title('10 debris-free param combinations')
#plt.title('25 debris-free param combinations')
ax.set_xlabel('a ice')
ax.set_ylabel('a snow')
ax.set_zlabel('MF')

plt.show()


###open and concatenate model runs  

MB_prefix = 'Mb'
MB_sufix = '.nc'

for i in range(0, len(aice_p)):

    earlist = [2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017, \
               2018]
    
    #create empty array with shape [x,y,z] = 12 years, MB_array dims
    N_MB_container = np.empty((12,218,328))
    
    for ear in earlist:
        
        File_MB_in = MB_prefix + str(ear) + str(i) + MB_sufix

        MB_ear = Dataset(File_MB_in, "r")
        MB_var = 'MB'
        MB_array = MB_ear.variables[MB_var][:]
        
        N_MB = np.nansum(MB_array, axis = 0) #nansum returns the sum of array elements over a given axis treating Not a Numbers (NaNs) as zero.
        
        #nanlocs = np.where(np.isnan(MB_array[0,:,:]) == 1)
        #N_MB[nanlocs] = np.nan
        
        N_MB_container[np.where(np.asarray(earlist) == ear),:,:] = N_MB
    
    out_file_MB = 'N_MB' + str(i) + '.npy'
    np.save(out_file_MB, N_MB_container)
    

###open and extract run NMB values and create map of mean NMB + std
overall_MB = np.empty((len(aice_p)-1,218,328))
for i in range(0, len(aice_p)):

    
    
    out_file_MB = 'N_MB' + str(i) + '.npy'

    run_MB = np.load(out_file_MB)
    
    run_mean_MB = np.mean(run_MB, axis = 0)
    run_mean_MB[nanlocs] = np.nan
    
    overall_MB[i-1,:,:] = run_mean_MB

#save output of run means
out_run_means = 'runMBs.npy'
np.save(out_run_means, overall_MB)

###Create maps
plt.contourf(np.flipud(np.mean(overall_MB, axis = 0)), cmap = 'RdYlBu', levels = np.linspace(-12,7, 20))
plt.colorbar()
plt.xlabel('test')
plt.contour(np.flipud(np.mean(overall_MB, axis = 0)), colors = 'black', levels = 0)
plt.title('Mean NMB')
plt.open()

plt.contourf(np.flipud(np.std(overall_MB, axis = 0)), cmap = 'cividis', levels = np.linspace(0, 1.5, 7))
plt.colorbar()
#plt.contour(np.flipud(np.mean(overall_MB, axis = 0)), colors = 'black', levels = 0)
plt.title('STD NMB')

p_overlap = (np.flipud(np.std(overall_MB, axis = 0))/np.flipud(np.mean(overall_MB, axis = 0)))
plt.contourf(p_overlap, cmap = 'rainbow', levels = np.linspace(-1, 1, 10))
plt.colorbar()
#plt.contour(np.flipud(np.mean(overall_MB, axis = 0)), colors = 'black', levels = 0)
plt.title('percent std/mean NMB')

###extract bf
master_mb = np.mean(overall_MB, axis = 0)
master_std = np.std(overall_MB, axis = 0)

kask = np.loadtxt('kaskonly_bf.txt', skiprows = 1, delimiter = ',')
BF_ids = kask[:,5]
Ix = kask[:,3]
Iy = kask[:,4]
bf_ids, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, BF_ids)

def balanceflux(upstream):
    return ((upstream*(200.**2.))/1000000000.)*1./0.9

bf_locs = []
for i in range(0, 10):
    bf_locs.append(np.where(bf_ids == i))

Balance_fluxes = []
BF_NA = balanceflux(np.nansum(master_mb[bf_locs[1]]))
Balance_fluxes.append(BF_NA)
BF_CA = balanceflux(np.nansum(master_mb[bf_locs[2]]))
Balance_fluxes.append(BF_CA)
BF_SW = balanceflux(np.nansum(master_mb[bf_locs[3]]))
Balance_fluxes.append(BF_SW)
BF_SA = balanceflux(np.nansum(master_mb[bf_locs[4]]))
Balance_fluxes.append(BF_SA)

BF_KW5 = balanceflux(np.nansum(master_mb[bf_locs[1]]) + np.nansum(master_mb[bf_locs[2]]) + np.nansum(master_mb[bf_locs[5]]))
Balance_fluxes.append(BF_KW5)
BF_KW4 = balanceflux(np.nansum(master_mb[bf_locs[1]]) + np.nansum(master_mb[bf_locs[2]]) + np.nansum(master_mb[bf_locs[5]]) + np.nansum(master_mb[bf_locs[6]]))
Balance_fluxes.append(BF_KW4)
BF_KW3 = balanceflux(np.nansum(master_mb[bf_locs[1]]) + np.nansum(master_mb[bf_locs[2]]) + np.nansum(master_mb[bf_locs[5]]) + np.nansum(master_mb[bf_locs[6]]) + np.nansum(master_mb[bf_locs[7]]) + np.nansum(master_mb[bf_locs[3]]))
Balance_fluxes.append(BF_KW3)
BF_KW2 = balanceflux(np.nansum(master_mb[bf_locs[1]]) + np.nansum(master_mb[bf_locs[2]]) + np.nansum(master_mb[bf_locs[5]]) + np.nansum(master_mb[bf_locs[6]]) + np.nansum(master_mb[bf_locs[7]]) + np.nansum(master_mb[bf_locs[3]]) + np.nansum(master_mb[bf_locs[8]]))
Balance_fluxes.append(BF_KW2)
BF_KW1 = balanceflux(np.nansum(master_mb[bf_locs[1]]) + np.nansum(master_mb[bf_locs[2]]) + np.nansum(master_mb[bf_locs[3]]) + np.nansum(master_mb[bf_locs[4]]) + np.nansum(master_mb[bf_locs[5]]) + np.nansum(master_mb[bf_locs[6]]) + np.nansum(master_mb[bf_locs[7]]) + np.nansum(master_mb[bf_locs[8]]) + np.nansum(master_mb[bf_locs[9]]))
Balance_fluxes.append(BF_KW1)


Balance_stds = []
BF_NA = balanceflux(np.nansum(master_std[bf_locs[1]]))
Balance_stds.append(BF_NA)
BF_CA = balanceflux(np.nansum(master_std[bf_locs[2]]))
Balance_stds.append(BF_CA)
BF_SW = balanceflux(np.nansum(master_std[bf_locs[3]]))
Balance_stds.append(BF_SW)
BF_SA = balanceflux(np.nansum(master_std[bf_locs[4]]))
Balance_stds.append(BF_SA)

BF_KW5 = balanceflux(np.nansum(master_std[bf_locs[1]]) + np.nansum(master_std[bf_locs[2]]) + np.nansum(master_std[bf_locs[5]]))
Balance_stds.append(BF_KW5)
BF_KW4 = balanceflux(np.nansum(master_std[bf_locs[1]]) + np.nansum(master_std[bf_locs[2]]) + np.nansum(master_std[bf_locs[5]]) + np.nansum(master_std[bf_locs[6]]))
Balance_stds.append(BF_KW4)
BF_KW3 = balanceflux(np.nansum(master_std[bf_locs[1]]) + np.nansum(master_std[bf_locs[2]]) + np.nansum(master_std[bf_locs[5]]) + np.nansum(master_std[bf_locs[6]]) + np.nansum(master_std[bf_locs[7]]) + np.nansum(master_std[bf_locs[3]]))
Balance_stds.append(BF_KW3)
BF_KW2 = balanceflux(np.nansum(master_std[bf_locs[1]]) + np.nansum(master_std[bf_locs[2]]) + np.nansum(master_std[bf_locs[5]]) + np.nansum(master_std[bf_locs[6]]) + np.nansum(master_std[bf_locs[7]]) + np.nansum(master_std[bf_locs[3]]) + np.nansum(master_std[bf_locs[8]]))
Balance_stds.append(BF_KW2)
BF_KW1 = balanceflux(np.nansum(master_std[bf_locs[1]]) + np.nansum(master_std[bf_locs[2]]) + np.nansum(master_std[bf_locs[3]]) + np.nansum(master_std[bf_locs[4]]) + np.nansum(master_std[bf_locs[5]]) + np.nansum(master_std[bf_locs[6]]) + np.nansum(master_std[bf_locs[7]]) + np.nansum(master_std[bf_locs[8]]) + np.nansum(master_std[bf_locs[9]]))
Balance_stds.append(BF_KW1)



B_M = []
B_S = []
for i in overall_MB:
    
    Balance_fluxes = []
    
    #BF_all = balanceflux(np.nansum(i))
    #Balance_fluxes.append(BF_all)
    
    BF_NA = balanceflux(np.nansum(i[bf_locs[1]]))
    Balance_fluxes.append(BF_NA)
    BF_CA = balanceflux(np.nansum(i[bf_locs[2]]))
    Balance_fluxes.append(BF_CA)
    BF_SW = balanceflux(np.nansum(i[bf_locs[3]]))
    Balance_fluxes.append(BF_SW)
    BF_SA = balanceflux(np.nansum(i[bf_locs[4]]))
    Balance_fluxes.append(BF_SA)
    
    BF_KW5 = balanceflux(np.nansum(i[bf_locs[1]]) + np.nansum(i[bf_locs[2]]) + np.nansum(i[bf_locs[5]]))
    Balance_fluxes.append(BF_KW5)
    BF_KW4 = balanceflux(np.nansum(i[bf_locs[1]]) + np.nansum(i[bf_locs[2]]) + np.nansum(i[bf_locs[5]]) + np.nansum(i[bf_locs[6]]))
    Balance_fluxes.append(BF_KW4)
    BF_KW3 = balanceflux(np.nansum(i[bf_locs[1]]) + np.nansum(i[bf_locs[2]]) + np.nansum(i[bf_locs[5]]) + np.nansum(i[bf_locs[6]]) + np.nansum(i[bf_locs[7]]) + np.nansum(i[bf_locs[3]]))
    Balance_fluxes.append(BF_KW3)
    BF_KW2 = balanceflux(np.nansum(i[bf_locs[1]]) + np.nansum(i[bf_locs[2]]) + np.nansum(i[bf_locs[5]]) + np.nansum(i[bf_locs[6]]) + np.nansum(i[bf_locs[7]]) + np.nansum(i[bf_locs[3]]) + np.nansum(i[bf_locs[8]]))
    Balance_fluxes.append(BF_KW2)
    BF_KW1 = balanceflux(np.nansum(i[bf_locs[1]]) + np.nansum(i[bf_locs[2]]) + np.nansum(i[bf_locs[3]]) + np.nansum(i[bf_locs[4]]) + np.nansum(i[bf_locs[5]]) + np.nansum(i[bf_locs[6]]) + np.nansum(i[bf_locs[7]]) + np.nansum(i[bf_locs[8]]) + np.nansum(i[bf_locs[9]]))
    Balance_fluxes.append(BF_KW1)
    
    B_M.append(Balance_fluxes)
    
STD_runs = np.std(B_M, axis = 0).tolist()
BFm_run = np.mean(B_M, axis = 0).tolist()
BFmed_run = np.median(B_M, axis = 0).tolist()

#line of code to get gate-by-gate STD
#one gate
std_l = []
for j in range(0,10):
    st = np.std(np.nansum(overall_MB[:,bf_locs[j][0],bf_locs[j][1]], axis = 1))
    std_l.append(st)
#multiple gate
np.std(np.nansum(overall_MB[:,bf_locs[1][0],bf_locs[1][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[2][0],bf_locs[2][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[5][0],bf_locs[5][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[6][0],bf_locs[6][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[7][0],bf_locs[7][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[3][0],bf_locs[3][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[8][0],bf_locs[8][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[4][0],bf_locs[4][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[9][0],bf_locs[9][1]], axis = 1) + \
       np.nansum(overall_MB[:,bf_locs[0][0],bf_locs[0][1]], axis = 1))

###extract annual MB differences

def mbs(upstream):
    return ((upstream*(200.**2.))/1095554179.)

NMB_list = []
for i in overall_MB:
    yearlyMB = mbs(np.nansum(i))
    NMB_list.append(yearlyMB)
    
mb_mean_year = []

###################################get each flux for each datasource#######################################
###get upstream MB a test

def mbs(upstream, area):
    #return (((upstream*(200.**2.))*1000000000)/area)
    return (upstream*1000000000.)/area
areas = [218454179., 318854179., 107254179., 262254179., 569654179., 591254179., 720354179., 749854179., 1048554179.]
fluxes = BFm_run

upstream_af = []
for j in range(0, len(areas)):
    upstream_af.append(mbs(fluxes[j], areas[j]))
   

### now perform for actual fluxes
def mbs(upstream, area):
    #return (((upstream*(200.**2.))*1000000000)/area)
    return (upstream*1000000000.)/area

areas = [218454179., 318854179., 107254179., 262254179., 569654179., 591254179., 720354179., 749854179., 1048554179.]
actual_fluxes = [0.221, 0.209, 0.042, 0.075, 0.341, 0.326, 0.326, 0.174, 0.195]

upstream_af = []
for j in range(0, len(areas)):
    upstream_af.append(mbs(actual_fluxes[j], areas[j]))
    
### now perform for actual fluxes
def mbs(upstream, area):
    #return (((upstream*(200.**2.))*1000000000)/area)
    
    return (upstream*1000000000.)/area
### now perform for dhdt
    #get fluxes
def balanceflux(upstream):
    return ((upstream*(20.**2.))/1000000000.)*1./0.9

areas = [218454179., 318854179., 107254179., 262254179., 569654179., 591254179., 720354179., 749854179., 1048554179.]
dhdt_indivsums = [-112776., -262124., -100681., -308093., -76639., -68087.,-77829., -101492., -159519.]
dhdt_cumsums = [-112776., -262124., -100681., -308093., -451539., -519626., -698138., -799628., -1267240.]

dhdt_bf = []
for j in range(0, len(areas)):
    dhdt_bf.append(balanceflux(dhdt_cumsums[j]))
    
    #get MB
def mbs(upstream, area):
    #return (((upstream*(200.**2.))*1000000000)/area)
    return (upstream*1000000000.)/area

dhdt_mb = []
for j in range(0, len(areas)):
    dhdt_mb.append(mbs(np.asarray(dhdt_bf)[j], areas[j]))
    
    ##############################################################################################################
    
    
   
###final tuning check to match geodetic MB slope
for j in overall_MB:
    flat_netMB = np.ndarray.flatten(j)
    flat_z = np.ndarray.flatten(Zgrid)

    MBmeans, MBstd, MBmedians = Check_out_that_MB(flat_netMB, flat_z, 2100.)
    
    mb_mean_year.append(MBmeans)
    
dhdt = np.genfromtxt('dhdt_table.txt', skip_header = 1, delimiter = ',')

dhdt_vals = dhdt[:,2]
z = dhdt[:,3]
MBmeans_dhdt, MBstd_dhdt, MBmedians_dhdt = Check_out_that_MB(flat_netMB, flat_z, 2100.)

zbins = np.linspace(0, 4000, 41)
for p in mb_mean_year:
    plt.plot(p, zbins)
    
plt.plot(MBmeans_dhdt, zbins, linewidth = 4, c ='k')
plt.ylabel('Elevation bin')
plt.xlabel('NMB')


np.nanmean()






    

    
    

