# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 09:54:13 2023

Script to calculate the linear regression coefficients used in the precipitation downscaling.

@author: katierobinson
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 10:08:16 2023

Calculate how many times each NARR gridcell is used in the precip downscaling

@author: katierobinson
"""
import numpy as np
from pyproj import Proj
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from Model_functions_ver3 import rainy_day_funk
from Model_functions_ver3 import regridXY_something
#from Model_functions_ver3 import closest_node
import sys
import os 

NARR_subregions_original = [8,9,10,14,15,16,20,21,22,26,27,28]
NARR_subregions_full = np.arange(0,36)
NARR_subregions_inner = [7,8,9,10,13,14,15,16,19,20,21,22,25,26,27,28]
NARR_subregions_kw = [14,15,19,20,21]
NARR_subregions_orographic = [9,10,14,15,16,19,20,21,22,25,26,27,28]

#NARR_subregions = NARR_subregions_original

glacier_outline = 'kask_catchment'
catchment = True #are we looking at the kaskawulsh only (False) or the whole catchment? (True)
static_surface = True
dynamic_surface = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Downscaling/Ih_1977.txt'

UTM = '7N'

rawnarr_inputs = 'F:/Mass Balance Model/CoarseNARR_KW'

# Load coarse NARR grid
File_CoarseElev_in = 'kaskCE.nc'
inCE = Dataset(File_CoarseElev_in, "r")
CE_array = inCE.variables['hgt']
elev = CE_array[0,:,:] # surface geopotential height, units = m, shape = (6,6)

year = 2020
File_elev_in = os.path.join(rawnarr_inputs,'kaskhgt.' + str(year) + '.nc')

inE = Dataset(File_elev_in, "r")
E_var = 'hgt'
sys.stdout.flush()
E_array = inE.variables[E_var][:]
#Get UTM values for each reanalysis grid point to cycle between nodes on larger 
#modelled glaciers
lons = inE.variables['lon'][:]
lats = inE.variables['lat'][:]

myProj = Proj('+proj=utm +zone=' + UTM + ', +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
UTMx, UTMy = myProj(lons, lats)  # converts lat/lons of coarse NARR grid to easting, northing on WGS84 projection.

#create list of array positions ungridded
UTMx_list = UTMx.ravel()
UTMy_list = UTMy.ravel()
grid_pts = np.stack((UTMy_list, UTMx_list))

File_glacier_in = glacier_outline + '.txt' 
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')

#Extract grid geometry, correcting for dem/shapefile non-overlaps
if static_surface == True:
    if catchment == False:
        IH = glacier[:,2]
        inval_loc = np.where(IH == -9999)
        Ih =np.delete(IH, inval_loc)
        Ix = glacier[:,3]
        Iy = glacier[:,4]
    else:
        IH = glacier[:,6]
        inval_loc = np.where(IH == -9999)
        Ih = np.delete(IH, inval_loc)
        Ix = glacier[:,4]
        Iy = glacier[:,5]
        sfc_type = glacier[:,8]
    print('Static Surface Loaded')
else:
    Ih = np.loadtxt(dynamic_surface)
    Ix = glacier[:,4]
    Iy = glacier[:,5]
    print('Dynamic Surface Loaded')

Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
Ygrid = np.flipud(Ygrid)
Sfc_grid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, sfc_type)

nanlocs = np.where(np.isnan(Zgrid))
offglacier = np.where(Sfc_grid==1)
onglacier = np.where(Sfc_grid==0)


print('Domain loaded')

def get_downscaling_coeffs(NARR_subregions,path,years=np.arange(1979,2023),normalize_coeffs=False):
    '''
    
    '''
    r2 = []
    b00 = []
    b1 = []
    b2 = []
    b3 = []
    b4 = []
    b5 = []
    b6 = []
    for year in years:
        print(year)
         
        File_precip_in = os.path.join(rawnarr_inputs,'kaskapcp.' + str(year) + '.nc')
        
        inP = Dataset(File_precip_in, "r")
        P_var = 'apcp'
        sys.stdout.flush()
        P_array = inP.variables[P_var][:]
        
        coeffs = np.zeros((len(P_array),9))
        for i in range(0,len(P_array)):
            #print('Day =',i+1)
            dailyP = P_array[i]/1000 # for units of m w.e.
            if normalize_coeffs == False:
                r_beta2, b_coeffs, b0 = rainy_day_funk(elev.ravel()[NARR_subregions], dailyP.ravel()[NARR_subregions], UTMy_list[NARR_subregions], UTMx_list[NARR_subregions])
            else:
                r_beta2, b_coeffs, b0 = rainy_day_funk(elev.ravel()[NARR_subregions]/np.max(elev.ravel()[NARR_subregions]), dailyP.ravel()[NARR_subregions], UTMy_list[NARR_subregions]/np.max(UTMy_list[NARR_subregions]), UTMx_list[NARR_subregions]/np.max(UTMx_list[NARR_subregions]))                
            
            r2.append(r_beta2)
            b00.append(b0)
            b1.append(b_coeffs[0])
            b2.append(b_coeffs[1])
            b3.append(b_coeffs[2])
            b4.append(b_coeffs[3])
            b5.append(b_coeffs[4])
            b6.append(b_coeffs[5])
            
            coeffs[i,0] = i+1
            coeffs[i,1] = r_beta2
            coeffs[i,2] = b0
            coeffs[i,3] = b_coeffs[0] #b1
            coeffs[i,4] = b_coeffs[1] #b2
            coeffs[i,5] = b_coeffs[2] #b3
            coeffs[i,6] = b_coeffs[3] #b4
            coeffs[i,7] = b_coeffs[4] #b5
            coeffs[i,8] = b_coeffs[5] #b6
        
        #file = os.path.join(path,str(year)+'_downscalingcoeffs.npy')
        #np.save(file,coeffs)
    
    print('done!')
    return r2, b00, b1, b2, b3, b4, b5, b6
    
original = get_downscaling_coeffs(NARR_subregions_original,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_original')
full =  get_downscaling_coeffs(NARR_subregions_full,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_full')  
inner = get_downscaling_coeffs(NARR_subregions_inner,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_inner')
kw = get_downscaling_coeffs(NARR_subregions_kw,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_kw')
orographic = get_downscaling_coeffs(NARR_subregions_orographic,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_orographic')

original_norm = get_downscaling_coeffs(NARR_subregions_original,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_original',normalize_coeffs=True)
full_norm =  get_downscaling_coeffs(NARR_subregions_full,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_full',normalize_coeffs=True)  
inner_norm = get_downscaling_coeffs(NARR_subregions_inner,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_inner',normalize_coeffs=True)
kw_norm = get_downscaling_coeffs(NARR_subregions_kw,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_kw',normalize_coeffs=True)
orographic_norm = get_downscaling_coeffs(NARR_subregions_orographic,'D:/Downscaled_files/Catchment/NARRsubregions_tests/subregions_orographic',normalize_coeffs=True)



meanlineprops = dict(linestyle='--', linewidth=2.5, color='purple')
meanpointprops = dict(marker='D',markeredgecolor='purple',markerfacecolor='purple')

plt.figure(figsize=(16,6))
plt.suptitle('For NARR Daily Precipitation (1979--2022)',fontsize=14,y=1.01)
labels=['Original','Full','Inner','KW','Orographic']
# 8 panels (1 per coefficient, w 5 boxes on each panel (one per subregions test))
plt.subplot(2,4,1)
plt.boxplot([original[0],full[0],inner[0],kw[0],orographic[0]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('r$^2_{\\beta}$',fontsize=14)
plt.subplot(2,4,2)
plt.boxplot([original[1],full[1],inner[1],kw[1],orographic[1]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_0$',fontsize=14)
plt.subplot(2,4,3)
plt.boxplot([original[2],full[2],inner[2],kw[2],orographic[2]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_1$',fontsize=14)
plt.subplot(2,4,4)
plt.boxplot([original[3],full[3],inner[3],kw[3],orographic[3]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_2$',fontsize=14)
plt.subplot(2,4,5)
plt.boxplot([original[4],full[4],inner[4],kw[4],orographic[4]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_3$',fontsize=14)
plt.subplot(2,4,6)
plt.boxplot([original[5],full[5],inner[5],kw[5],orographic[5]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_4$',fontsize=14)
plt.subplot(2,4,7)
plt.boxplot([original[6],full[6],inner[6],kw[6],orographic[6]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_5$',fontsize=14)
plt.subplot(2,4,8)
plt.boxplot([original[7],full[7],inner[7],kw[7],orographic[7]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_6$',fontsize=14)
#plt.legend()
plt.tight_layout()
#plt.savefig('Precipitation_downscaling_coefficients.png',bbox_inches='tight')

plt.figure(figsize=(16,6))
plt.suptitle('For NARR Daily Precipitation (1979--2022)',fontsize=14,y=1.01)
labels=['Original','Full','Inner','KW','Orographic']
# 8 panels (1 per coefficient, w 5 boxes on each panel (one per subregions test))
plt.subplot(2,4,1)
plt.boxplot([original[0],full[0],inner[0],kw[0],orographic[0]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('r$^2_{\\beta}$',fontsize=14)
plt.subplot(2,4,2)
plt.boxplot([original[1],full[1],inner[1],kw[1],orographic[1]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_0$',fontsize=14)
plt.subplot(2,4,3)
plt.boxplot([original[2],full[2],inner[2],kw[2],orographic[2]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_1$',fontsize=14)
plt.subplot(2,4,4)
plt.boxplot([original[3],full[3],inner[3],kw[3],orographic[3]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_2$',fontsize=14)
plt.subplot(2,4,5)
plt.boxplot([original[4],full[4],inner[4],kw[4],orographic[4]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_3$',fontsize=14)
plt.subplot(2,4,6)
plt.boxplot([original[5],full[5],inner[5],kw[5],orographic[5]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_4$',fontsize=14)
plt.subplot(2,4,7)
plt.boxplot([original[6],full[6],inner[6],kw[6],orographic[6]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_5$',fontsize=14)
plt.subplot(2,4,8)
plt.boxplot([original[7],full[7],inner[7],kw[7],orographic[7]],labels=labels,showfliers=True,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_6$',fontsize=14)
#plt.legend()
plt.tight_layout()
#plt.savefig('Precipitation_downscaling_coefficients_withfliers.pdf',bbox_inches='tight')

plt.figure(figsize=(16,6))
plt.suptitle('For NARR Daily Precipitation (1979--2022)',fontsize=14,y=1.01)
labels=['Original','Inner','Orographic','Full']
# 8 panels (1 per coefficient, w 5 boxes on each panel (one per subregions test))
plt.subplot(2,4,1)
plt.boxplot([original[0],inner[0],orographic[0],full[0]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('r$^2_{\\beta}$',fontsize=14)
plt.subplot(2,4,2)
plt.boxplot([original[1],inner[1],orographic[1],full[1]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_0$',fontsize=14)
plt.subplot(2,4,3)
plt.boxplot([original[2],inner[2],orographic[2],full[2]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_1$',fontsize=14)
plt.subplot(2,4,4)
plt.boxplot([original[3],inner[3],orographic[3],full[3]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_2$',fontsize=14)
plt.subplot(2,4,5)
plt.boxplot([original[4],inner[4],orographic[4],full[4]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_3$',fontsize=14)
plt.subplot(2,4,6)
plt.boxplot([original[5],inner[5],orographic[5],full[5]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_4$',fontsize=14)
plt.subplot(2,4,7)
plt.boxplot([original[6],inner[6],orographic[6],full[6]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_5$',fontsize=14)
plt.subplot(2,4,8)
plt.boxplot([original[7],inner[7],orographic[7],full[7]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_6$',fontsize=14)
#plt.legend()
plt.tight_layout()
#plt.savefig('Precipitation_downscaling_coefficients_noKW.png',bbox_inches='tight')

plt.figure(figsize=(16,6))
plt.suptitle('NARR Daily Precipitation (1979--2022) coeffs, normalized X,Y,Z',fontsize=14,y=1.01)
labels=['Original','Inner','Orographic','Full']
# 8 panels (1 per coefficient, w 5 boxes on each panel (one per subregions test))
plt.subplot(2,4,1)
plt.boxplot([original_norm[0],inner_norm[0],orographic_norm[0],full_norm[0]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('r$^2_{\\beta}$',fontsize=14)
plt.subplot(2,4,2)
plt.boxplot([original_norm[1],inner_norm[1],orographic_norm[1],full_norm[1]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_0$',fontsize=14)
plt.subplot(2,4,3)
plt.boxplot([original_norm[2],inner_norm[2],orographic_norm[2],full_norm[2]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_1$',fontsize=14)
plt.subplot(2,4,4)
plt.boxplot([original_norm[3],inner_norm[3],orographic_norm[3],full_norm[3]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_2$',fontsize=14)
plt.subplot(2,4,5)
plt.boxplot([original_norm[4],inner_norm[4],orographic_norm[4],full_norm[4]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_3$',fontsize=14)
plt.subplot(2,4,6)
plt.boxplot([original_norm[5],inner_norm[5],orographic_norm[5],full_norm[5]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_4$',fontsize=14)
plt.subplot(2,4,7)
plt.boxplot([original_norm[6],inner_norm[6],orographic_norm[6],full_norm[6]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_5$',fontsize=14)
plt.subplot(2,4,8)
plt.boxplot([original_norm[7],inner_norm[7],orographic_norm[7],full_norm[7]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_6$',fontsize=14)
#plt.legend()
plt.tight_layout()
#plt.savefig('Precipitation_downscaling_coefficients_normalizedXYZ.png',bbox_inches='tight')


# plot x,y,z vs P

def plot_P_regression(downscaling_coeffs,ymin,ymax):
    r2 = np.mean(downscaling_coeffs[0])
    b0 = np.mean(downscaling_coeffs[1])
    b1 = np.mean(downscaling_coeffs[2])
    b2 = np.mean(downscaling_coeffs[3])
    b3 = np.mean(downscaling_coeffs[4])
    b4 = np.mean(downscaling_coeffs[5])
    b5 = np.mean(downscaling_coeffs[6])
    b6 = np.mean(downscaling_coeffs[7])
    
    P = []
    for z in range(0,len(Ih)):
        Plocal = (b0 + (b1*Ix[z]) + (b2*Iy[z]) + (b3*Ix[z]*Iy[z]) + (b4*(Ix[z]**2)) + (b5*(Iy[z]**2)) + (b6*Ih[z]))
        P.append(Plocal)
        
    plt.subplot(1,3,1)
    plt.scatter(Ix,P,alpha=0.2,s=5)
    plt.xlabel('Easting (m)')
    plt.ylabel('Downscaled Precip (m w.e.)')
    plt.ylim(ymin,ymax)
    plt.subplot(1,3,2)
    plt.scatter(Iy,P,alpha=0.2,s=5)
    plt.xlabel('Northing (m)')
    plt.ylabel('Downscaled Precip (m w.e.)')
    plt.ylim(ymin,ymax)
    plt.subplot(1,3,3)
    plt.scatter(Ih,P,alpha=0.2,s=5)
    plt.xlabel('Elevation (m a.s.l.)')
    plt.ylabel('Downscaled Precip (m w.e.)')
    plt.ylim(ymin,ymax)
    plt.tight_layout()
    #return P
  
plt.figure(figsize=(10,3))
plt.suptitle('Young et al. (2021) NARR nodes',y=1.01)
plot_P_regression(original,0.001,0.00225)
#plt.savefig('originalnodes_regression.png',bbox_inches='tight')

plt.figure(figsize=(10,3))
plt.suptitle('All NARR nodes',y=1.01)
plot_P_regression(full,-0.001,0.0015)
#plt.savefig('fullnodes_regression.png',bbox_inches='tight')

plt.figure(figsize=(10,3))
plt.suptitle('Inner NARR nodes',y=1.01)
plot_P_regression(inner,0.001,0.00225)
#plt.savefig('innernodes_regression.png',bbox_inches='tight')

plt.figure(figsize=(10,3))
plt.suptitle('Kaskawulsh NARR nodes',y=1.01)
plot_P_regression(kw,-0.001,0.001)
#plt.savefig('kwnodes_regression.png',bbox_inches='tight')

plt.figure(figsize=(10,3))
plt.suptitle('Orographic NARR nodes',y=1.01)
plot_P_regression(orographic,0.001,0.00225)
#plt.savefig('orographicnodes_regression.png',bbox_inches='tight')

# plot each b value multiplied by the average value of the regressor 
# (e.g., x, y, xx, yy, xy, z) so we could see the approximate magnitude 
# of the terms in the first equation for P

# get avg val of each regressor
X = np.nanmean(Xgrid)
Y = np.nanmean(Ygrid)
XY = X*Y
XX = X**2
YY = Y**2
Z = np.nanmean(Zgrid)

plt.figure(figsize=(16,6))
plt.suptitle('Precipitation Regression Coefficients (Regular X,Y,Z)',fontsize=14,y=1.02)
labels=['Original','Inner','Orographic','Full']
# 8 panels (1 per coefficient, w 5 boxes on each panel (one per subregions test))
plt.subplot(2,4,1)
plt.boxplot([original[0],inner[0],orographic[0],full[0]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('r$^2_{\\beta}$',fontsize=14)
plt.ylim(0,1.1)
plt.subplot(2,4,2)
plt.boxplot([original[1],inner[1],orographic[1],full[1]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_0$',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,3)
plt.boxplot([np.array(original[2])*X,np.array(inner[2])*X,np.array(orographic[2])*X,np.array(full[2])*X],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_1$X',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,4)
plt.boxplot([np.array(original[3])*Y,np.array(inner[3])*Y,np.array(orographic[3])*Y,np.array(full[3])*Y],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_2$Y',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,5)
plt.boxplot([np.array(original[4])*XY,np.array(inner[4])*XY,np.array(orographic[4])*XY,np.array(full[4])*XY],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_3$XY',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,6)
plt.boxplot([np.array(original[5])*XX,np.array(inner[5])*XX,np.array(orographic[5])*XX,np.array(full[5])*XX],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_4$X$^2$',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,7)
plt.boxplot([np.array(original[6])*YY,np.array(inner[6])*YY,np.array(orographic[6])*YY,np.array(full[6])*YY],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_5$Y$^2$',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,8)
plt.boxplot([np.array(original[7])*Z,np.array(inner[7])*Z,np.array(orographic[7])*Z,np.array(full[7])*Z],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_6$Z',fontsize=14)
plt.ylim(-40,40)
#plt.legend()
plt.tight_layout()
#plt.savefig('Precipitation_downscaling_coeffs_times_regressor_v2.pdf',bbox_inches='tight')

def normalizedXYZ(NARR_subregions):
    X = np.nanmean(Xgrid)/np.max(UTMx_list[NARR_subregions])
    Y = np.nanmean(Ygrid)/np.max(UTMy_list[NARR_subregions])
    XY = X*Y
    XX = X**2
    YY = Y**2
    Z = np.nanmean(Zgrid)/np.max(elev.ravel()[NARR_subregions])
    return X,Y,XY,XX,YY,Z

Xog,Yog,XYog,XXog,YYog,Zog = normalizedXYZ(NARR_subregions_original)
Xin,Yin,XYin,XXin,YYin,Zin = normalizedXYZ(NARR_subregions_inner)
Xor,Yor,XYor,XXor,YYor,Zor = normalizedXYZ(NARR_subregions_orographic)
Xfu,Yfu,XYfu,XXfu,YYfu,Zfu = normalizedXYZ(NARR_subregions_full)
    
plt.figure(figsize=(16,6))
plt.suptitle('Precipitation Regression Coefficients (Normalized X,Y,Z)',fontsize=14,y=1.02)
labels=['Original','Inner','Orographic','Full']
# 8 panels (1 per coefficient, w 5 boxes on each panel (one per subregions test))
plt.subplot(2,4,1)
plt.boxplot([original_norm[0],inner_norm[0],orographic_norm[0],full_norm[0]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('r$^2_{\\beta}$',fontsize=14)
plt.ylim(0,1.1)
plt.subplot(2,4,2)
plt.boxplot([original_norm[1],inner_norm[1],orographic_norm[1],full_norm[1]],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_0$',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,3)
plt.boxplot([np.array(original_norm[2])*Xog,np.array(inner_norm[2])*Xin,np.array(orographic_norm[2])*Xor,np.array(full_norm[2])*Xfu],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_1$X',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,4)
plt.boxplot([np.array(original_norm[3])*Yog,np.array(inner_norm[3])*Yin,np.array(orographic_norm[3])*Yor,np.array(full_norm[3])*Yfu],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_2$Y',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,5)
plt.boxplot([np.array(original_norm[4])*XYog,np.array(inner_norm[4])*XYin,np.array(orographic_norm[4])*XYor,np.array(full_norm[4])*XYfu],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_3$XY',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,6)
plt.boxplot([np.array(original_norm[5])*XXog,np.array(inner_norm[5])*XXin,np.array(orographic_norm[5])*XXor,np.array(full_norm[5])*XXfu],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_4$X$^2$',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,7)
plt.boxplot([np.array(original_norm[6])*YYog,np.array(inner_norm[6])*YYin,np.array(orographic_norm[6])*YYor,np.array(full_norm[6])*YYfu],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_5$Y$^2$',fontsize=14)
plt.ylim(-40,40)
plt.subplot(2,4,8)
plt.boxplot([np.array(original_norm[7])*Zog,np.array(inner_norm[7])*Zin,np.array(orographic_norm[7])*Zor,np.array(full_norm[7])*Zfu],labels=labels,showfliers=False,showmeans=True,meanprops=meanpointprops)
plt.ylabel('b$_6$Z',fontsize=14)
plt.ylim(-40,40)
#plt.legend()
plt.tight_layout()
#plt.savefig('Precipitation_downscaling_coeffs_times_regressor_normalizedcoeffs.pdf',bbox_inches='tight')
