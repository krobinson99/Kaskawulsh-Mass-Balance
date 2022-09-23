# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:39:48 2022

Read coordinate info from CoordinateTablev1.csv and map Rounce (2021) debris map
to the 200 m mass balance model domain

@author: katierobinson
"""

# steps to generate CoordinateTablev1.csv in ArcMap:
# 1. upload .tiff file containing debris thickness info
# 2. use Arc Toolbox Raster to Point tool
# 3. Open attribute table in ArcMap, add fields Lat / Lon
# 4. Open editor, select "calculate geometry"
# 5. calculate x,y coordinates for each column
# 6. export shapefile as .csv
# Adding elevation layer to the ArcMap:
# 1. Add Data > Add Basemap > Topographic World Map

import numpy as np
import os
import sys
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import regridXY_something
import utm
import matplotlib.pyplot as plt
import PIL
from PIL import Image


# READ IN COORDINATE TABLE (coords for Rounce et al. (2021) debris map)########
coordinate_file = 'CoordinateTablev1.csv' #this one contains only pixels with actual data (only debris covered cells)
coordstable = np.loadtxt(coordinate_file,delimiter=",",skiprows=1) #namelist!!
OG_debristhickness_list = coordstable[:,1]
lat = coordstable[:,2]
lon = coordstable[:,3]
x_meters = coordstable[:,4]
y_meters = coordstable[:,5]

# TRY TO MAKE ARRAY OF EASTING NORTHING FROM DR DEBRIS MAP #####################
thicknessmap = 'HMA_DTE_1.16201_hdts_m.tif'
debristhickness = Image.open(thicknessmap)
#debristhickness.show()
debristhickness_array = np.array(debristhickness)
debristhickness_array.shape

nanlocs = np.where(debristhickness_array >= 1e20 )
#zerolocs = np.where(debristhickness_array == 0 )
debristhickness_array[nanlocs] = 0

coordinate_file_v2 = 'CoordinateTablev2.csv' #this one contains data for ALL pixels
coordstable_v2 = np.loadtxt(coordinate_file_v2,delimiter=",",skiprows=1) #namelist!!
lat_v2 = coordstable_v2[:,4]
lon_v2 = coordstable_v2[:,3]

#covert to utm easting/northing
easting_DR_v2 = []
northing_DR_v2 =[]
for i in range(0,len(lat_v2)):
    x = utm.from_latlon(lat_v2[i],lon_v2[i])
    easting_DR_v2.append(x[0])
    northing_DR_v2.append(x[1])

#try filling the coords into the DR array TOP TO BOTTOM, LEFT TO RIGHT:
DR_Xgrid = np.reshape(easting_DR_v2,(debristhickness_array.shape))    
DR_Ygrid = np.reshape(northing_DR_v2,(debristhickness_array.shape)) 
   
# GET COORDINATES FOR MODEL DOMAIN#############################################
Path2Model = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel'
File_glacier_in = os.path.join(Path2Model,'Kaskonly_deb.txt')
glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
Ix = glacier[:,3] 
Iy = glacier[:,4] 
Ih = glacier[:,2] 
debris_array = glacier[:,6]

#-------Turn vectors into 3D gridded inputs--------------------

Zgrid, EY_Xgrid, EY_Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
nanlocs = np.where(np.isnan(Zgrid))

#Setup debris mask for use in radiation parameters
debris_grid, EY_Xgrid, EY_Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, debris_array)
debris_m = np.zeros(debris_grid.shape)
debris_m[np.where(debris_grid > 100)] = 1
debris_m[np.where(debris_grid <= 100)] = np.nan
debris_m[nanlocs] = np.nan

EY_Ygrid_flipped = np.flipud(EY_Ygrid) # EY_Ygrid is originally upside down (values decreasing south to north)

easting_EY = []
northing_EY = []
elevation_EY = []
for i in range(0,len(debris_m)):
    for j in range(0,len(debris_m[1])):
        if debris_m[i,j] == 1:
            easting_EY.append(EY_Xgrid[i,j])
            northing_EY.append(EY_Ygrid_flipped[i,j])
            elevation_EY.append(Zgrid[i,j])
        else:
            pass
        
# CONVERT EASTING/NORTHING TO LAT/LON using python module utm #################
easting_DR = []
northing_DR =[]
for i in range(0,len(lat)):
    x = utm.from_latlon(lat[i],lon[i])
    easting_DR.append(x[0])
    northing_DR.append(x[1])

    
# PLOT THE TWO DIFFERENT GRIDS! COMPARE
plt.figure(figsize=(8,5))
plt.scatter(easting_EY,northing_EY,color='red')
#plt.scatter(easting_EY,northing_EY,c=elevation_EY, cmap="RdYlBu")
#legend = plt.colorbar()
plt.scatter(easting_DR,northing_DR,color='dodgerblue')
plt.legend(['EY Debris Cells','DR Debris Cells'],fontsize=12)
plt.xlabel('Easting (m)',fontsize=12)
plt.ylabel('Northing (m)',fontsize=12)
#plt.savefig('EY_DR_debriscomparison.png')

DR_cellarea = 0.1*0.1
EY_cellarea = 0.2*0.2

DR_debrisarea = len(easting_DR)*DR_cellarea
EY_debrisarea = len(easting_EY)*EY_cellarea

print('Original DR debris area = ' + str(np.round(DR_debrisarea,2)) + ' km2')
print('Original EY debris area = ' + str(np.round(EY_debrisarea,2)) + ' km2')


# CALCULATE THE TOTAL VOLUME OF DEBRIS IN THE DR MAP ##########################
# Volume per cell = DR_cellarea * debris thickness
def CalculateDebrisVolume(debris_map,cellsize_m):
    
    Vol_per_cell = (cellsize_m*cellsize_m) * debris_map
    Total_deb_vol = np.nansum(Vol_per_cell)
    print('Total volume of debris in map = ' + str(Total_deb_vol) + ' m^3')
    return Total_deb_vol

CalculateDebrisVolume(OG_debristhickness_list,100)

# CALCULATE ELEVATION OF THE DR DEBRIS CELLS FROM EY ZGRID. ###################
# ASSIGN Z VALUE BASED ON NEAREST NEIGHBOUR EY CELL

# only let nearest neighbour be somewhere on the actual glacier (insert nanlocs into EY_Xgrid and EY_Ygrid)
EY_Xgrid[nanlocs] = np.nan #EY_Xgrid looks normal --> increases West to East
#EY_Ygrid[nanlocs] = np.nan # EY_Ygrid not normal --> increases South to North
EY_Ygrid_flipped[nanlocs] = np.nan #normal South to North increasing, but Kask is upside down, nanlocs are wrong. 

OG_elevation_DR = []
distance_to_NN = []
for i in range(0,len(easting_DR)):
    x_dist = EY_Xgrid - easting_DR[i]
    y_dist = EY_Ygrid_flipped - northing_DR[i]
    distance = np.sqrt((x_dist**2)+(y_dist**2))
    loc = np.where(distance == np.nanmin(distance))
    distance_to_NN.append(np.nanmin(distance))
    OG_elevation_DR.append(Zgrid[loc[0][0]][loc[1][0]])
# distance to nearest neighbour (NN) is 1-275 meters, mean = 78m away

# PLOT HISTOGRAM OF DEBRIS VOLUME VS ELEVATION ################################
# --> This is to check conservation of mass after interpolating to the 200 m domain
#working with: OG_elevation_DR and OG_debristhickness_list

#range is 700 - 2500 m (18 bins of 100 m each)
def DebrisVol_by_Elevation(debrislist,elevationlist,resolution_m):
    b700, b800, b900, b1000, b1100, b1200, b1300, b1400, b1500, b1600, b1700, b1800, b1900, \
    b2000, b2100, b2200, b2300, b2400 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    for i in range(0,len(elevationlist)):
        z = elevationlist[i]
        if z < 1600:
            if z >= 700 and z < 800:
                b700.append(debrislist[i])
            elif z >= 800 and z < 900:
                b800.append(debrislist[i])
            elif z >= 900 and z < 1000:
                b900.append(debrislist[i])
            elif z >= 1000 and z < 1100:
                b1000.append(debrislist[i])
            elif z >= 1100 and z < 1200:
                b1100.append(debrislist[i])
            elif z >= 1200 and z < 1300:
                b1200.append(debrislist[i])
            elif z >= 1300 and z < 1400:
                b1300.append(debrislist[i])
            elif z >= 1400 and z < 1500:
                b1400.append(debrislist[i])
            elif z >= 1500 and z < 1600:
                b1500.append(debrislist[i])
        elif z >= 1600:
            if z >= 1600 and z < 1700:
                b1600.append(debrislist[i])
            elif z >= 1700 and z < 1800:
                b1700.append(debrislist[i])
            elif z >= 1800 and z < 1900:
                b1800.append(debrislist[i])
            elif z >= 1900 and z < 2000:
                b1900.append(debrislist[i])
            elif z >= 2000 and z < 2100:
                b2000.append(debrislist[i])
            elif z >= 2100 and z < 2200:
                b2100.append(debrislist[i])
            elif z >= 2200 and z < 2300:
                b2200.append(debrislist[i])
            elif z >= 2300 and z < 2400:
                b2300.append(debrislist[i])
            elif z >= 2400 and z <= 2500:
                b2400.append(debrislist[i])
    
    #check that all cells were accounted for:
    numcells = np.sum([len(b700),len(b800),len(b900),len(b1000),len(b1100),len(b1200),len(b1300) \
                       ,len(b1400),len(b1500),len(b1600),len(b1700),len(b1800),len(b1900)\
                       ,len(b2000),len(b2100),len(b2200),len(b2300),len(b2400)])
    
    mass_per_bin = ([np.sum(b700),np.sum(b800),np.sum(b900),np.sum(b1000),np.sum(b1100),np.sum(b1200)\
                       ,np.sum(b1300),np.sum(b1400),np.sum(b1500),np.sum(b1600),np.sum(b1700),np.sum(b1800)\
                       ,np.sum(b1900),np.sum(b2000),np.sum(b2100),np.sum(b2200),np.sum(b2300),np.sum(b2400)]) #units = m
        
    volume_per_bin = np.array(mass_per_bin) * (resolution_m*resolution_m) / 1e9 #units = km^3
    
    mean_thickness_per_bin = ([np.mean(b700),np.mean(b800),np.mean(b900),np.mean(b1000),np.mean(b1100),np.mean(b1200)\
                       ,np.mean(b1300),np.mean(b1400),np.mean(b1500),np.mean(b1600),np.mean(b1700),np.mean(b1800)\
                       ,np.mean(b1900),np.mean(b2000),np.mean(b2100),np.mean(b2200),np.mean(b2300),np.mean(b2400)])
    
    i = 700
    zlabels = []
    while i < 2500:
        label = str(i+50) 
        zlabels.append(label)
        i += 100  
    
    return volume_per_bin,mean_thickness_per_bin,zlabels
        
# plot VOLUME of debris in each elevation band AND AVERAGE debris thickness at each elevation + range 
# also plot histogram of debris thicknesses? ie. debris covered-area vs thickness bin  
OG_DR_volume_per_bin = DebrisVol_by_Elevation(OG_debristhickness_list,OG_elevation_DR,100)[0]
OG_DR_thickness_per_bin = DebrisVol_by_Elevation(OG_debristhickness_list,OG_elevation_DR,100)[1]
zlabels = DebrisVol_by_Elevation(OG_debristhickness_list,OG_elevation_DR,100)[2]


x = np.arange(len(zlabels))
width = 0.75
plt.figure(figsize=(10,6))
plt.suptitle('Original DR Debris Map, Elevation interpolated from EY DEM',fontsize=14,y=1.01)
plt.subplot(1,2,1)
plt.title('Volume of Debris vs Elevation', fontsize=14)
#plt.bar(x-width, tmean_shift30, width,color='gold')
plt.barh(x, OG_DR_volume_per_bin, width,color='turquoise')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.yticks(x,zlabels, fontsize=14)
plt.ylabel('Elevation (m)',fontsize=14)
plt.xlabel('Volume of Debris (km$^3$)',fontsize=14)
plt.subplot(1,2,2)
plt.title('Mean Debris Thickness vs Elevation', fontsize=14)
#plt.bar(x-width, tmean_shift30, width,color='gold')
plt.barh(x, OG_DR_thickness_per_bin, width,color='turquoise')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.yticks(x,zlabels, fontsize=14)
plt.ylabel('Elevation (m)',fontsize=14)
plt.xlabel('Mean Debris Thickness (m)',fontsize=14)
plt.tight_layout()
#plt.savefig('DRdebris_volume_thickness_vs_z.png',bbox_inches='tight')

def debris_area_histogram(debrislist,resolution_km):
    
    #calculate histogram of debris thicknesses (ie. y-axis = thickness, x-axis = covered area)
    full_deb_hist = np.histogram(debrislist,bins=20)
    full_deb_labels = np.round(full_deb_hist[1][1:],2)
    full_deb_area = full_deb_hist[0]*(resolution_km*resolution_km)

    partial_deb_hist = np.histogram(debrislist,range=(0,0.1))
    partial_deb_labels = np.round(partial_deb_hist[1][1:],2)
    partial_deb_area = partial_deb_hist[0]*(resolution_km*resolution_km)
    
    return full_deb_area,full_deb_labels,partial_deb_area,partial_deb_labels

OG_DR_fulldebrisarea = debris_area_histogram(OG_debristhickness_list,0.1)[0]
OG_DR_fulldebrislabels = debris_area_histogram(OG_debristhickness_list,0.1)[1]
OG_DR_partialdebrisarea = debris_area_histogram(OG_debristhickness_list,0.1)[2]
OG_DR_partialdebrislabels = debris_area_histogram(OG_debristhickness_list,0.1)[3]


x2 = np.arange(len(OG_DR_fulldebrislabels))
x3 = np.arange(len(OG_DR_partialdebrislabels))

plt.figure(figsize=(10,6))
plt.suptitle('Original DR Debris Map \n Area vs Debris Thickness', fontsize=14,y=1.1)
plt.subplot(1,2,1)
#plt.bar(x-width, tmean_shift30, width,color='gold')
plt.barh(x2, OG_DR_fulldebrisarea, width,color='turquoise')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.yticks(x2,OG_DR_fulldebrislabels, fontsize=14)
plt.ylabel('Debris Thickness (m)',fontsize=14)
plt.xlabel('Area (km$^2$)',fontsize=14)
plt.subplot(1,2,2)
#plt.bar(x-width, tmean_shift30, width,color='gold')
plt.barh(x3, OG_DR_partialdebrisarea, width,color='turquoise')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.yticks(x3,OG_DR_partialdebrislabels, fontsize=14)
plt.ylabel('Debris Thickness (m)',fontsize=14)
plt.xlabel('Area (km$^2$)',fontsize=14)
plt.tight_layout()
plt.savefig('DRdebris_thickness_vs_area.png',bbox_inches='tight')

def Get_debris_lists_for_plotting(debrisarray,elevationarray):
    # make list of elevation and debris
    debristhickness_list = []
    debriselevation_list = []
    for i in range(0,len(debrisarray)):
        for j in range(0,len(debrisarray[0])):
            if np.isnan(debrisarray[i][j]):
                #print('naaan')
                pass
            else:
                debristhickness_list.append(debrisarray[i][j])
                debriselevation_list.append(elevationarray[i][j])
                
    return debristhickness_list, debriselevation_list
    
#now that this data is ready for quality checking, I can try a few different interpolations
# (100m grid to 200m grid) and check that they match the data well. 

# INTERPOLATION METHOD 1: NEAREST NEIGHBOUR ###################################
def NearestNeighbourInterp(original_domain='CoordinateTablev2.csv',target_domain=File_glacier_in,debthickness_array=thicknessmap):
    """
    this function takes in the original domain (coordinate table with debris cells, lat, long)
    (100m resolution debris map from Rounce et al. (2021)) and returns an "upscaled"
    debris map with same domain as 'target domain' (200 m).
    """
    #load the tiff file with debris thicknessess
    debristhickness = Image.open(debthickness_array)
    debristhickness_array = np.array(debristhickness)    
    nanlocs = np.where(debristhickness_array >= 1e20 )
    #zerolocs = np.where(debristhickness_array == 0 )
    debristhickness_array[nanlocs] = 0
    
    # load the original (100 m) domain
    coordstable_v2 = np.loadtxt(original_domain,delimiter=",",skiprows=1) #namelist!!
    lat_v2 = coordstable_v2[:,4]
    lon_v2 = coordstable_v2[:,3]
    
    #covert to utm easting/northing
    easting_DR_v2 = []
    northing_DR_v2 =[]
    for i in range(0,len(lat_v2)):
        x = utm.from_latlon(lat_v2[i],lon_v2[i])
        easting_DR_v2.append(x[0])
        northing_DR_v2.append(x[1])

    #fill the coords into the DR array TOP TO BOTTOM, LEFT TO RIGHT:
    DR_Xgrid = np.reshape(easting_DR_v2,(debristhickness_array.shape))    
    DR_Ygrid = np.reshape(northing_DR_v2,(debristhickness_array.shape)) 
        
    # load the target (200 m) domain
    glacier = np.genfromtxt(target_domain, skip_header=1, delimiter=',')
    Ix = glacier[:,3] 
    Iy = glacier[:,4] 
    Ih = glacier[:,2] 
    
    #-------Turn vectors into 3D gridded inputs--------------------
    
    Zgrid, EY_Xgrid, EY_Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
    nanlocs = np.where(np.isnan(Zgrid))
    
    EY_Ygrid_flipped = np.flipud(EY_Ygrid) # EY_Ygrid is originally upside down (values decreasing south to north)
    EY_Xgrid[nanlocs] = np.nan
    EY_Ygrid_flipped[nanlocs] = np.nan
    
    NN_map = np.empty(EY_Xgrid.shape)
    distance_to_NN = []
    for i in range(0,len(EY_Xgrid)):
        for j in range(0,len(EY_Xgrid[0])):
            if np.isnan(EY_Xgrid[i][j]):
                pass
            else:
                x_dist = DR_Xgrid - EY_Xgrid[i][j]
                y_dist = DR_Ygrid - EY_Ygrid_flipped[i][j]
                distance = np.sqrt((x_dist**2)+(y_dist**2))
                closest_cell = np.where(distance == np.nanmin(distance))
                distance_to_NN.append(np.nanmin(distance))
                nearestneighbour_debris = debristhickness_array[closest_cell]
                #try adding a contigency that if the 1st closest neighbour is NaN, look for the second closest neighbour:
                if nearestneighbour_debris == 0:
                    distance[closest_cell] = 9999
                    closest_cell = np.where(distance == np.nanmin(distance)) #new closest cell
                    print('searching for second closest cell: distance = ' + str(np.nanmin(distance)))
                    if np.nanmin(distance) < 75:
                        nearestneighbour_debris = debristhickness_array[closest_cell]
                        print('Pass! second closest cell: distance = ' + str(np.nanmin(distance)))
                    else:
                        nearestneighbour_debris = 0
                        print('Fail! second closest cell: distance = ' + str(np.nanmin(distance)))
                    
                    
                NN_map[i][j] = nearestneighbour_debris
    
    NN_map[nanlocs] = np.nan
        
    return NN_map,distance_to_NN
        
NNdebmap = NearestNeighbourInterp()[0]
nodeblocs = np.where(NNdebmap == 0)
NNdebmap[nodeblocs] = np.nan

distance_to_nearneighbour =  NearestNeighbourInterp()[1]

plt.figure(figsize=(8,5))
plt.contourf(np.flipud(NNdebmap[:180,100:]), cmap = 'PuOr', levels = np.round(np.linspace(0,1,10),1))
legend = plt.colorbar()
legend.ax.set_ylabel('Debris Thickness (m)', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
plt.title('Nearest Neighbour Resampling',fontsize=14)
#plt.savefig('NN_resampled_debris.png',bbox_inches = 'tight')

volumeratio = CalculateDebrisVolume(OG_debristhickness_list,100)/CalculateDebrisVolume(NNdebmap,200)

NN_debristhickness_list = Get_debris_lists_for_plotting(NNdebmap,Zgrid)[0]
NN_debriselevation_list = Get_debris_lists_for_plotting(NNdebmap,Zgrid)[1]

NN_volume_per_bin = DebrisVol_by_Elevation(NN_debristhickness_list,NN_debriselevation_list,200)[0]
NN_thickness_per_bin = DebrisVol_by_Elevation(NN_debristhickness_list,NN_debriselevation_list,200)[1]

NN_fulldebrisarea = debris_area_histogram(NN_debristhickness_list,0.2)[0]
NN_partialdebrisarea = debris_area_histogram(NN_debristhickness_list,0.2)[2]


# INTERPOLATION METHOD 2: INVERSE DISTANCE WEIGHTED AVERAGE ###################
def InverseDistance_Interp(searchradius,p,original_domain='CoordinateTablev2.csv',target_domain=File_glacier_in,debthickness_array=thicknessmap):
    """
    this function uses inverse distane weighted average to interpolate the DR debris thickness map to the
    200 m mass balance model domain, using the formula: https://en.wikipedia.org/wiki/Inverse_distance_weighting.
    
    The value of each cell in the model domain is obtained by weighing the DR debris thickness values,
    NaNs are treated as zeros (ie. clean ice). 
    """
    #load the tiff file with debris thicknessess
    debristhickness = Image.open(debthickness_array)
    debristhickness_array = np.array(debristhickness)    
    nanlocs = np.where(debristhickness_array >= 1e20 )
    #zerolocs = np.where(debristhickness_array == 0 )
    debristhickness_array[nanlocs] = 0
    
    # load the original (100 m) domain
    coordstable_v2 = np.loadtxt(original_domain,delimiter=",",skiprows=1) #namelist!!
    lat_v2 = coordstable_v2[:,4]
    lon_v2 = coordstable_v2[:,3]
    
    #covert to utm easting/northing
    easting_DR_v2 = []
    northing_DR_v2 =[]
    for i in range(0,len(lat_v2)):
        x = utm.from_latlon(lat_v2[i],lon_v2[i])
        easting_DR_v2.append(x[0])
        northing_DR_v2.append(x[1])

    #fill the coords into the DR array TOP TO BOTTOM, LEFT TO RIGHT:
    DR_Xgrid = np.reshape(easting_DR_v2,(debristhickness_array.shape))    
    DR_Ygrid = np.reshape(northing_DR_v2,(debristhickness_array.shape)) 
        
    # load the target (200 m) domain
    glacier = np.genfromtxt(target_domain, skip_header=1, delimiter=',')
    Ix = glacier[:,3] 
    Iy = glacier[:,4] 
    Ih = glacier[:,2] 
    
    #-------Turn vectors into 3D gridded inputs--------------------
    
    Zgrid, EY_Xgrid, EY_Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
    nanlocs = np.where(np.isnan(Zgrid))
    
    EY_Ygrid_flipped = np.flipud(EY_Ygrid) # EY_Ygrid is originally upside down (values decreasing south to north)
    EY_Xgrid[nanlocs] = np.nan
    EY_Ygrid_flipped[nanlocs] = np.nan
    
    #get all nearest neighbours within the search radius
    IDWA_map = np.empty(EY_Xgrid.shape)
    for i in range(0,len(EY_Xgrid)):
        for j in range(0,len(EY_Xgrid[0])):
            if np.isnan(EY_Xgrid[i][j]):
                pass
            else:
                x_dist = DR_Xgrid - EY_Xgrid[i][j]
                y_dist = DR_Ygrid - EY_Ygrid_flipped[i][j]
                distance = np.sqrt((x_dist**2)+(y_dist**2))
                within_search_radius = np.where(distance <= searchradius)
                
                ui = debristhickness_array[within_search_radius] #deb thicknesses within search radius
                dx = distance[within_search_radius] #distances within search radius
                
                wi = 1/((dx)**p)
                ux = (np.sum(wi*ui))/(np.sum(wi))
                
                IDWA_map[i,j] = ux
                
    IDWA_map[nanlocs] = np.nan

    return IDWA_map

IDWAdebmap = InverseDistance_Interp(70,1)

nodeblocs = np.where(IDWAdebmap == 0)
IDWAdebmap[nodeblocs] = np.nan

diff = NNdebmap-IDWAdebmap

plt.figure(figsize=(8,5))
#plt.contourf(debristhickness_array, cmap = 'PuOr', levels = np.round(np.linspace(0,3,10),1))
plt.contourf(np.flipud(IDWAdebmap[:180,100:]), cmap = 'PuOr',levels = np.round(np.linspace(0.0001,1,10),3))
legend = plt.colorbar()
legend.ax.set_ylabel('Debris Thickness (m)', rotation=270,fontsize=14,labelpad=20)
plt.tight_layout()
plt.title('Inverse Distance Weighted Average Resampling',fontsize=14)
#plt.savefig('IDWA_resampled_debris.png',bbox_inches = 'tight')

IDWA_debristhickness_list = Get_debris_lists_for_plotting(IDWAdebmap,Zgrid)[0]
IDWA_debriselevation_list = Get_debris_lists_for_plotting(IDWAdebmap,Zgrid)[1]

IDWA_volume_per_bin = DebrisVol_by_Elevation(IDWA_debristhickness_list,IDWA_debriselevation_list,200)[0]
IDWA_thickness_per_bin = DebrisVol_by_Elevation(IDWA_debristhickness_list,IDWA_debriselevation_list,200)[1]

IDWA_fulldebrisarea = debris_area_histogram(IDWA_debristhickness_list,0.2)[0]
IDWA_partialdebrisarea = debris_area_histogram(IDWA_debristhickness_list,0.2)[2]

# INTERPOLATION METHOD 3: SPLINE INTERPOLATION ###############################


###############################################################################




# MASTER COMPARISON PLOT !!
x = np.arange(len(zlabels))
width = 0.25
plt.figure(figsize=(10,6))
plt.suptitle('Original DR Debris Map vs Nearest Neighbour Resampling',fontsize=14,y=1.01)
plt.subplot(1,2,1)
plt.title('Volume of Debris vs Elevation', fontsize=14)
#plt.bar(x-0, tmean_shift30, width,color='gold')
plt.barh(x-width, OG_DR_volume_per_bin, width,color='turquoise')
plt.barh(x+width, NN_volume_per_bin, width,color='red')
plt.barh(x, IDWA_volume_per_bin, width, color='orange')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.legend(['Orignal debris map (Rounce et al.)','Nearest neighbour resampling','Inverse distance weighted \n average resampling'])
plt.yticks(x,zlabels, fontsize=14)
plt.ylabel('Elevation (m)',fontsize=14)
plt.xlabel('Volume of Debris (km$^3$)',fontsize=14)
plt.subplot(1,2,2)
plt.title('Mean Debris Thickness vs Elevation', fontsize=14)
#plt.bar(x-width, tmean_shift30, width,color='gold')
plt.barh(x-width, OG_DR_thickness_per_bin, width,color='turquoise')
plt.barh(x+width, NN_thickness_per_bin, width,color='red')
plt.barh(x, IDWA_thickness_per_bin, width, color='orange')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.yticks(x,zlabels, fontsize=14)
plt.ylabel('Elevation (m)',fontsize=14)
plt.xlabel('Mean Debris Thickness (m)',fontsize=14)
plt.tight_layout()
#plt.savefig('resampleddebris_volume_thickness_vs_z.png',bbox_inches='tight')

x2 = np.arange(len(OG_DR_fulldebrislabels))
x3 = np.arange(len(OG_DR_partialdebrislabels))
width = 0.25
plt.figure(figsize=(10,6))
plt.suptitle('Area vs Debris Thickness', fontsize=14,y=1.1)
plt.subplot(1,2,1)
#plt.bar(x-width, tmean_shift30, width,color='gold')
plt.barh(x2-width, OG_DR_fulldebrisarea, width,color='turquoise')
plt.barh(x2+width, NN_fulldebrisarea, width,color='red')
plt.barh(x2, IDWA_fulldebrisarea, width, color='orange')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.legend(['Orignal debris map (Rounce et al.)','Nearest neighbour resampling','Inverse distance weighted \n average resampling'])
plt.yticks(x2,OG_DR_fulldebrislabels, fontsize=14)
plt.ylabel('Debris Thickness (m)',fontsize=14)
plt.xlabel('Area (km$^2$)',fontsize=14)
plt.subplot(1,2,2)
#plt.bar(x-width, tmean_shift30, width,color='gold')
plt.barh(x3-width, OG_DR_partialdebrisarea, width,color='turquoise')
plt.barh(x3+width, NN_partialdebrisarea, width,color='red')
plt.barh(x3, IDWA_partialdebrisarea, width, color='orange')
#plt.bar(x+width, tmean_shift90, width,color='crimson')
plt.legend(['Orignal debris map (Rounce et al.)','Nearest neighbour resampling','Inverse distance weighted \n average resampling'])
plt.yticks(x3,OG_DR_partialdebrislabels, fontsize=14)
plt.ylabel('Debris Thickness (m)',fontsize=14)
plt.xlabel('Area (km$^2$)',fontsize=14)
plt.tight_layout()
#plt.savefig('resampleddebris_thicknessvsarea.png',bbox_inches='tight')

#PLOT THE SPATIAL EXTENT OF THE RESAMPLED DEBRIS MAPS

easting_NN = []
northing_NN = []
elevation_NN = []
for i in range(0,len(NNdebmap)):
    for j in range(0,len(NNdebmap[1])):
        if np.isnan(NNdebmap[i,j]):
            pass
        else:
            easting_NN.append(EY_Xgrid[i,j])
            northing_NN.append(EY_Ygrid_flipped[i,j])
            elevation_NN.append(Zgrid[i,j])

# PLOT THE TWO DIFFERENT GRIDS! COMPARE
plt.figure(figsize=(8,5))
plt.scatter(easting_EY,northing_EY,color='k',s=10)
#plt.scatter(easting_EY,northing_EY,c=elevation_EY, cmap="RdYlBu")
#legend = plt.colorbar()
plt.scatter(easting_DR,northing_DR,color='turquoise',s=10)
plt.scatter(easting_NN,northing_NN,color='red',s=1)
plt.legend(['EY Debris Cells','DR Debris Cells','NN resampled cells'],fontsize=12)
plt.xlabel('Easting (m)',fontsize=12)
plt.ylabel('Northing (m)',fontsize=12)
#plt.savefig('EY_DR_NN_debriscomparison.png')

easting_IDWA = []
northing_IDWA = []
elevation_IDWA = []
for i in range(0,len(IDWAdebmap)):
    for j in range(0,len(IDWAdebmap[1])):
        if np.isnan(IDWAdebmap[i,j]):
            pass
        else:
            easting_IDWA.append(EY_Xgrid[i,j])
            northing_IDWA.append(EY_Ygrid_flipped[i,j])
            elevation_IDWA.append(Zgrid[i,j])

# PLOT THE TWO DIFFERENT GRIDS! COMPARE
plt.figure(figsize=(8,5))
plt.scatter(easting_EY,northing_EY,color='k',s=10)
#plt.scatter(easting_EY,northing_EY,c=elevation_EY, cmap="RdYlBu")
#legend = plt.colorbar()
plt.scatter(easting_DR,northing_DR,color='turquoise',s=10)
plt.scatter(easting_IDWA,northing_IDWA,color='orange',s=10)
plt.legend(['EY Debris Cells','DR Debris Cells','IDWA resampled cells (radius = 70m)'],fontsize=12)
plt.xlabel('Easting (m)',fontsize=12)
plt.ylabel('Northing (m)',fontsize=12)
#plt.savefig('EY_DR_IDWA_debriscomparison.png')


# make boolean version of the chosen debris map ########
#which map is better? NN or IDWA:
#answer = NN (IDWA has holes in debris cover at the terminus)
def booleandebmap(debmap):
    '''
    make a debris map that has the same exact shape and debris covered cells
    as the variable thickness map, except all debris covered cells are boolean
    0 (debris) or 1 (no debris)
    '''
    
    debmapbool = np.ones(debmap.shape)
    deblocs = np.where(debmap >= 0)
    debmapbool[deblocs] = 0
    
    return debmapbool

booleanmap = booleandebmap(NNdebmap)

#np.save('debmap_variableh.npy',NNdebmap)
#np.save('debmap_boolean.npy',booleanmap)

