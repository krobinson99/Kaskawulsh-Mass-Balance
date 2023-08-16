import numpy as np
#from scipy.spatial import distance
from sklearn import linear_model
#from scipy.stats import linregress
#import math
from scipy import interpolate
from netCDF4 import Dataset
import netCDF4
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd


###Downscale temperature based on atmospheric structure and inversions###
def T_downscale_funkfest(T, E, UTMx_list, UTMy_list):
        '''
        E has shape 29,6,6 (pressure level, x,y)
        '''
        
        ###Get values for geopotential and temperature
        # KR_note: these would probably work just as well as 6,6 arrays
        xi_list = np.empty_like(T[-1][:][:]).tolist()
        yi_list = np.empty_like(T[-1][:][:]).tolist()
        xi_list_inver = np.empty_like(T[-1][:][:]).tolist()
        yi_list_inver = np.empty_like(T[-1][:][:]).tolist()        
        funclist = np.empty_like(T[-1][:][:]).tolist()
        funclist_inver = np.empty_like(T[-1][:][:]).tolist()
        inversion_list = np.empty_like(T[-1][:][:]).tolist()
        
        L_list = np.empty_like(T[-1][:][:]).tolist()
        L_list_inver = np.empty_like(T[-1][:][:]).tolist()
        y0_list = np.empty_like(T[-1][:][:]).tolist()
        y0_list_inver = np.empty_like(T[-1][:][:]).tolist()
        
        #set threshold for dropping high altitude values
        th = -9

        #interpolate elev and temp, get interp function for each grid point
        #also save T and geopotential points for extrapolation if needed   
        for u in range(0,T[0].shape[0]):
            for w in range(0, T[0].shape[0]): 
                inversion_list[u][w] = 0
                j = 0
                
                #scan for inversions in reanalysis data, stop at height v to 
                #avoid including upper atm inversions
                while j == 0:
                    for v in range(2, len(T)-1):
                        if v > 10:
                            j = 1
                        elif T[v+1][u][w] < T[v-2][u][w]:
                            j = 0
                        else:
                            inversion_list[u][w] = v+1
                            j = 1
                            #print u
                            #print w
                            #print v
                            
                
                #If no inversions present, use full values distribution. Save 
                #linear polyline function for this distribution
                if inversion_list[u][w] == 0:
                    xi = T[0:th, u, w]
                    yi = E[0:th, u, w]
                    function_param = np.polyfit(yi, xi, 1) 
                    function = np.poly1d(function_param)
                    funclist[u][w] = function
                    funclist_inver[u][w] = np.nan
                    xi_list[u][w] = xi
                    yi_list[u][w] = yi
                
                #If inversion present, use its position to truncate into two 
                #value distributions {upper and lower}
                else:
                    yi_up = E[inversion_list[u][w]:th, u, w]
                    yi_down = E[0:inversion_list[u][w], u, w]
                    xi_up = T[inversion_list[u][w]:th, u, w]
                    xi_down = T[0:inversion_list[u][w], u, w]
                    
                    function_up_param = np.polyfit(yi_up , xi_up, 1)
                    function_up = np.poly1d(function_up_param)
                    function_down_param = np.polyfit(yi_down , xi_down, 1) 
                    function_down = np.poly1d(function_down_param)
                    funclist[u][w] = function_up
                    funclist_inver[u][w] = function_down
                    xi_list[u][w] = xi_up
                    yi_list[u][w] = yi_up      
                    xi_list_inver[u][w] = xi_down
                    yi_list_inver[u][w] = yi_down 
                    
         
        #Get spatially interpolated lapse rates and intercepts
        for u in range(0, T[-1][:][:].shape[0]):
            for w in range(0, T[-1][:][:].shape[1]):
                
                L_list[u][w] = (funclist[u][w](4000) - funclist[u][w](2000))/(4000 - 2000)
                
                if np.nanmean(np.isnan(funclist_inver[u][w])) == 1:
                    L_list_inver[u][w] = L_list[u][w]
                else:
                    L_list_inver[u][w] = (funclist_inver[u][w](4000) - funclist_inver[u][w](2000))/(4000 - 2000)
                    
                y0_list[u][w] = funclist[u][w](0)
                
                if np.nanmean(np.isnan(funclist_inver[u][w])) == 1:
                    y0_list_inver[u][w] = y0_list[u][w]
                else:
                    y0_list_inver[u][w] = funclist_inver[u][w](0)
          
        y0func = interpolate.bisplrep(UTMy_list, UTMx_list, y0_list)
        y0func_inver = interpolate.bisplrep(UTMy_list, UTMx_list, y0_list_inver)
        Lfunc = interpolate.bisplrep(UTMy_list, UTMx_list, L_list)
        Lfunc_inver = interpolate.bisplrep(UTMy_list, UTMx_list, L_list_inver)
                          
                    
        return      xi_list, yi_list, xi_list_inver, yi_list_inver, \
                        funclist, funclist_inver, inversion_list, y0func, \
                            y0func_inver, Lfunc, Lfunc_inver
                        
                        
#convert liquid precip to snow based on T###  
def Precip_2_Snow(P,T,snowfac):
            
        #calculate snow
        if 273.75 >= T:
            S = P * snowfac         #use kienzle discussed factor of snow accumulation x10 precip
        elif 273.75 < T < 276.75:
            S = P*(1-((((T-273.15)/3)-0.2))) * snowfac #trouble!! use T in celcius for T/3
        else:
            S = 0.  
            
        return S
        
def Precip_2_Snow_thresh(P,T,snowfac, thresh):
            
        #calculate snow
        if thresh >= T:
            S = P * snowfac         #use kienzle discussed factor of snow accumulation x10 precip
        else:
            S = 0.  
            
        return S
        
        
        
        ###melt model for glaciers###
def MB_341(Thour, Phour, Pmelt, Pmelt_SRarc, Pmelt_SRre, SRhour_arc, SRhour_re, MFi, MFs, asnow, aice, MF):
    

        Melt_list = []
        Melt_listSR_arc = []
        Melt_listSR_re = []
        
        Leftover_list = []
        Leftover_listSR_arc = []
        Leftover_listSR_re = []
        
        MBhour = []
        MBhour_SRarc = []
        MBhour_SRre = []

        #DD melt model and snow tracker        
        for t in range(0, len(Thour)):
            if Thour[t] <= 0:
                Melt = 0.
                
                Melt_list.append(Melt)
                Melt_listSR_arc.append(Melt)
                Melt_listSR_re.append(Melt)
                
                snow = Pmelt[t]
                snowSRarc = Pmelt_SRarc[t]
                snowSRre = Pmelt_SRre[t]    
                
                Leftover_list.append(snow)
                Leftover_listSR_arc.append(snowSRarc)
                Leftover_listSR_re.append(snowSRre)
                               
            else:
                snow = Pmelt[t]
                snowSRarc = Pmelt_SRarc[t]
                snowSRre = Pmelt_SRre[t]  
                
                Melt_snow = MFs * Thour[t]
                Melt_snowSR_arc = (MF + asnow * SRhour_arc[t]) * Thour[t]
                Melt_snowSR_re = (MF + asnow * SRhour_re[t]) * Thour[t]
                
                #DD model
                if Melt_snow < snow:
                    Melt = Melt_snow
                    Melt_list.append(Melt)
                    Leftover = snow - Melt_snow
                    Leftover_list.append(Leftover)
                    
                else:
                    Melt_snow = snow
                    Tice = Melt_snow/MFs
                    Melt_ice = MFi * (Thour[t] - Tice)
                    Melt = Melt_snow + Melt_ice
                    Melt_list.append(Melt)
                    Leftover = 0.
                    Leftover_list.append(Leftover)
                
                #EDD model using I arc    
                if Melt_snowSR_arc < snowSRarc:
                    Melt = Melt_snowSR_arc
                    Melt_listSR_arc.append(Melt)
                    Leftover = snowSRarc - Melt_snowSR_arc
                    Leftover_listSR_arc.append(Leftover)
                    
                else:
                    Melt_snowSR_arc = snowSRarc
                    Tice = Melt_snowSR_arc/(MF + asnow * SRhour_arc[t])
                    Melt_iceSR_arc = (MF + aice * SRhour_arc[t]) * (Thour[t] - Tice)
                    Melt = Melt_snowSR_arc + Melt_iceSR_arc
                    Melt_listSR_arc.append(Melt)
                    Leftover = 0.
                    Leftover_listSR_arc.append(Leftover)
                    
                #EDD model using I re    
                if Melt_snowSR_re < snowSRre:
                    Melt = Melt_snowSR_re
                    Melt_listSR_re.append(Melt)
                    Leftover = snow - Melt_snowSR_re
                    Leftover_listSR_re.append(Leftover)
                    
                else:
                    Melt_snowSR_re = snowSRre
                    Tice = Melt_snowSR_re/(MF + asnow * SRhour_re[t])
                    Melt_iceSR_re = (MF + aice * SRhour_re[t]) * (Thour[t] - Tice)
                    Melt = Melt_snowSR_re + Melt_iceSR_re
                    Melt_listSR_re.append(Melt)
                    Leftover = 0.
                    Leftover_listSR_re.append(Leftover)
            
            
            
        Mass_in = Phour
        Mass_out = Melt_list
        MB = np.asarray(Mass_in) - np.asarray(Mass_out)
        Mass_in = Phour
        Mass_out = Melt_listSR_arc
        MB_SR_arc = np.asarray(Mass_in) - np.asarray(Mass_out)
        Mass_in = Phour
        Mass_out = Melt_listSR_re
        MB_SR_re = np.asarray(Mass_in) - np.asarray(Mass_out)
        
        MBhour.append(MB) 
        MBhour_SRarc.append(MB_SR_arc)
        MBhour_SRre.append(MB_SR_re)
        
        return MBhour, MBhour_SRarc, MBhour_SRre, Leftover_list, \
                Leftover_listSR_arc, Leftover_listSR_re, Melt_list, \
                    Melt_listSR_re, Melt_listSR_arc
                
                

                ###melt model for debris covered glaciers###
def MB_341_debris(Thour, Phour, Pmelt, Pmelt_SRarc, SRhour_arc, MFi, MFs, asnow, aice, MF, deb):
    

        Melt_list = []
        Melt_listSR_arc = []
        # Melt_listSR_re = []
        
        Leftover_list = []
        Leftover_listSR_arc = []
        # Leftover_listSR_re = []
        
        MBhour = []
        MBhour_SRarc = []
        # MBhour_SRre = []

        #DD melt model and snow tracker        
        for t in range(0, len(Thour)):
            if Thour[t] <= 0:
                Melt = 0.
                
                Melt_list.append(Melt)
                Melt_listSR_arc.append(Melt)
                # Melt_listSR_re.append(Melt)
                
                snow = Pmelt[t]
                snowSRarc = Pmelt_SRarc[t]
                # snowSRre = Pmelt_SRre[t]    
                
                Leftover_list.append(snow)
                Leftover_listSR_arc.append(snowSRarc)
                # Leftover_listSR_re.append(snowSRre)
                               
            else:
                snow = Pmelt[t]
                snowSRarc = Pmelt_SRarc[t]
                # snowSRre = Pmelt_SRre[t]  
                
                Melt_snow = MFs * Thour[t]
                Melt_snowSR_arc = (MF + asnow * SRhour_arc[t]) * Thour[t]
                # Melt_snowSR_re = (MF + asnow * SRhour_re[t]) * Thour[t]
                
                #DD model
                if Melt_snow < snow:
                    Melt = Melt_snow
                    Melt_list.append(Melt)
                    Leftover = snow - Melt_snow
                    Leftover_list.append(Leftover)
                    
                else:
                    Melt_snow = snow
                    Tice = Melt_snow/MFs
                    Melt_ice = MFi * (Thour[t] - Tice) * deb[t]
                    Melt = Melt_snow + Melt_ice
                    Melt_list.append(Melt)
                    Leftover = 0.
                    Leftover_list.append(Leftover)
                
                #EDD model using I arc    
                if Melt_snowSR_arc < snowSRarc:
                    Melt = Melt_snowSR_arc
                    Melt_listSR_arc.append(Melt)
                    Leftover = snowSRarc - Melt_snowSR_arc
                    Leftover_listSR_arc.append(Leftover)
                    
                else:
                    Melt_snowSR_arc = snowSRarc
                    Tice = Melt_snowSR_arc/(MF + asnow * SRhour_arc[t])
                    Melt_iceSR_arc = (MF + aice * SRhour_arc[t]) * (Thour[t] - Tice)  * deb[t]
                    Melt = Melt_snowSR_arc + Melt_iceSR_arc
                    Melt_listSR_arc.append(Melt)
                    Leftover = 0.
                    Leftover_listSR_arc.append(Leftover)
                    
                # #EDD model using I re    
                # if Melt_snowSR_re < snowSRre:
                #     Melt = Melt_snowSR_re
                #     Melt_listSR_re.append(Melt)
                #     Leftover = snow - Melt_snowSR_re
                #     Leftover_listSR_re.append(Leftover)
                #     
                # else:
                #     Melt_snowSR_re = snowSRre
                #     Tice = Melt_snowSR_re/(MF + asnow * SRhour_re[t])
                #     Melt_iceSR_re = (MF + aice * SRhour_re[t]) * (Thour[t] - Tice)  * deb[t]
                #     Melt = Melt_snowSR_re + Melt_iceSR_re
                #     Melt_listSR_re.append(Melt)
                #     Leftover = 0.
                #     Leftover_listSR_re.append(Leftover)
            
            
            
        Mass_in = Phour
        Mass_out = Melt_list
        MB = np.asarray(Mass_in) - np.asarray(Mass_out)
        Mass_in = Phour
        Mass_out = Melt_listSR_arc
        MB_SR_arc = np.asarray(Mass_in) - np.asarray(Mass_out)
        # Mass_in = Phour
        # Mass_out = Melt_listSR_re
        # MB_SR_re = np.asarray(Mass_in) - np.asarray(Mass_out)
        
        MBhour.append(MB) 
        MBhour_SRarc.append(MB_SR_arc)
        # MBhour_SRre.append(MB_SR_re)
        
        return MBhour, MBhour_SRarc, Leftover_list, \
                Leftover_listSR_arc, Melt_list, \
                 Melt_listSR_arc                                
                                                                
                
                            ###melt model for office snow covered surfaces###            
def MB_341_office(Thour, Phour, Pmelt, Pmelt_SRarc, SRhour_arc, MFi, MFs, asnow, aice, MF):

        Melt_list = []
        Melt_listSR_arc = []
        # Melt_listSR_re = []
        
        Leftover_list = []
        Leftover_listSR_arc = []
        # Leftover_listSR_re = []
        
        MBhour = []
        MBhour_SRarc = []
        # MBhour_SRre = []

        #DD melt model and snow tracker        
        for t in range(0, len(Thour)):
            if Thour[t] <= 0:
                Melt = 0.
                
                Melt_list.append(Melt)
                Melt_listSR_arc.append(Melt)
                # Melt_listSR_re.append(Melt)
                
                snow = Pmelt[t]
                snowSRarc = Pmelt_SRarc[t]
                # snowSRre = Pmelt_SRre[t]    
                
                Leftover_list.append(snow)
                Leftover_listSR_arc.append(snowSRarc)
                # Leftover_listSR_re.append(snowSRre)
                               
            else:
                snow = Pmelt[t]
                snowSRarc = Pmelt_SRarc[t]
                # snowSRre = Pmelt_SRre[t]  
                
                Melt_snow = MFs * Thour[t]
                Melt_snowSR_arc = (MF + asnow * SRhour_arc[t]) * Thour[t]
                # Melt_snowSR_re = (MF + asnow * SRhour_re[t]) * Thour[t]
                
                #DD model
                if Melt_snow < snow:
                    Melt = Melt_snow
                    Melt_list.append(Melt)
                    Leftover = snow - Melt_snow
                    Leftover_list.append(Leftover)
                    
                else:
                    Melt_snow = snow
                    Melt = Melt_snow
                    Melt_list.append(Melt)
                    Leftover = 0.
                    Leftover_list.append(Leftover)
                
                #EDD model using I arc    
                if Melt_snowSR_arc < snowSRarc:
                    Melt = Melt_snowSR_arc
                    Melt_listSR_arc.append(Melt)
                    Leftover = snowSRarc - Melt_snowSR_arc
                    Leftover_listSR_arc.append(Leftover)
                    
                else:
                    Melt_snowSR_arc = snowSRarc
                    Melt = Melt_snowSR_arc
                    Melt_listSR_arc.append(Melt)
                    Leftover = 0.
                    Leftover_listSR_arc.append(Leftover)
                    
                # #EDD model using I re    
                # if Melt_snowSR_re < snowSRre:
                #     Melt = Melt_snowSR_re
                #     Melt_listSR_re.append(Melt)
                #     Leftover = snow - Melt_snowSR_re
                #     Leftover_listSR_re.append(Leftover)
                #     
                # else:
                #     Melt_snowSR_re = snowSRre
                #     Melt = Melt_snowSR_re
                #     Melt_listSR_re.append(Melt)
                #     Leftover = 0.
                #     Leftover_listSR_re.append(Leftover)
            
            
            
        Mass_in = Phour
        Mass_out = Melt_list
        MB = np.asarray(Mass_in) - np.asarray(Mass_out)
        Mass_in = Phour
        Mass_out = Melt_listSR_arc
        MB_SR_arc = np.asarray(Mass_in) - np.asarray(Mass_out)
        # Mass_in = Phour
        # Mass_out = Melt_listSR_re
        # MB_SR_re = np.asarray(Mass_in) - np.asarray(Mass_out)
        
        MBhour.append(MB) 
        MBhour_SRarc.append(MB_SR_arc)
        # MBhour_SRre.append(MB_SR_re)
        
        return MBhour, MBhour_SRarc, Leftover_list, \
                Leftover_listSR_arc, Melt_list, \
                 Melt_listSR_arc
                
                
                

def regridXY_something(Ix, Iy, something):

    ###Regrid array into 3d, like meshgrid but supports third variable###
    
    #check for irregular gridding
    xr = np.unique(Ix)
    yr = np.unique(Iy)
    grid_spacex = np.diff(xr)
    grid_spacey = np.diff(yr)
    norm_gridx = np.mean(grid_spacex).round(-2)
    norm_gridy = np.mean(grid_spacey).round(-2) #Get grid spacing (= 200m)
    
    if np.mean(grid_spacex) == norm_gridx: #check that delta x = delta y (200 m each)
        pass
    else: # if not equal, make a correction (this whole loop is skipped since all cells have 200x200m res.)
        holes = np.where(grid_spacex > norm_gridx)
        d = grid_spacex[holes]
        vals = xr[holes]
        for val in range(0, len(vals)):
            v = vals[val]
            D = d[val]/norm_gridx
            for i in range(1, int(D)):
                new_val = v + (i * norm_gridx)
                xr = np.insert(xr, 0, new_val)
        
    if np.mean(grid_spacey) == norm_gridy: #they are equal, so this loop is also skipped
        pass
    else:
        holes = np.where(grid_spacey > norm_gridy)
        d = grid_spacey[holes]
        vals = yr[holes]
        for val in range(0, len(vals)):
            v = vals[val]
            D = d[val]/norm_gridy
            for i in range(1, int(D)):
                new_val = v + (i * norm_gridy)
                yr = np.insert(yr, 0, new_val)
        
        
    
    x = np.linspace(np.min(Ix), np.max(Ix), len(xr))
    y = np.linspace(np.min(Iy), np.max(Iy), len(yr))
    
    X, Y = np.meshgrid(x, y) #X is Xgrid, Y is Ygrid
        
    hold = np.empty(X.shape)
    hold[:] = np.nan
        
    XYZ = np.empty(X.shape)    
    XYZ[:] = np.ones(X.shape) * hold # an array of all nans the size of the domain
    
        
    for i in range(0, len(X)):

        for j in range(0, len(X[i])):

            locx = np.where(Ix == X[i,j])            
            locy = np.where(Iy[locx] == Y[i,j])
            loc = locx[0][locy]
            #print i, j
        
            if len(loc) > 0:
            
                XYZ[i,j] = something[loc][0]
            
            else:
                pass
                
    xyz = np.flip(XYZ, 0)
    #print(xyz)

                
    return xyz, X, Y, x, y


    
    

#get closest node to a point from a list of nodes
def closest_node(downscaledcell, NARRx, NARRy):
    
    x_dist = downscaledcell[0] - NARRx
    y_dist = downscaledcell[1] - NARRy
    dist = np.sqrt(x_dist**2 + y_dist**2)
    node_pos = np.argmin(dist)
    
    return node_pos
    

def precip_downscaling(elev, precip, lon, lat):

    #define list of relationships to be included in MLRM
    op_list = ['b0', 'b1X', 'b2Y', 'b3XY', 'b4X2', 'b5Y2', 'b6Z']

    #find rainy pixels
    rainy_px = np.where(precip > 0)
    
    #set variable list
    Y = lat[rainy_px]
    X = lon[rainy_px]
    XY = (lat[rainy_px])*(lon[rainy_px])
    Z = elev[rainy_px]
    Y2 = (lat[rainy_px])**2
    X2 = (lon[rainy_px])**2
    
    #check that enough pixels are rainy
    if len(rainy_px[0]) > len(op_list):
    
        clf = linear_model.LinearRegression()
        clf.fit(np.asarray((X,Y,XY,X2,Y2,Z)).transpose(), precip[rainy_px])
        
        W1 = clf.score(np.asarray((X,Y,XY,X2,Y2,Z)).transpose(), precip[rainy_px]) 
        coefs = clf.coef_
        I0 = clf.intercept_     
        
    else:
        W1 = 0.
        coefs = np.zeros(len(op_list))
        I0 = 0.
            
    return W1, coefs, I0
    

#bilinear interpolation functionn obtained from 
#https://stackoverflow.com/questions/8661537/how-to-perform-bilinear-interpolation-in-python
#python has trouble with proper bilinear interpolation instead opting for 
def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)
           

###create empty netcdf container, need ref files in directory
def save_to_netcdf(MB, var_n, File_name, year, Xgrid, Ygrid):

    ###Create empty .nc file for MB values###
    f = Dataset(File_name, 'w', format='NETCDF4') #'w' stands for write
        
    time = f.createDimension('time', len(MB))
    #z = f.createDimension('z', len(Ih))
    y = f.createDimension('y', len(MB[0]))
    x = f.createDimension('x', len(MB[0][0]))
    #f.createDimension('date', None)
    
    #load subset into new netcdffile
    times = f.createVariable('time', np.float64, 'time')
    #Zs = f.createVariable('z', np.int32, 'z')
    Ys = f.createVariable('y', np.float32, 'y')
    Xs = f.createVariable('x', np.float32, 'x')  
    mb = f.createVariable(var_n, np.float32, ('time', 'y', 'x'))
        
    # Global Attributes 
    f.description = 'Mass-Balance Model Input for Kaskawulsh Glacier (K. Robinson MSc Thesis)'   
    f.source = File_name
    # Variable Attributes  
    Ys.units = 'm'  
    Xs.units = 'm'  
    #Zs.units = 'm' 
    mb.units = 'm we' 
    times.units = 'hours since 1800-01-01 00:00:00' 
    times.calendar = 'gregorian' 
    times.delta_t = '0000-00-00 03:00:00'
    #time.actual_range = [1835688. , 1836405]
        
    #insert metadata into .nc file 
    XX = Xgrid[0]
    YY = Ygrid[:,0]
    #Z = Ih
    #time_nc = fh.variables['time'][:]
    #b = fh.variables['time']
    #dtime = netCDF4.num2date(b[:], b.units)
    dtime = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='3H').to_pydatetime()
    
    Xs[:] = XX #The "[:]" at the end of the variable instance is necessary
    Ys[:] = YY
    #Zs[:] = Z
    times[:] = netCDF4.date2num(dtime, units = times.units, calendar = times.calendar)
    mb[:] = np.zeros(mb[:].shape) * np.nan
    mb[:] = MB
    
    f.close()
    
    #return(mb[:], f)


def get_meanSP(year_list,Glacier_ID,R2S,Climate_inputs):
    '''
    Calculates the mean annual total accumulation
    '''
    DH_list = []
    for year in year_list:
        # Load Temp and Precip inputs for each year:
        inT = Dataset(os.path.join(Climate_inputs,'Temperature_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
        T_array = inT.variables['Temperature'][:]
        sys.stdout.flush()
        
        inP = Dataset(os.path.join(Climate_inputs,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
        P_array = inP.variables['Precipitation'][:]
        sys.stdout.flush()
        
        # Get indices for all rainy pixels:
        rainlocs = np.where(T_array > R2S)        
        
        # Get array with accumulation only (by zero-ing any pixels with rain instead of snow)
        Accumulation_array = np.array(P_array)
        Accumulation_array[rainlocs] = 0

        
        DH = np.sum(Accumulation_array, axis = 0)
        DH_list.append(DH)
        
        inT.close()
        inP.close()
        
        
    DHval = np.mean(DH_list, axis = 0)  # Think that DHval should be an array of shape (x,y)?   

    return DHval          
                                                
def cold_content(This_year_is, year_list, P_array, T_array, Glacier_id, DHval, BC_P,P_path):

#function called for calculating the cold content of each cell in the model
#domain. Applied for pre-processed model.

    #first define parameters
    c = 2009 #J kg^-1 K^-1
    L = 335000. #J kg^-1
    d = np.ones(P_array[0,:,:].shape)*2 #m
        
    #mean annual temperature
    Tval = np.mean(T_array, axis = 0)
    for t in range(0, len(Tval)):
        for tt in range(0, len(Tval[t])):
            if Tval[t][tt] > 0.:
                Tval[t][tt] = 0.
    
    #Total snow pack
    if BC_P == True:
        File_precip_in_future = 'Snow_' + 'kaskawulsh' + '_BC_' + str(This_year_is) + '.nc'
    else:
        File_precip_in_future = 'netSnow' + Glacier_id + str(This_year_is + 1) + '.nc'
    
    File_p_in = os.path.join(P_path,File_precip_in_future)
    
    infut = Dataset(File_p_in, "r")
    if BC_P == True:
        P_var = 'Snow'
    else:
        P_var = 'Temperature'
    P_array_fut = infut.variables[P_var][:]
    
    past = P_array[(int(-len(P_array_fut)/3)):,:,:]
    future = P_array_fut[:(int(len(P_array)*0.41)),:,:]
    SP = np.sum(past,axis = 0) + np.sum(future, axis = 0)
       
    #calculate snowpack to melt
    proportion = (c/L)*np.abs(Tval)*(d/DHval) 
    CC = proportion*SP
    
    return CC
    
def cold_content_simplified(year_range, MT, MSP, SPs):
    
    CCs = []
    ps = []
    
    #yearly snow pack
    for i in range(0, len(year_range)):
    
        #first define parameters
        c = 2009. #kJ kg^-1 K^-1
        L = 335000. #kJ kg^-1
        d = 2. #m
        
        proportion = (c/L)*(np.abs(MT[i]))*(d/MSP) 
        CC = proportion*SPs[i]
        
        CCs.append(CC)
        ps.append(proportion)
    
    return CCs, ps
    
def mean_snowpacks_pts(year_range):
    
    SPs = []
    
    #Total snow pack
    for i in year_range:
        
        precip = np.genfromtxt('netSnowtuning' + str(i) + '.txt', delimiter = ',')
        total_precip = np.sum(precip, axis = 0)
        SPs.append(total_precip)
    
    MSP = np.mean(SPs, axis = 0)
    
    return MSP, SPs
    
    
def mean_postemps_pts(year_range):
    
    MT = []
    
    #Total snow pack
    for i in year_range:
        
        temp = np.genfromtxt('Temptuning' + str(i) + '.txt', delimiter = ',')
        #temp[temp > 0] = 0
        mean_temp = np.mean(temp, axis = 0)
        mean_temp[mean_temp > 0] = 0
        MT.append(mean_temp)
    
    return MT
    


        ###melt model for glaciers with tuning###
def MB_simplified(Thour, Phour, SRhour, Leftover_in, asnow, aice, MF):
    

        Melt_list = []
        Leftover_list = []
        MBhour = []

        #ETIM melt model and snow tracker        
        for t in range(0, len(Thour)):
            if Thour[t] <= 0:
                Melt = 0.
                
                Melt_list.append(Melt)
                snow = Leftover_in[t]   
                Leftover_list.append(snow)
                               
            else:
                snow = Leftover_in[t]
                Melt_snow = (MF + asnow * SRhour[t]) * Thour[t]
                
                #EDD model using I arc    
                if Melt_snow < snow:
                    Melt = Melt_snow
                    Melt_list.append(Melt)
                    Leftover = snow - Melt_snow
                    Leftover_list.append(Leftover)
                    
                else:
                    Melt_snow = snow
                    Tice = Melt_snow/(MF + asnow * SRhour[t])
                    Melt_ice = (MF + aice * SRhour[t]) * (Thour[t] - Tice)
                    Melt = Melt_snow + Melt_ice
                    Melt_list.append(Melt)
                    Leftover = 0.
                    Leftover_list.append(Leftover)
            
            
            
        Mass_in = Phour
        Mass_out = Melt_list
        MB = np.asarray(Mass_in) - np.asarray(Mass_out)
        
        MBhour.append(MB) 
        
        return MBhour, Melt_list, Leftover_list
        
        
# def MB_vectorized(Thour, Phour, SRhour, Leftover_in, asnow, aice, MF, CC_in):   
# 
# 
#     Melt_list = np.empty(Thour.shape)
#     Leftover_list = np.empty(Thour.shape)
#     MBhour = np.empty(Thour.shape)
#     nan_locs = np.isnan(Thour)
#     CC_out = np.empty(Thour.shape)
#     #ETIM melt model and snow tracker  
#     
#     #get indices of ice and snowmelt based on snow melt delta and T
#     Melt_snow = (MF + asnow * SRhour) * Thour
#     DDsnow = (Leftover_in - CC_in)/(MF + asnow * SRhour)
#     ##############################################################################################################CHECK THIS
#     snowmelt_ind = np.where(np.logical_and(Thour > 0., (Melt_snow - CC_in) < Leftover_in, Melt_snow > CC_in))
#     icemelt_ind = np.where(np.logical_and(Thour > 0., (Melt_snow - CC_in) >= Leftover_in, Melt_snow > CC_in))
#     
#     #update CC lists
#     CC_0_ind = np.where(np.logical_and(Thour > 0., Melt_snow > CC_in))
#     CC_ind = np.where(np.logical_and(Thour > 0., Melt_snow < CC_in))
#     
#     CC_out[CC_0_ind] = np.zeros(Thour.shape)[CC_0_ind] 
#     CC_out[CC_ind] = CC_in[CC_ind] - Melt_snow[CC_ind]
#     #update melt list where melt is all CC reduction
#     Melt_list[CC_ind] = np.zeros(Thour.shape)[CC_ind]
#     
#     #update snow melt arrays
#     Leftover_list[snowmelt_ind] = Leftover_in[snowmelt_ind] - Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind]
#     Melt_list[snowmelt_ind] = Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind]
#     #MBhour[snowmelt_ind] = np.zeros(Thour.shape)[snowmelt_ind]  
#     
#     #get indices for icemelt
#     #snow melt is Melt_snow[icemelt_ind]
#     DDsnow = (Leftover_in - CC_in)/(MF + asnow * SRhour)
#     Melt_ice = (MF + aice * SRhour) * (Thour - DDsnow)
#     Melt_list[icemelt_ind] = Melt_ice[icemelt_ind] + Leftover_in[icemelt_ind] - CC_in[icemelt_ind]
#     Leftover_list[icemelt_ind] = np.zeros(Thour.shape)[icemelt_ind]
#     
#     #get indices of no melt, and update outputs
#     nomelt_ind = np.where(Thour <= 0.) 
#     Leftover_list[nomelt_ind] = Leftover_in[nomelt_ind]
#     Melt_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
#     CC_out[nomelt_ind] = CC_in[nomelt_ind]
#     #MBhour[noimelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
#            
#     #calculate mass balance           
#     #Mass_in = Phour
#     #Mass_out = Melt_list
#     MBhour = Phour - Melt_list
#     
#     MBhour[nan_locs] = np.nan
#     Melt_list[nan_locs] = np.nan
#     Leftover_list[nan_locs] = np.nan
#     CC_out[nan_locs] = np.nan
#     
#     
#     
#     return MBhour, Melt_list, Leftover_list, CC_out
    
    
def MB_vectorized(Thour, Phour, SRhour, Leftover_in, asnow, aice, MF, CC_in):  
    
    #create empty arrays for outputs
    Melt_list = np.empty(Thour.shape)
    Leftover_list = np.empty(Thour.shape)
    MBhour = np.empty(Thour.shape)
    nan_locs = np.isnan(Thour)
    CC_out = np.empty(Thour.shape)
   
    #Get Snow melt/icemelt/etc
    Melt_snow = (MF + asnow * SRhour) * Thour
    DDsnow = (Leftover_in + CC_in)/(MF + asnow * SRhour)
    Melt_ice = (MF + aice * SRhour) * (Thour - DDsnow)
   
    #get indices of ice and snowmelt based on snow melt delta and T
    snowmelt_ind = np.where(np.logical_and(Thour > 0., (Melt_snow - CC_in) < Leftover_in, Melt_snow > CC_in))
    icemelt_ind = np.where(np.logical_and(Thour > 0., (Melt_snow - CC_in) >= Leftover_in, Melt_snow > CC_in))
    #update CC lists
    CC_0_ind = np.where(np.logical_and(Thour > 0., Melt_snow > CC_in))
    CC_ind = np.where(np.logical_and(Thour > 0., Melt_snow < CC_in))
   
    #Update CC lists
    CC_out[CC_0_ind] = np.zeros(Thour.shape)[CC_0_ind]
    CC_out[CC_ind] = CC_in[CC_ind] - Melt_snow[CC_ind] 
   
    #update snow melt arrays
    Leftover_list[snowmelt_ind] = Leftover_in[snowmelt_ind] - (Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind])
    Melt_list[snowmelt_ind] = Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind] 
   
    #get indices for icemelt
    #snow melt is Melt_snow[icemelt_ind]
    Melt_list[icemelt_ind] = Melt_ice[icemelt_ind] + (Leftover_in[icemelt_ind] - CC_in[icemelt_ind])
    Leftover_list[icemelt_ind] = np.zeros(Thour.shape)[icemelt_ind]
    
    #get indices of no melt, and update outputs
    nomelt_ind = np.where(Thour <= 0.)
    Leftover_list[nomelt_ind] = Leftover_in[nomelt_ind]
    Melt_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
    CC_out[nomelt_ind] = CC_in[nomelt_ind]
    #update melt list where melt is all CC reduction
    Melt_list[CC_ind] = np.zeros(Thour.shape)[CC_ind]
          
    #calculate mass balance          
    #Mass_in = Phour
    #Mass_out = Melt_list
    MBhour = Phour - Melt_list
   
    #Fill in with nans
    MBhour[nan_locs] = np.nan
    Melt_list[nan_locs] = np.nan
    Leftover_list[nan_locs] = np.nan
    CC_out[nan_locs] = np.nan
   
   
   
    return MBhour, Melt_list, Leftover_list, CC_out 
    


    
def MB_vectorized_discreteSnI(Thour, Phour, SRhour, Leftover_in, asnow, aice, MF, Topo, CC_in, meltfactors, debtreatment):   


    Melt_list = np.zeros(Thour.shape) #will track NET ablation (total melt - refreezing)
    #Total_ablation_list = np.empty(Thour.shape) #adding a new list to track TOTAL MELT
    Snowmelt_list = np.zeros(Thour.shape)
    Icemelt_list = np.zeros(Thour.shape)
    Leftover_list = np.zeros(Thour.shape)
    MBhour = np.zeros(Thour.shape)
    nan_locs = np.isnan(Thour)
    CC_out = np.zeros(Thour.shape)

    #ETIM melt model and snow tracker  
    
    #get indices of ice and snowmelt based on snow melt delta and T
    Melt_snow = (MF + asnow * SRhour) * Thour
    DDsnow = (Leftover_in + CC_in)/(MF + asnow * SRhour)
    if debtreatment == 'Boolean':
        Melt_ice = ((MF + aice * SRhour) * (Thour - DDsnow)) #aice should already be zero in debris covered cells
    else:
        Melt_ice = ((MF + aice * SRhour) * (Thour - DDsnow)) * meltfactors
    #Melt_ice_totalablation = (MF + aice * SRhour) * (Thour)
    
    #get indices of ice and snowmelt based on snow melt delta and T
    snowmelt_ind = np.where(np.logical_and(Thour > 0., (Melt_snow - CC_in) < Leftover_in, Melt_snow > CC_in))
    icemelt_ind = np.where(np.logical_and(Thour > 0., (Melt_snow - CC_in) >= Leftover_in, Melt_snow > CC_in))
    
    # snow and icemelt indices which neglect refreezing!! --> use to calculate total ablation
    #snowmelt_ind_totalablation = np.where(np.logical_and(Thour > 0., (Melt_snow) < Leftover_in, Melt_snow >= 0))
    #icemelt_ind_totalablation = np.where(np.logical_and(Thour > 0., (Melt_snow) >= Leftover_in, Melt_snow >= 0))
    
    #update CC lists
    CC_0_ind = np.where(np.logical_and(Thour > 0., Melt_snow > CC_in))
    CC_ind = np.where(np.logical_and(Thour > 0., Melt_snow < CC_in))
    CC_out[CC_0_ind] = np.zeros(Thour.shape)[CC_0_ind]
    CC_out[CC_ind] = CC_in[CC_ind] - Melt_snow[CC_ind]      
    
    #update snow melt arrays
    Leftover_list[snowmelt_ind] = Leftover_in[snowmelt_ind] - (Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind])
    Melt_list[snowmelt_ind] = Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind] 
    Snowmelt_list[snowmelt_ind] = Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind]
    #Total_ablation_list[snowmelt_ind_totalablation] = Melt_snow[snowmelt_ind_totalablation]
    
    #update ice melt arrays
    Melt_list[icemelt_ind] = Melt_ice[icemelt_ind] + (Leftover_in[icemelt_ind] - CC_in[icemelt_ind])
    Icemelt_list[icemelt_ind] = Melt_ice[icemelt_ind] + (Leftover_in[icemelt_ind] - CC_in[icemelt_ind])
    Leftover_list[icemelt_ind] = np.zeros(Thour.shape)[icemelt_ind]
    
    #Total_ablation_list[icemelt_ind_totalablation] = Melt_ice_totalablation[icemelt_ind_totalablation] + Leftover_in[icemelt_ind_totalablation]
    
    #get indices of no melt, and update outputs
    nomelt_ind = np.where(Thour <= 0.)
    no_icemelt_ind = np.where(Topo == 0)
    Leftover_list[nomelt_ind] = Leftover_in[nomelt_ind]
    Melt_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
    Snowmelt_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
    Icemelt_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
    Icemelt_list[no_icemelt_ind] = 0
    
    CC_out[nomelt_ind] = CC_in[nomelt_ind]
    #update melt list where melt is all CC reduction
    Melt_list[CC_ind] = np.zeros(Thour.shape)[CC_ind]
    Snowmelt_list[CC_ind] = np.zeros(Thour.shape)[CC_ind]
    Icemelt_list[CC_ind] = np.zeros(Thour.shape)[CC_ind]
    
    #Total_ablation_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]

    #Update both ice and snow outputs
    #add snowmelt for snowmelting
    #Snowmelt_list[snowmelt_ind] = Melt_snow[snowmelt_ind] - CC_in[snowmelt_ind]
    #and icemelting where snowpack is consumed
    #Snowmelt_list[icemelt_ind] = Leftover_in[icemelt_ind]
    #add icemelt
    #Icemelt_list[icemelt_ind] = Melt_ice[icemelt_ind] + (Leftover_in[icemelt_ind] - CC_in[icemelt_ind])
    #Snowmelt_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
    #Icemelt_list[nomelt_ind] = np.zeros(Thour.shape)[nomelt_ind]
    #setup CC indices in both arrays
    #Snowmelt_list[CC_ind] = np.zeros(Thour.shape)[CC_ind]
    #Icemelt_list[CC_ind] = np.zeros(Thour.shape)[CC_ind]
    #and set ice melt to 0 at office gridecells: assume off ice = 0 and onice = 1
    Icemelt_list = Icemelt_list * Topo

           
    #calculate mass balance           
    #Mass_in = Phour
    #Mass_out = Melt_list
    MBhour = Phour - Melt_list
    
    MBhour[nan_locs] = np.nan
    Melt_list[nan_locs] = np.nan
    Leftover_list[nan_locs] = np.nan
    Icemelt_list[nan_locs] = np.nan
    Snowmelt_list[nan_locs] = np.nan
    CC_out[nan_locs] = np.nan
    
    #Total_ablation_list[nan_locs] = np.nan
    
    
    return MBhour, Melt_list, Leftover_list, Icemelt_list, Snowmelt_list, CC_out
    
    
    
# Polynomial Regression
#https://stackoverflow.com/questions/893657/how-do-i-calculate-r-squared-using-python-and-numpy
def polyfit_homebrew(x, y, degree):
    results = []

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    #results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    #results['determination'] = ssreg / sstot
    results.append(ssreg / sstot)

    return results
    
    
#function to pull mass balance curve from distributed mass balance field
#inputs are distributed net mass balance for a given time period, grid cell 
#elevations, and ELA for cutoff
def Check_out_that_MB(flat_netMB, flat_z, ELA):
    
    Z, MASSB = zip(*sorted(zip(flat_z, flat_netMB)))

    zbins = np.linspace(0, 4000, 41)

    MBvals = [] 
    
    for i in range(0, len(zbins)):
        
        MBvals_bin = []
            
        for j in range(0, len(flat_netMB)):
            
            if zbins[i] < Z[j] < zbins[i+1]:
                MBvals_bin.append(flat_netMB[j])
            else:
                pass
                
        MBvals.append(MBvals_bin)
        
    means = []
    stds = []
    medians = []

    for k in range(0, len(MBvals)):
        means.append(np.nanmean(MBvals[k]))
        stds.append(np.nanstd(MBvals[k]))
        medians.append(np.nanmedian(MBvals[k]))
        
    return means, stds, medians

def KWtributaries():
    """
    balance flux codes:
    0 = KW0
    1 = NA
    2 = CA
    3 = SW
    4 = SA
    5 = KW5
    6 = KW4
    7 = KW3
    8 = KW2
    9 = KW1
    tributary codes:
    SA = 1
    SW = 2
    CA = 3
    NA = 4
    Trunk = 5
    """
    Path2KWoutline = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel'
    File_glacier_in = os.path.join(Path2KWoutline,'kaskonly_deb.txt')
    glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
    Ix = glacier[:,3] 
    Iy = glacier[:,4] 
    balancefluxes = glacier[:,5]

    #-------Turn vectors into 2D gridded inputs--------------------

    BalFluxes, Xgrid_KW, Ygrid_KW, xbounds, ybounds = regridXY_something(Ix, Iy, balancefluxes)
    nanlocs = np.where(np.isnan(BalFluxes))
    tribarray = np.zeros(BalFluxes.shape)
    tribarray[nanlocs] = np.nan
    
    tribarray[np.where(BalFluxes >= 5)] = 5 # Trunk
    tribarray[np.where(BalFluxes == 0)] = 5 # Trunk
    tribarray[np.where(BalFluxes == 1)] = 4 # NA
    tribarray[np.where(BalFluxes == 2)] = 3 # CA
    tribarray[np.where(BalFluxes == 3)] = 2 # SW
    tribarray[np.where(BalFluxes == 4)] = 1 # SA
    
    return tribarray

    
def generate_meltfactors(debrismask,cleaniceM,peakM,peakM_thickness,transition_thickness,b0 = 11.0349260206858,k = 1.98717418666925):
    '''
    the final version of this function should always be in the Model Functions script
    #function should take in:
        # debris mask,
        # cleanicemelt and peak melt value
        # transition thickness and thickness at peak melt 
        # Rounce et al. model params for the kaskawulsh
    #function should return:
        # 2D melt factors map with same shape as input debris mask. should be 1 in non debris areas
    so that it doesnt change the melt of cleanice when multiplied with the melt array
    '''
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    #make an empty array of same shape as debris mask
    #loop through each gridcell in debris mask
    #if nan (ie no debris) place a 1 in the meltfactors array
    #sort the debris cells by thickness and calculate corresponding MF, then place in melt array
    
    meltfactors = np.empty(debrismask.shape)
    for x in range(len(debrismask)):
        for y in range(len(debrismask[0])):
            #print(debrismask[x,y])
            if np.isnan(debrismask[x,y]):
                meltfactors[x,y] = 1
            elif debrismask[x,y] <= peakM_thickness:
                Melt = (((peakM - cleaniceM)/peakM_thickness)*debrismask[x,y]) + cleaniceM
                Meltfact = Melt/cleaniceM
                meltfactors[x,y] = Meltfact
                #print(str(debrismask[x,y]) + ' ,,, ' + str(Meltfact))
                #calculate the melt factor on a linear line between clean ice melt and peak melt:
            elif (debrismask[x,y] > peakM_thickness) and (debrismask[x,y] <= transition_thickness):
                #print(debrismask[x,y])
                yintercept = peakM - (((cleaniceM-peakM)/(transition_thickness-peakM_thickness))*peakM_thickness)
                Melt = (((cleaniceM-peakM)/(transition_thickness-peakM_thickness))*debrismask[x,y]) + yintercept
                Meltfact = Melt/cleaniceM
                meltfactors[x,y] = Meltfact
            elif debrismask[x,y] > transition_thickness:
                #determine how much this part of the curve needs to be shifted over to meet the linear parts
                h_curvestart = (b0-cleaniceM)/(cleaniceM*k*b0)
                horizontal_shift = transition_thickness - h_curvestart
                Melt = b0/(1+(k*b0*(debrismask[x,y]-horizontal_shift))) #calculate the melt on the original curve(take away the shift) so that you're actually getting the melt at the desired thickness
                Meltfact = Melt/cleaniceM
                meltfactors[x,y] = Meltfact
            else:
                print('debris thickness not accounted for in if statements')

    return meltfactors 

def write_config_file(outputpath,scriptname):
    with open(scriptname) as f:
        data = f.read()
        f.close()

    txtfile = os.path.join(outputpath,"namelist.txt")    
    with open(txtfile, mode="w") as f:
        f.write(data)
        f.close()
        
def model_domain(catchment):
    if catchment == True:
        File_glacier_in = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kask_catchment.txt'
        glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
                
        Ix = glacier[:,4] 
        Iy = glacier[:,5] 
        Ih = glacier[:,6]  
        sfc_type = glacier[:,8] 
        Sfc, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, sfc_type)
    else:
        File_glacier_in_KW = 'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel\kaskonly_deb.txt'
        glacier_KW = np.genfromtxt(File_glacier_in_KW, skip_header=1, delimiter=',')
                    
        Ix = glacier_KW[:,3] 
        Iy = glacier_KW[:,4] 
        Ih = glacier_KW[:,2]    
        Sfc = np.nan
        
    Zgrid, Xgrid, Ygrid, xbounds, ybounds = regridXY_something(Ix, Iy, Ih)
    
    return Zgrid, Xgrid, Ygrid, xbounds, ybounds, Sfc

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

