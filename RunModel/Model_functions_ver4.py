import numpy as np
from scipy.spatial import distance
from sklearn import linear_model
from scipy.stats import linregress
import math
from scipy import interpolate
from netCDF4 import Dataset
import netCDF4
import os
import pandas as pd


###Downscale temperature based on atmospheric structure and inversions###
def T_downscale_funkfest(T, E, UTMx_list, UTMy_list):

        
        ###Get values for geopotential and temperature
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
        for u in range(0, T[-1][:][:].shape[0]):
            for w in range(0, T[-1][:][:].shape[1]): 
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
          
     #use scipy.inter2d      
        # y0func = interpolate.interp2d(UTMx_list, UTMy_list, y0_list, kind='linear')
        # y0func_inver = interpolate.interp2d(UTMx_list, UTMy_list, y0_list_inver, kind='linear')
        # Lfunc = interpolate.interp2d(UTMx_list, UTMy_list, L_list, kind='linear')
        # Lfunc_inver = interpolate.interp2d(UTMx_list, UTMy_list, L_list_inver, kind='linear')
     #use scipy.Rbf   
        # y0func = interpolate.Rbf(UTMx_list, UTMy_list, y0_list, function='linear', smooth = 0)
        # y0func_inver = interpolate.Rbf(UTMx_list, UTMy_list, y0_list_inver, function='linear', smooth = 0)
        # Lfunc = interpolate.Rbf(UTMx_list, UTMy_list, L_list, function='linear', smooth = 0)
        # Lfunc_inver = interpolate.Rbf(UTMx_list, UTMy_list, L_list_inver, function='linear', smooth = 0)
     #use scipy.griddata
        y0func = interpolate.bisplrep(UTMx_list, UTMy_list, y0_list)
        y0func_inver = interpolate.bisplrep(UTMx_list, UTMy_list, y0_list_inver)
        Lfunc = interpolate.bisplrep(UTMx_list, UTMy_list, L_list)
        Lfunc_inver = interpolate.bisplrep(UTMx_list, UTMy_list, L_list_inver)
                          
                    
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
    norm_gridy = np.mean(grid_spacey).round(-2)
    
    if np.mean(grid_spacex) == norm_gridx:
        pass
    else:
        holes = np.where(grid_spacex > norm_gridx)
        d = grid_spacex[holes]
        vals = xr[holes]
        for val in range(0, len(vals)):
            v = vals[val]
            D = d[val]/norm_gridx
            for i in range(1, int(D)):
                new_val = v + (i * norm_gridx)
                xr = np.insert(xr, 0, new_val)
        
    if np.mean(grid_spacey) == norm_gridy:
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
    
    X, Y = np.meshgrid(x, y)
        
    hold = np.empty(X.shape)
    hold[:] = np.nan
        
    XYZ = np.empty(X.shape)    
    XYZ[:] = np.ones(X.shape) * hold
    
        
    for i in range(0, len(X)):
        #print 'i = '
        #print i
        for j in range(0, len(X[i])):
            #print 'j = '
            #print j
         
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
def closest_node(node, nodes):
    
    dist_2 = np.sum((nodes - node)**2, axis=0)
    node_pos = np.argmin(dist_2)
    #closest_index = distance.cdist([node], nodes).argmin()
    
    return node_pos
    

def rainy_day_funk(elev, precip, lat, lon):

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
def netcdf_container_gen(melt_model, var_n, File_name, File_sufix, ybounds, xbounds, This_year_is):

    ###Create empty .nc file for MB values###
    ref_file_in = 'ref' + str(This_year_is) + '.nc'    
    fh = Dataset(ref_file_in, "r")
    f = Dataset(File_name + File_sufix, 'w', format='NETCDF4') #'w' stands for write
        
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
    mb = f.createVariable(var_n, np.float32, ('time', 'y', 'x'))
        
    # Global Attributes 
    f.description = 'Model Output' + File_name   
    f.source = 'Melt model for Kaskawulsh' 
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
    mb[:] = np.zeros(mb[:].shape) * np.nan
    mb[:] = melt_model
    
    f.close()
    
    #return(mb[:], f)


def get_meanSP(year_list, Glacier_id,R2S,BiasCorrectionT,BiasCorrectionP,T_path,P_path):

    #mean annual total precipitation
    DH_list = []
    for y in year_list:
        
        if BiasCorrectionT == True:
            File_temp_in_future = 'Temp' + Glacier_id + '_BC_' + str(y) + '.nc'
        else:
            File_temp_in_future = 'Temp' + Glacier_id + str(y) + '.nc'
            
        if BiasCorrectionP == True:
            File_precip_in_future = 'Prcp' + Glacier_id + '_BC_' + str(y) + '.nc'
            P_var = 'Net snow'
        else:
            File_precip_in_future = 'netSnow' + Glacier_id + str(y) + '.nc'
            P_var = 'Temperature'

        File_t_in = os.path.join(T_path,File_temp_in_future)
        File_p_in = os.path.join(P_path,File_precip_in_future)
        
        infut = Dataset(File_p_in, "r")
        tinfut = Dataset(File_t_in, "r")
        T_var = 'Temperature'
        T_array = tinfut.variables[T_var][:]
        
        snowlocs = np.where(T_array > R2S)        
        
        #P_var = 'Net snow'
        P_array = infut.variables[P_var][:]
        P_array[snowlocs] = 0.
        
        DH = np.sum(P_array, axis = 0)
        DH_list.append(DH)
        
        infut.close()
        tinfut.close()
        
        
    DHval = np.mean(DH_list, axis = 0)     

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
        File_precip_in_future = 'Prcp' + Glacier_id + '_BC_' + str(This_year_is) + '.nc'
    else:
        File_precip_in_future = 'netSnow' + Glacier_id + str(This_year_is + 1) + '.nc'
    
    File_p_in = os.path.join(P_path,File_precip_in_future)
    
    infut = Dataset(File_p_in, "r")
    if BC_P == True:
        P_var = 'Net snow'
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
    


    
def MB_vectorized_discreteSnI(Thour, Phour, SRhour, Leftover_in, asnow, aice, MF, Topo, CC_in):   


    Melt_list = np.empty(Thour.shape) #will track NET ablation (total melt - refreezing)
    #Total_ablation_list = np.empty(Thour.shape) #adding a new list to track TOTAL MELT
    Snowmelt_list = np.empty(Thour.shape)
    Icemelt_list = np.empty(Thour.shape)
    Leftover_list = np.empty(Thour.shape)
    MBhour = np.empty(Thour.shape)
    nan_locs = np.isnan(Thour)
    CC_out = np.empty(Thour.shape)

    #ETIM melt model and snow tracker  
    
    #get indices of ice and snowmelt based on snow melt delta and T
    Melt_snow = (MF + asnow * SRhour) * Thour
    DDsnow = (Leftover_in + CC_in)/(MF + asnow * SRhour)
    Melt_ice = (MF + aice * SRhour) * (Thour - DDsnow)
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

def KWtributaries(coordinatefile,Zgrid):
    """
    codes:
    nontrib = 0
    SA = 1
    SW = 2
    CA = 3
    NA = 4
    Trunk = 5
    """
    #get nanlocs
    #nanlocs = np.where(np.isnan(Zgrid))
    
    #load .csv file with tributary codes
    df = pd.read_csv(coordinatefile)
    #df['Tributary'].astype('Int64')
    arr = np.array(df)
    
    tribs = arr[:,3]
    
    ones = np.where(tribs == 'SA')
    twos = np.where(tribs == 'SW')
    threes = np.where(tribs == 'CA')
    fours = np.where(tribs == 'NorthArm')
    fives = np.where(tribs == 'Tr') # trunk
    zeros = np.where(pd.isna(tribs))
    
    tribs[ones] = 1
    tribs[twos] = 2
    tribs[threes] = 3
    tribs[fours] = 4
    tribs[fives] = 5
    tribs[zeros] = 0
        
    #reshape tribcodes list in Zgrid.shape
    tribarray = np.reshape(tribs,(Zgrid.shape))
    
    tribarray_f = np.flipud(tribarray)
    #tribarray[nanlocs] = np.nan
    
    return tribarray_f

    
