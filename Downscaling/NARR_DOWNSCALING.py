        ###Mass Balance model adapeted to netcdf files, current use###
        ###Ice-marginal lake drainage basin (4x ice sfc and 1x topo)###
        
######IMPORTANT################
#                             #
# RUN DOWNSCALINGnamelist.py  # 
#    BEFORE running this      #
#          script!            #
#                             #
###############################
        
import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate
import netCDF4
from pyproj import Proj
from netCDF4 import Dataset
import datetime as dt
from Model_functions_ver3 import T_downscale_funkfest
from Model_functions_ver3 import Precip_2_Snow
from Model_functions_ver3 import MB_341_debris
from Model_functions_ver3 import MB_341_office
from Model_functions_ver3 import rainy_day_funk
from Model_functions_ver3 import closest_node
import sys
import os 

#Import parameters from config file
#from DOWNSCALINGnamelist import *
from DOWNSCALINGnamelist import start_year
from DOWNSCALINGnamelist import end_year
from DOWNSCALINGnamelist import glacier_id
from DOWNSCALINGnamelist import Topo_param
from DOWNSCALINGnamelist import Debris
from DOWNSCALINGnamelist import UTM
#from DOWNSCALINGnamelist import ELA
from DOWNSCALINGnamelist import SNOW_START
from DOWNSCALINGnamelist import gen_snow
from DOWNSCALINGnamelist import NARR_subregions
from DOWNSCALINGnamelist import D_T
from DOWNSCALINGnamelist import D_P
from DOWNSCALINGnamelist import delta_t
from DOWNSCALINGnamelist import cfactor
from DOWNSCALINGnamelist import snowfactor
from DOWNSCALINGnamelist import downscaled_outputs
from DOWNSCALINGnamelist import glacier_outline
from DOWNSCALINGnamelist import considering_kaskonly
from DOWNSCALINGnamelist import rawnarr_inputs
from DOWNSCALINGnamelist import solar_in

#set path for downscaled outputs(Temp, Precip, & Radiation) to go into 
OUTPUT_PATH = downscaled_outputs


## set up time range ###
years = []
year = start_year
while year <= end_year:
  years.append(year)
  year = year + 1
#print(years)

for year in years:
                 #######INPUT Parameters#######


    print ("Glacier ID is: " + glacier_outline)
    sys.stdout.flush()
    
    #Boolean for whether on ice or off ice
    Topo = Topo_param

    if Topo == True:
        print ("Off ice grid,")
    else:
        print ("On ice grid,")
    
    #Boolean for whether to include debris cover map and no melt condition
    Debris_cover = Debris
    
    if Debris_cover == False:
        print ("No debris cover considered,")
    else:
        print ("Debris cover considered")
        
    #UTM zone used for reprojection of downscaled grid
    UTM_zone = UTM
    
    #ELA value used for area
    #ela = ELA
    #print(ela)
    
    #Boolean for whether initial snow profile exists or start with no snow condition
    #Followed by boolean for whether to generate snow.txt output
    snow_start = SNOW_START
    Generate_snow = gen_snow
    
    if snow_start == True:
        print ("Snow initial condition is 0 m.w.e.,")
    else:
        print ("Snow initial condition is carried over from previous year") 
    
    #Indices of subregions used for precip downscaling, picked manually to omit 
    #points on opposite side of divide
    NARR_indices = NARR_subregions
    
    #Set year for simulation run
    current_year = year
    
    print('it is ' +str(current_year))
    
    #Implement bias correction; given array are DT values developped indepedently of
    #model for time period, interpolated linearly between months at a daily 
    #resolution 
    DT = D_T
    DTval = interp1d(np.asarray([1.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305. \
        ,335., 366.]), np.asarray([-4.86,-3.98,-2.1,0.03,0.89,1.42,0.94,0.32,-0.32,-1.72,-4.37,-5.22, -4.86]), kind = 'linear')   
    
        
    #ranges for DP correction. More for reference than actual algorithm
    DP = D_P
    DPval = interp1d(np.asarray([500, 900, 1100, 1300, 1500, 1700, 2100, 2300, 2500, \
        2700, 3100, 3700]), np.asarray([ 1.28530635,  1.26180458,  1.23360247,  1.25240269,  1.28293036,
        1.26320929,  1.24449735,  1.30082501,  2.46677552,  4.24804321,
        5.45333788,  5.45333788]), kind = 'linear')
        
                        #######INPUTS#######
    
    ###Load data from NETcdf files, see NARR documentation
    ###for further info, of netcdf_spatial.py for concatenation
    ###steps performed on original datasets, use cdo in terminal
    ###for easy navigation
    
    
    File_elev_in = os.path.join(rawnarr_inputs,'kaskhgt.' + str(current_year) + '.nc')
    File_temp_in = os.path.join(rawnarr_inputs,'kaskair.' + str(current_year) + '.nc')
    File_precip_in = os.path.join(rawnarr_inputs,'kaskapcp.' + str(current_year) + '.nc')
    File_CoarseElev_in = 'kaskCE.nc'
    File_glacier_in = glacier_outline + '.txt' 
    if Debris_cover == True:
        File_debris_in =  glacier_outline + '.txt'
    else:
        pass
    if snow_start == False:
        File_snow_in = 'kaskonlysnow' + str(current_year) + '.txt'
    else:
        pass

    
    inE = Dataset(File_elev_in, "r")
    E_var = 'hgt'
    print("Geopotential in...")
    sys.stdout.flush()
    
    inT = Dataset(File_temp_in, "r")
    T_var = 'air'
    print("Temperature in...")
    sys.stdout.flush()
    
    inP = Dataset(File_precip_in, "r")
    P_var = 'apcp'
    print("Precipitation in...")
    sys.stdout.flush()
    
    inCE = Dataset(File_CoarseElev_in, "r")
    CE_var = 'hgt'
    print("Coarse DEM in...")
    sys.stdout.flush()
    
    #Solar data prefixes for current run
    if considering_kaskonly == True:
        solar_prefix = 'fixed' + glacier_id + '_' #GLright/left/mid, Kaskice, kasktopo
        solar_suffix = 'DS.txt'
    else:
        solar_prefix = 'fixed' + 'kask' + '_' #GLright/left/mid, Kaskice, kasktopo
        solar_suffix = 'DS.txt'
    
    #Get values arrays from binary files
    E_array = inE.variables[E_var][:]
    T_array = inT.variables[T_var][:]
    P_array = inP.variables[P_var][:]
    #get coarse NARR grid
    CE_array = inCE.variables[CE_var]
    elev = CE_array[0,:,:]
    
    #grbindexGR = pygrib.index(File_GR_in, 'date')
    print("Netcdf values loaded....")
    sys.stdout.flush()
    
    #import glacier grid(s)
    glacier = np.genfromtxt(File_glacier_in, skip_header=1, delimiter=',')
    if snow_start == False:
        snow_ini = np.genfromtxt(File_snow_in, delimiter=',')
    else:
        pass
    
    if Topo == 0:
        if Debris_cover == True:
            debris_mask = np.genfromtxt(File_debris_in, skip_header=1, delimiter=',')
        else:
            pass
    else:
        pass
    
    
    
                        #######OUTPUTS#######
    
    ###Prep output file and variables
    
    #Prep non-averaged time step outputs
    File_out1 = open( 'Melt' + glacier_id + str(current_year) + '_dd.txt', 'w')
    File_out2 = open( 'Melt' + glacier_id + str(current_year) + '_edd.txt', 'w')
    
    #Prep non-averaged time step outputs
    File_out4 = open( 'MBG' + glacier_id + str(current_year) + '_dd.txt', 'w')
    File_out5 = open( 'MBG' + glacier_id + str(current_year) + '_edd.txt', 'w')
    
    #Outputs of downscaled params
    #TEMPERATURE:
    Temp_outfile_name = 'Temp' + glacier_id + str(current_year) + '.txt'
    Temp_output_path = os.path.join(OUTPUT_PATH,Temp_outfile_name)
    File_out3 = open(Temp_output_path, 'w')
    
    #SOLAR RADIATION:
    Srad_outfile_name = 'Srad' + glacier_id + str(current_year) + '.txt'
    Srad_output_path = os.path.join(OUTPUT_PATH,Srad_outfile_name)
    File_out7 = open(Srad_output_path, 'w')
    
    #ACCUMULATION:
    netSnow_outfile_name = 'netSnow' + glacier_id + str(current_year) + '.txt'
    netSnow_output_path = os.path.join(OUTPUT_PATH,netSnow_outfile_name)
    File_out11 = open(netSnow_output_path, 'w')
    
    File_out8 = open('AccSnow' + glacier_id + str(current_year) + '_dd.txt', 'w')
    File_out9 = open('AccSnow' + glacier_id + str(current_year) + '_edd.txt', 'w')
    File_out10 = open('precip' + glacier_id + str(current_year) + '.txt', 'w')
    
                #######SETUP Model domain#######
    
    ###Prep time and space resolution for model, and arrays used for
    ###Operations during MB calculations
    
    #set up time steps in a day...
    time_delta = delta_t
    t_st = 0
    t_nd = 24
    timeinterv = np.linspace(t_st, t_nd, t_nd/time_delta, endpoint=False)
    
    #...and set up timesteps through time, with days limited array for indexing
        #perform twice for hourly values and daily precip
    units = inE.variables['time'].units
    dates = netCDF4.num2date(inE.variables['time'][:], units)
    dates_p = netCDF4.num2date(inP.variables['time'][:], units)
    date_list = []
    date_list_p = []
    for d in dates:
        date_el = d.replace(hour=0, minute=0, second=0, microsecond=0)
        date_list.append(date_el)
    for d in dates_p:
        date_el = d.replace(hour=0, minute=0, second=0, microsecond=0)
        date_list_p.append(date_el)
    date_list = np.array(date_list)
    date_list_p = np.array(date_list_p)
    
    #get pressure level dimensions for pl variables
    plinterv = inE.variables['level'][:]
    
    #Get interv for solar data, and set t = -1 day solar values for PCT-UTC gap
    solhour_interv = 0.5
    ini_DS = np.zeros(3)
    
    #Extract grid geometry, correcting for dem/shapefile non-overlaps
    if considering_kaskonly == True:
        IH = glacier[:,6]
        #Ih[Ih==-9999]=np.nan
        inval_loc = np.where(IH == -9999)
        Ih =np.delete(IH, inval_loc)
        Ix = glacier[:,4]
        Iy = glacier[:,5]
        TOPO = glacier[:,8]
        ELA = glacier[:,10] #ELA is just an array with same shape as Zgrid, where all non-NaN cells are 0. 
    else:
        IH = glacier[:,2]
        #Ih[Ih==-9999]=np.nan
        inval_loc = np.where(IH == -9999)
        Ih = np.delete(IH, inval_loc)
        Ix = glacier[:,3]
        Iy = glacier[:,4]
        #TOPO = glacier[:,8] #TOPO is just an array with same shape as Ih, all vals are 0 
        #ELA = glacier[:,10] #ELA is also an array with same shape as Ih, all vals are 0 
        TOPO = np.zeros(Ih.shape)
        ELA = np.zeros(Ih.shape)
        
    
    #Create ELA array
    ELAs = []
    for j in range(0, len(ELA)):
        #0 means ice
        if TOPO[j] == 0:
            #when greater than ELA, 1 means add snow
            if Ih[j] > ELA[j]:
                ELAs.append(int(1))
            else:
                ELAs.append(int(0))
        else:
            ELAs.append(int(1))
    
        
    #setup factors and params for MB calculations
    cfac = cfactor
    snowfac = snowfactor
    MFi = (0.0062/len(timeinterv)) * cfac 
    MFs = (0.0054/len(timeinterv)) * cfac
    MF = (0.0026/len(timeinterv)) * cfac
    asnow = 0.00000056 * time_delta * cfac 
    aice = 0.00000086 * time_delta * cfac
    
    #arrays for dumping carry-over values for snow in  each model, set initial snow 
    #condition here: use previous years sow profile
    
    if year == years[0]:
        snow_ini_a = np.zeros(Ih.shape) # set snow profile intitial condition to zero for first year
    else:
        if snow_start == True:
            snow_ini_a = np.zeros(Ih.shape)
        else:    
            snow_ini_a = snow_ini
    
    Leftover_list = np.zeros(np.size(Ih)) + snow_ini_a
    Leftover_listSR_arc = np.zeros(np.size(Ih)) + snow_ini_a
    
    #Control for iffy debris cover input, courtesy of esri
    if Topo == 0:
        if Debris_cover == 1:
            debris_mask[:,1], debris_mask[:,5] = zip(*sorted(zip(debris_mask[:,1], debris_mask[:,5])))
            debris_raw = debris_mask[:,5]
            debris_m = np.zeros(debris_raw.shape) + debris_raw
            debris_m[np.where(debris_raw > 100)] = 0.
            debris_m[np.where(debris_raw <= 100)] = 1.
        else:
            debris_m = np.zeros(Ih.shape) + 1.
    else:
        pass
        
    #number of counts per day for inputs (might be irrelevant in netcdf design)
    #dayT = 32
    #dayP = 2
    
    #set start day for model run
    st_day = 1
    it_day = 0
    
    #Get UTM values for each reanalysis grid point to cycle between nodes on larger 
    #modelled glaciers
    lons = inE.variables['lon'][:]
    lats = inE.variables['lat'][:]
    
    myProj = Proj('+proj=utm +zone=' + UTM_zone + ', +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
    UTMy, UTMx = myProj(lons, lats)
    UTMx[UTMx >100000000000]=np.nan
    UTMy[UTMy >100000000000]=np.nan
    #create list of array positions ungridded
    UTMx_list = UTMx.ravel()
    UTMy_list = UTMy.ravel()
    grid_pts = np.stack((UTMx_list, UTMy_list))
    
    ###Setup starting positions for each input
    fruit = dates[(st_day*8) - 8]
    print("Start year:") 
    print(fruit.year)
    sys.stdout.flush()
    
    
    


                      #######MODEL#######

    ###Begin MB calculations, where dates are used to toggle while loop which
    ###iterates through every day in the record, then through every hour in the
    ###record, in order to obtain both daily and hourly (as per timestep) MB
    ###values across the study area. 
    
    
    while (fruit < dates[-1]):
        
        
        ###Get values for current iteration  
                
        #Current position in time of interest          
        cur_day = fruit.day
        cur_month = fruit.month
        cur_year = fruit.year
        today = np.where(date_list == np.datetime64(dt.datetime(cur_year, cur_month, cur_day)))
        today_p = np.where(date_list_p == np.datetime64(dt.datetime(cur_year, cur_month, cur_day)))
    
        #select values for this day 
        selectedE = E_array[today][:][:][:]
        selectedT = T_array[today][:][:][:]
        selectedP = P_array[today_p][:][:][:]
    
        #get daily rad values from arc files
        today_num = it_day + st_day    
        #UNCOMMENT
        curSR_arc = np.genfromtxt(os.path.join(solar_in,solar_prefix + str(today_num) + solar_suffix), delimiter=' ')
        
        
        for i in range(0, len(timeinterv)):
            
            ###Get values for current hour, indexing is yucky, but straighforward
            ###Note that units for NARR precip are mm we, conversion is done here  
            
            curE = selectedE[i][:][:][:]        
            curT = selectedT[i][:][:][:]
            #units conversion
            curP = selectedP[0][:][:]/1000
            
            #all sorts of empty arrays for calculations
            Thour = []
            Phour = []
            SRhour_arc = []
            SRhour_re = []
            Dhour = []
            precipList = []
            
      
            #Apply funkfest funktion to get inversion tracking downscaled T 
            #functions for every reanalysis grid point
            
            xi_list, yi_list, xi_list_inver, yi_list_inver, funclist, \
                funclist_inver, inversion_list, y0func, y0func_inver, Lfunc, \
                    Lfunc_inver = \
                        T_downscale_funkfest(curT, curE, UTMx_list, UTMy_list) 
                
            #Apply precip funk to get multivariate based downscaling of precip values
            if i == 0:
                W1, Pfunk, I0 = rainy_day_funk(elev.ravel()[NARR_indices], \
                    curP.ravel()[NARR_indices], UTMx_list[NARR_indices], \
                        UTMy_list[NARR_indices]) 
            else:
                pass           
            
            ###For large Kask type grid use if loop to assign reanalysis grid position                             
            #Get T value for every point on glacier
            #use appropriate interp functon and extrap values
            for z in range(0, len(Ih)):
                
                #Get closest NARR grid point for appropriate downscaling T values
                grid_pt = np.asarray(([Iy[z]], [Ix[z]]))
                #print(grid_pt)
                grid_loc = closest_node(grid_pt, grid_pts)
                #print(grid_loc)
            
                #use index to get nearest grid point in u, w notation
                u = int(np.where(UTMx == grid_pts[0][grid_loc])[0])
                w = int(np.where(UTMy == grid_pts[1][grid_loc])[1])
                
                #find downscaled value at every grid point. If outside bounds, 
                #extrapolate, if within bounds... New functionuses interpolated
                #Lapse rate and y0, meaning no extrapolation routine is needed 
                
                #implement bias correction
                if DT == True:
                    dToday = DTval(today_num)
                else:
                    dToday = 0.
                
                                            
                if inversion_list[u][w] == 0:
                    k = interpolate.bisplev(Iy[z], Ix[z], Lfunc) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func)                         
                        
                else:
                    
                    if Ih[z] < yi_list[u][w][0]:
                        k = interpolate.bisplev(Iy[z], Ix[z], Lfunc_inver) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func_inver) 
                    else:
                        k = interpolate.bisplev(Iy[z], Ix[z], Lfunc) * Ih[z] + interpolate.bisplev(Iy[z], Ix[z], y0func)  
                        
                K = (k - 273.3) + dToday  
                
                #Snag date for output for lapse rate, not currently implemented         
                D = timeinterv[i]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                    
                
                #get precip values at start of day
                if i == 0:
                    P1 = (I0 + (Pfunk[0] * Ix[z]) + (Pfunk[1] * Iy[z]) + (Pfunk[2] \
                        * Ix[z] * Iy[z]) + (Pfunk[3] * Ix[z]**2) + (Pfunk[4] * \
                            Iy[z]**2) + (Pfunk[5] * Ih[z]))*W1 
                    
                    #account for negative precip values in areas where statistical
                    #model predicts less than zero value on mx + b regression curve
                    if P1 < 0:
                        P1 = 0.
                    else:
                        pass
                    
                    P2 = curP[u][w]*(1-W1)
                    Pdval = P1 + P2
                    ##vals: Z = Ih[z], Y = Iy[z], X = Ix[z]
                    S = Precip_2_Snow(Pdval, K, snowfac)
                else:
                    S = 0.
                    Pdval = 0.
                
                #if statement to calculate the Precip bias correction factor:
                if DP == True:
                    deltaC = DPval(Ih[z])
                else:
                    deltaC = 1
                
                #Calculate SR/GR enhancement for both reanalysis dataset
                #add sum function when using hourly arc SR. NOTE: code exists for
                #averaged SR scheme in previous verison of this model
                #UNCOMMENT
                SR_arc = curSR_arc[z , :]
                #UNCOMMENT
                SR = SR_arc[int(timeinterv[i]/solhour_interv)]
                
                #reset late solar values
                #UNCOMMENT
                ini_DS = SR_arc[-4:-1]
                
                
                #Update lists  
                Thour.append(K)
                Dhour.append(D)
                Phour.append(S * deltaC)
                #UNCOMMENT
                SRhour_arc.append(SR)
                precipList.append(Pdval)
                
            
            ###Calculate and append hourly MB
            
            #get snow cover from previous timestep, but preserve snow above ELA        
            Pmelt = np.asarray(Phour) + np.asarray(Leftover_list) + np.asarray(ELAs)
            Pmelt[np.where(np.asarray(ELAs) == 1.)] = 1.
            Pmelt_SRarc = np.asarray(Phour) + np.asarray(Leftover_listSR_arc) + np.asarray(ELAs)
            Pmelt_SRarc[np.where(np.asarray(ELAs) == 1.)] = 1.
            
            #get Srad updated to units per hour
            #UNCOMMENT
            SRhour_ARC = np.asarray(SRhour_arc)/solhour_interv
    
            if TOPO[z] == 0:
                MBhour, MBhour_SRarc, Leftover_list, \
                    Leftover_listSR_arc, Melt_list, \
                            Melt_listSR_arc = \
                        MB_341_debris(Thour, Phour, Pmelt, Pmelt_SRarc, \
                            SRhour_ARC, MFi, MFs, asnow, aice, MF, debris_m)
            else:
                MBhour, MBhour_SRarc, Leftover_list, \
                    Leftover_listSR_arc, Melt_list, \
                        Melt_listSR_arc = \
                    MB_341_office(Thour, Phour, Pmelt, Pmelt_SRarc, SRhour_ARC, \
                        MFi, MFs, asnow, aice, MF)      
            
            #Get hourly MB saved
            # 
            # MBwrite = Melt_list
            # for im in MBwrite:
            #     File_out1.write("%s," %im) #prepare for next row/day by next line   
            # File_out1.write("\n")
            # 
            # MBwrite = Melt_listSR_arc
            # for im in MBwrite:
            #     File_out2.write("%s," %im) #prepare for next row/day by next line 
            # File_out2.write("\n")
            
            MBwrite = Thour
            for im in MBwrite:
                File_out3.write("%s," %im) #prepare for next row/day by next line 
            File_out3.write("\n")
            
            # MBwrite = MBhour
            # for im in MBwrite[0]:
            #     File_out4.write("%s," %im) #prepare for next row/day by next line   
            # File_out4.write("\n")
            # 
            # MBwrite = MBhour_SRarc
            # for im in MBwrite[0]:
            #     File_out5.write("%s," %im) #prepare for next row/day by next line 
            # File_out5.write("\n")
            
            
            #UNCOMMENT
            MBwrite = SRhour_ARC
            for im in MBwrite:
                File_out7.write("%s," %im) #prepare for next row/day by next line 
            File_out7.write("\n")
            
            # MBwrite = Leftover_list
            # for im in MBwrite:
            #     File_out8.write("%s," %im) #prepare for next row/day by next line 
            # File_out8.write("\n")
            # 
            # MBwrite = Leftover_listSR_arc
            # for im in MBwrite:
            #     File_out9.write("%s," %im) #prepare for next row/day by next line 
            # File_out9.write("\n")
                
            # MBwrite = precipList
            # for im in MBwrite:
            #     File_out10.write("%s," %im) #prepare for next row/day by next line 
            # File_out10.write("\n")
            
            MBwrite = Phour
            for im in MBwrite:
                File_out11.write("%s," %im) #prepare for next row/day by next line 
            File_out11.write("\n")
    
        
        #Move to next day in netcdf files
        print(fruit)
        sys.stdout.flush()
        
        move = dt.timedelta(days=1)
        fruit = fruit + move
        it_day = it_day + 1
        
    
    
    File_out1.close()
    File_out2.close()
    File_out3.close()
    File_out4.close()
    File_out5.close()
    File_out7.close()
    File_out8.close()
    File_out9.close()
    File_out10.close()
    File_out11.close()
    inT.close()
    inE.close()
    inP.close()
    
    if Generate_snow == True:
        np.savetxt(glacier_id + 'snow' + str(current_year+1) + '.txt', Leftover_listSR_arc)
        print('Snow profile is saved as:' + glacier_id + 'snow' + str(current_year+1) + '.txt')
        
    print(glacier_id + 'Complete')
    
