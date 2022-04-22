######Mass balance model version used to generate melt and radiation factors
######by running a picking random factor values from normal distributions based 
######on literature values. 

########Inputs: Temperature, Precipitation, and geographical information########
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import random
import netCDF4
#from pyproj import Proj
from netCDF4 import Dataset
import datetime as dt
from datetime import datetime
from Model_functions_ver3 import T_downscale_funkfest
from Model_functions_ver3 import Precip_2_Snow
from Model_functions_ver3 import MB_simplified
from Model_functions_ver3 import rainy_day_funk
from Model_functions_ver3 import closest_node
from Model_functions_ver3 import mean_snowpacks_pts
from Model_functions_ver3 import mean_postemps_pts
from Model_functions_ver3 import cold_content_simplified
from Model_functions_ver3 import polyfit_homebrew
import sys
import os
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
import scipy.stats as stats


###########################BEGIN MONTE CARLO HERE###############################

#Set simulation output lists
sim=0
rmse_l = []
mae_l = []
aice_l = []
asnow_l = []
MF_l = []
re_out = []
mre_out = []
r21_l = []
r22_l = []
r23_l = []

year_l = []
duration_l = []
sub1_3_l = []
sub2_3_l = []
loc_l = []
locs_l = []
re_all = []
mre_all = []
#elev_l = []

it_num = 1

###setup cross validation frequency loop
XV = interp1d(np.linspace(5,35,30), np.asarray([   5.8,    6. ,    6. ,   \
        7.7,    9.8,    9.6,   12.2,   16.8, \
         14.7,   16.8,   19.2,   28.4,   27.3,   38. ,   48. ,   50.2,
         54.1,   60.3,   84.2,   99.9,   94.3,  118.1,  148.4,  165.4,
        168.7,  220.9,  276.8,  291.7,  387.2,  588.4]), kind = 'linear')

#run n simulation
I = np.load('Aice_8.npy')
S = np.load('Asnow_8.npy')
Mf = np.load('Mf_8.npy')
Mb = np.load('Mb_8.npy')
ELA = np.load('ELA_8.npy')*100


GMB = -0.46
iv = 0.17
GELA = 2281.
ela_up = 2477.
ela_dn = 1927.

index_list_mb = np.where(np.logical_and(GMB-iv<Mb, Mb<GMB+iv))
index_list_ela = np.where(np.logical_and(ela_dn<ELA, ELA<ela_up))

index_list = np.intersect1d(index_list_mb, index_list_ela)

#####save params that pass the first round of tuning:
#pass_1st_round_list = index_list.tolist()
#np.savetxt('test_param_noBCacc_v2.csv', [np.transpose(I[pass_1st_round_list]), np.transpose(S[pass_1st_round_list]), np.transpose(Mf[pass_1st_round_list])])


#for i in range(0, len(index_list[0])):
#    if ELA
    
    
    

for index in index_list:

    print('\rRun ' + str(sim) + ' started:') 

    sim+=1

    #Model parameters, force positive values and aice greater than asnow

    aice = I[index]
    asnow = S[index]
    MF = Mf[index]
    
    #print 'MF is ' + str(MF) + ', aice is ' + str(aice) + ', asnow is ' + \
     #   str(asnow)
    
    aice_l.append(aice)
    asnow_l.append(asnow)
    MF_l.append(MF)


    ###Output filenames
    run_num = "_aice=" + str(aice) + "_asnow=" + str(asnow) + "_MF=" + str(MF)
    File_melt_out = 'Tuning' + str(run_num) + '.txt'
    run_errors = []
    run_merrors = []
    run_year = []
    run_duration = []
    run_sub1_3 = []
    run_sub2_3 = []
    run_loc = []
    run_loc_s = []
    run_re_all = []
    run_mre_all = []
    #run_elev = []

#######################Usefule time manipulation funktions######################
#these functions are homeade to create time vectors

    str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%m/%d')
        
    def date_linspace(start, end, steps):
        delta = (end - start) / steps
        increments = start + (range(0, steps) * np.array(delta))
        return increments


###########################Timerange############################################

    #run model for ful time period
    earlist = [2006, 2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017, \
        2018]
    #earlist = [2016, 2017]


###########################Get cold content#####################################
    #three functions for getting cold content ONLY WORKS FOR POINT BASED LOCATIONS
    Mean_temp = mean_postemps_pts(earlist[1:])
    Mean_snow, proportions = mean_snowpacks_pts(earlist[1:])
    Cold_contents = cold_content_simplified(earlist[1:], Mean_temp, Mean_snow, \
        proportions)
    CC_4_years = np.asarray(Cold_contents)[0,:,:-1]
    CC_4_years = np.append(CC_4_years, [np.zeros(np.shape(CC_4_years[0,:]))], \
        axis = 0)
    

###########################Begin model HERE#####################################

    for ear in earlist:
    
        This_year_is = ear
    
        ###Set file names###
        ##### EDITS FOR TUNING THE MODEL FOR NO ACCUMULATION BIAS CORRECTION###
        T_inputs = 'F:\\Mass Balance Model\\BiasCorrectedInputs\\Kaskonly_R2S=1'
        P_inputs = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'
        SR_inputs = 'F:\\Mass Balance Model\\Kaskonly_Downscaled_NoBC'
        File_temp_name = 'Temp' + 'kaskonly' + '_BC_' + str(This_year_is) + '.nc'
        File_precip_name = 'netSnow' + 'kaskonly' + str(This_year_is) + '.nc'
        File_PDCSR_name = 'Srad' + 'kaskonly' + str(This_year_is) + '.nc'
        File_temp_in = os.path.join(T_inputs,File_temp_name)
        File_precip_in = os.path.join(P_inputs,File_precip_name)
        File_PDCSR_in = os.path.join(SR_inputs,File_PDCSR_name)
    
        File_tuning_in = 'Tuninglocs' + str(This_year_is) + '.csv'
        #File_temp_in = 'Temptuning' + str(This_year_is) + '.txt'
        #File_precip_in = 'netSnowtuning' + str(This_year_is) + '.txt'
        #File_PDCSR_in = 'PDCSRtuning' + str(This_year_is) + '.txt'
        File_locations_in = 'tuning.csv'
        File_snow_in = 'snowini' + str(This_year_is) + '.txt'
    
        ###Load data: vector form###
        
        tuning = np.genfromtxt(File_tuning_in,  delimiter=',', skip_header = 1)
        TF = tuning[:,7]
        Err = tuning[:,6]
        Tst = np.genfromtxt(File_tuning_in, dtype = None, names = True, \
            delimiter=',', usecols = 4, converters = {4: str2date})
        Tnd = np.genfromtxt(File_tuning_in, dtype = None, names = True, \
            delimiter=',', usecols = 5, converters = {5: str2date})
        #ID locations where GL1 = 1, GL2 = 2 and KASK = 3
        ID = np.zeros(np.shape(TF[:-1]))
        ID[0:17] = 1.; ID[17:35] = 2.; ID[35:44] = 3.
        ID_specific = np.linspace(0, len(ID)-1, len(ID))
    
        for d in Tst:
            d[0] = d[0].replace(year=This_year_is)
        for d in Tnd:
            d[0] = d[0].replace(year=This_year_is)
        
        locations = np.genfromtxt(File_locations_in, delimiter=',')
        Txx = tuning[:,2]
        Tyy = tuning[:,1]
        Tzz = tuning[:,3]
    
        ###Load data: climate###
        
        #T_array = np.genfromtxt(File_temp_in, delimiter = ',')
        #P_array = np.genfromtxt(File_precip_in, delimiter = ',')
        #S_array = np.genfromtxt(File_PDCSR_in, delimiter = ',')
        
        
        inT = Dataset(File_temp_in, "r")
        T_var = 'Temperature'
        T_array = inT.variables[T_var][:]
        #print "Temperature in..."
        sys.stdout.flush()
        
        inP = Dataset(File_precip_in, "r")
        #P_var = 'Net snow'
        P_var = 'Temperature'
        P_array = inP.variables[P_var][:]
        #print "Precipitation in..."
        sys.stdout.flush()
        
        inS = Dataset(File_PDCSR_in, "r")
        #S_var = 'Potential direct clear sky radiation'
        S_var = 'Temperature'
        S_array = inS.variables[S_var][:]
        #print "Potential direct clear-sky radiation in..."
        sys.stdout.flush()
    
    
    
###########################SETUP Model domain###################################

###Prep time and space resolution for model, and arrays used for
###Operations during MB calculations

#arrays for dumping carry-over values for snow in  each model, set initial snow 
#condition here: use previous years sow profile. Also generate leftover lists for
#tracking snowpack through time
        
        #Leftover_list = np.zeros(np.size(Tzz[:-1])) + np.genfromtxt( \
            #File_snow_in, delimiter = ',') 
        #add_snow_ELA = np.where(Tzz[:-1] > 2500)
        #Leftover_list[add_snow_ELA] = 10.
        
        #prep cold content for this year
        CC_4_thisyear = CC_4_years[np.where(np.asarray(earlist) == ear)]
        
        Melthour = []
    
        #set start day for model run
    
        st_day = 1
        it_day = 0
        
        #set up time steps in a day...
        time_delta = 3
        t_st = 0
        t_nd = 24
        timeinterv = np.linspace(t_st, t_nd, t_nd/time_delta, endpoint=False)
    
        #...and set up timesteps through time, with days limited array for indexing
        #perform twice for hourly values and daily precip
        
        ###Date vectors for different datasets###  
        st = dt.datetime(ear, 1, 1, 0, 0, 0)
        nd = dt.datetime(ear+1, 1, 1, 0, 0, 0)
        if ear == 2008 or ear == 2012 or ear == 2016:
            ly = 1
        else: 
            ly = 0
            
        inv = (8 * 365) + (ly * 8)
        delta = dt.timedelta(hours = 3)
        date_list = date_linspace(st, nd, inv)
        dates = date_list
        fruit = dates[(st_day*8) - 8]
        #print "\r", "Start year: ",  
        #print fruit.year,
        sys.stdout.flush()
        
###########################Calculate Mass Balance###############################
        #run for each timestep
        i = 0
        while (fruit < (dates[-1] + dt.timedelta(hours=3))):
            i+=1
            #add cold content in late october to snow pack to prevent winter 
            #melt events
            
            #if fruit.timetuple().tm_yday == 274:
            #    Leftover_list = Leftover_list + CC_4_thisyear[0]
            
            pos = np.where(dates == fruit)[0][0]
                    
            #curT = T_array[pos,:-1]        
            #curS = S_array[pos,:-1]
            #units conversion to m.w.e. for precipitation
            #curP = P_array[pos,:-1]/1000
            curT = T_array[pos,:,:]        
            curS = S_array[pos,:,:]
            ###units conversion to m.w.e. for precipitation done in DS step
            curP = P_array[pos,:,:]
            ###update snowpack with acumulation
            #LefLeftotover_list = Leftover_list + curP
            
            snowlocs = np.where(curT > 1.)
            curP[snowlocs] = 0.
            
            
            #determine current temporal position
            #pos = np.where(dates == fruit)[0][0]
            Topo = np.ones(T_array[1,:,:].shape)
            CC_list = np.zeros(T_array[1,:,:].shape)
            
            if i == 1:
                Leftover_list = np.zeros(T_array[1,:,:].shape)
            else:
                Leftover_list = np.asarray(Leftover_list_out)

            ###calculate mass balance
            print(fruit)
            MBhour, Melt, Leftover_list_out, Icemelt, Snowmelt, CC_out = MB_vectorized_discreteSnI(curT, curP, curS, Leftover_list, asnow, aice, MF, Topo, CC_list)
            
                    
            #MBhour, Melt, Leftover_list_out = MB_simplified(curT[0][0], curP[0][0], curS[0][0], Leftover_list, asnow, aice, MF)
                
            #MBhour, Melt, Leftover_list_out = MB_simplified(curT, curP, curS, Leftover_list, asnow, aice, MF)
                        
                            
            Melthour.append(Melt)
            Leftover_list = np.asarray(Leftover_list_out)
            
            #print "\r", fruit,
                    
            move = dt.timedelta(hours=3)
            fruit = fruit + move
        
        #generate snow content for next year
        np.savetxt('snowini' + str(This_year_is+1) + '.txt', Leftover_list)
        MELTS = np.asarray(Melthour)

        
            
    ################Calculate error for each melt factor combination################
        #For the year
        print('calculating errors')
        errors = []
        merrors = []
        
        year_y = []
        duration_y = []
        sub1_3_y = []
        sub2_3_y = []
        loc_y = []
        loc_s_y = []
        #elev_y = []
        
    
        tunes = np.where(TF[:-1] == 1) 
        
        #drop year if measured melt n < 3 because x validation needs integer 
        #value of at least 3
        if np.sum(TF[:-1]) < 3:
            pass
        else:
            #for cross validation loop, n = 50 x validation runs per year:
             #XV(len(tunes[0]*0.66)) for xv loop len scaled to sample size
            for r in range(0, 1):
                
                #randomly select 2/3 of melt targets    
                Melt_indexes = np.random.choice(tunes[0], int(len(tunes[0])), replace=False)
                Holdout_indices = np.delete(tunes[0], np.in1d(tunes[0], Melt_indexes).nonzero())
    
                mini_errors = []
                mini_merrors = [] 
                mini_year_y = []
                mini_duration_y = []
                #mini_elev_y = []
                      
                #for each selected target
                for t in Melt_indexes:
                    
                    #set cumulative melt and get time intervals for target
                    #Melt_cum = np.cumsum(MELTS[:,t]) - P_array[:,t] #OG
                    Melt_cum = np.cumsum(np.nanmean(MELTS[:,t],axis=1)) - np.cumsum(np.nanmean(P_array[:,t],axis=1))
                    #Melt_cum = np.cumsum(MELTS[:,t] - P_array[:,t]) #NEW
                    
                    
                    stday = Tst[t][0]
                    ndday = Tnd[t][0]
                    m_err = Err[t]
                    
                    Melt_st = np.where(dates == stday)
                    Melt_nd = np.where(dates == ndday)
            
                    Melt_interv = Melt_cum[Melt_nd] - Melt_cum[Melt_st]
            
                    Error = Melt_interv - m_err
                    
                    #segment to test index errors in date vector
                    # w = 0
                    # if Melt_interv == 0:
                    #     print stday
                    #     print ndday
                    #     print Melt_nd
                    #     print ID_specific[t]
                    #     w = w + 1
                        
                        
                    
                    #saves the melt values at each interval 1: modelled, 2: obs              
                    mini_errors.append(Melt_interv)
                    mini_merrors.append(m_err)
                    mini_year_y.append(This_year_is)
                    mini_duration_y.append((ndday - stday).days) 
                    #mini_elev_y.append(Tzz[t])
                    
                #Get error here from x validation run and append to year errors
                errors.append(mini_errors)
                merrors.append(mini_merrors)
                year_y.append(mini_year_y)
                duration_y.append(mini_duration_y)
                sub1_3_y.append(Holdout_indices)
                sub2_3_y.append(Melt_indexes)
                loc_y.append(ID[Melt_indexes])
                loc_s_y.append(ID_specific[Melt_indexes])
                #elev_y.append(mini_elev_y)
                
            #Get overall error from year here and save for overall run
            run_errors.append(errors)
            run_merrors.append(merrors)
            run_year.append(year_y)
            run_duration.append(duration_y)
            run_sub1_3.append(sub1_3_y)
            run_sub2_3.append(sub2_3_y)
            run_loc.append(loc_y)
            run_loc_s.append(loc_s_y)
            #run_elev.append(elev_y)
            
    #get error for the run here
    #mega nested loop monstrosity to flatten n-dimensional lists
    re_l = []
    mre_l = []
    y_l = []
    d_l = []
    s13 = []
    s23 = []
    loc = []
    loc_s = []
    #ELEV = []
    for j in range(0, len(run_errors)):
        for jj in range(0, len(run_errors[j])):
            for jjj in range(0, len(run_errors[j][jj])):
                re_l.append(run_errors[j][jj][jjj])
                mre_l.append(run_merrors[j][jj][jjj])
                y_l.append(int(run_year[j][jj][jjj]))
                d_l.append(int(run_duration[j][jj][jjj]))
                s23.append(int(run_sub2_3[j][jj][jjj]))
                loc.append(int(run_loc[j][jj][jjj]))
                loc_s.append(int(run_loc_s[j][jj][jjj]))
                #ELEV.append(int(run_elev[j][jj][jjj]))
                
    for j in range(0, len(run_sub1_3)):
        for jj in range(0, len(run_sub1_3[j])):
            for jjj in range(0, len(run_sub1_3[j][jj])):
                s13.append(int(run_sub1_3[j][jj][jjj]))

    
    #get mean modelled and mean obs melt for the run for out
    re = np.asarray(re_l)
    mre = np.asarray(mre_l)
    re_out.append(np.mean(re))
    mre_out.append(np.mean(mre))
    
    #calc rmse and other metrics for the run for out
    rmse = np.sqrt(mean_squared_error(re, mre))
    mae = mean_absolute_error(re, mre)
    rmse_l.append(rmse)
    mae_l.append(mae) 
    
    #get r2 for a given degree, note re needs some coaxing to be 1D array: [:,0]
    #degree = [1,2,3]
    #for d in degree: 
    r2 = polyfit_homebrew(re[:,0], mre, 1)  
    r21_l.append(r2)
    r2 = polyfit_homebrew(re[:,0], mre, 2)  
    r22_l.append(r2)
    r2 = polyfit_homebrew(re[:,0], mre, 3)  
    r23_l.append(r2)
    
    #generate a scatterplot to see
    plt.figure()
    plt.scatter(mre, re, c = np.asarray(loc))
    plt.plot([0,4] , [0,4], '-k' )
    plt.xlim([0,6])
    plt.ylim([0,6])
    plt.xlabel('', fontsize = 18)
    plt.xlabel('Obs melt (m.w.e.)', fontsize = 18)
    plt.ylabel('Mod melt (m.w.e.)', fontsize = 18)
    plt.savefig(str(run_num) + '.png')
    plt.close()
    
    #track the mega metrics
    re_all.append(re_l)
    mre_all.append(mre_l)
    year_l.append(y_l)
    duration_l.append(d_l)
    sub1_3_l.append(s13)
    sub2_3_l.append(s23)
    loc_l.append(loc) 
    locs_l.append(loc_s)  
    #elev_l.append(ELEV)

#outputs small enough to just write n outputs x len = 500 lists using np.savetxt 
#m prefix is observed!!! Measured!!!! whatever not modelled. Its stupid I know.       
np.savetxt('rmse' + str(it_num) + '.txt', rmse_l)  
np.savetxt('mae' + str(it_num) + '.txt', mae_l)    
np.savetxt('MF' + str(it_num) + '.txt', MF_l) 
np.savetxt('aice' + str(it_num) + '.txt', aice_l) 
np.savetxt('asnow' + str(it_num) + '.txt', asnow_l) 
np.savetxt('re' + str(it_num) + '.txt', re_out) 
np.savetxt('mre' + str(it_num) + '.txt', mre_out) 
np.savetxt('r21' + str(it_num) + '.txt', r21_l) 
np.savetxt('r22' + str(it_num) + '.txt', r22_l) 
np.savetxt('r23' + str(it_num) + '.txt', r23_l) 
np.savetxt('re_all' + str(it_num) + '.txt', np.asarray(re_all)[:,:,0])
np.savetxt('mre_all' + str(it_num) + '.txt', mre_all)
np.savetxt('year' + str(it_num) + '.txt', year_l)
np.savetxt('duration' + str(it_num) + '.txt', duration_l)
np.savetxt('sub_1_3' + str(it_num) + '.txt', sub1_3_l)
np.savetxt('sub2_3' + str(it_num) + '.txt', sub2_3_l)
np.savetxt('loc' + str(it_num) + '.txt', loc_l)
np.savetxt('locs' + str(it_num) + '.txt', locs_l)
#np.savetxt('elev.txt', elev_l)

#############################PLOT THE WHOLE SHEBANG#############################

rmse_l = np.loadtxt('rmse1.txt')
MF_l = np.loadtxt('MF.txt')
aice_l = np.loadtxt('aice.txt')
asnow_l = np.loadtxt('asnow.txt')
re_l = np.loadtxt('re1.txt')
mre_l = np.loadtxt('mre1.txt')
r21_l = np.loadtxt('r21.txt')
r22_l = np.loadtxt('r22.txt')
r23_l = np.loadtxt('r23.txt')


#rmse_l, MF_l, aice_l, asnow_l, re_l, mre_l, r21_l, r22_l, r23_l = zip(*sorted(zip(rmse_l, MF_l, aice_l, asnow_l, re_l, mre_l, r21_l, r22_l, r23_l)))
    

#aice_ln = aice_l * (MF_l/aice_l)

########################LINE PLOPT######################
# 
fig = plt.figure(figsize=(10,10))
#f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10,10))

ax1 = fig.add_subplot(311)
ax1.plot(rmse_l, MF_l, '-o', c = 'r', linewidth = 2, label = 'MF')
ax1.set_ylabel('Melt Factor (m.w.e. 3hr-1 C-1)', fontsize = 18, color = 'r')
ax1.tick_params(axis='both', which='major', labelsize=18)
ax1.tick_params(axis='both', which='major', labelsize=18)

for tl in ax1.get_yticklabels():
    tl.set_color('r')


ax2 = ax1.twinx()
ax2.plot(rmse_l, aice_l, '-o', c = 'b', linewidth = 2, label = 'aice')
ax2.plot(rmse_l, asnow_l, '-o', c = 'm', linewidth = 2, label = 'asnow')
ax2.set_ylabel('Radiation factor (m.w.e. 3hr-1 C-1 m2 W-1)', color='b', \
    fontsize = 18)
#ax2.set_xticklabels([])
ax2.tick_params(axis='both', which='major', labelsize=18)

for tl in ax2.get_yticklabels():
    tl.set_color('b')
    
#ax2.set_xlim(0.8, 0.9)
    
ax3 = fig.add_subplot(312)
ax3.fill_between(rmse_l, mre_l, re_l, color = 'cadetblue', alpha = 0.4,label = 'Melt residual')
ax3.plot(rmse_l, mre_l, color = 'r', label = 'Measured melt')
ax3.plot(rmse_l, re_l, color = 'b', label = 'Modelled melt')
ax3.set_xlabel('RMSE (m.w.e.)', fontsize = 18)
ax3.set_ylabel('Mean melt (m.w.e.)', fontsize = 18)
ax3.tick_params(axis='both', which='major', labelsize=18)
#ax3.set_xlim(0.8, 0.9)
plt.legend(fontsize = 18)

ax4 = fig.add_subplot(313)
ax4.plot(rmse_l, r21_l, color = 'darkseagreen', label = 'deg = 1')
ax4.plot(rmse_l, r22_l, color = 'lightgreen', label = 'deg = 2')
ax4.plot(rmse_l, r23_l, color = 'limegreen', label = 'deg = 3')
ax4.set_xlabel('RMSE (m.w.e.)', fontsize = 18)
ax4.set_ylabel('r2', fontsize = 18)
ax4.tick_params(axis='both', which='major', labelsize=18)
#ax3.set_xlim(0.8, 0.9)


   
#ax1.set_xticklabels(fontsize = 18)
ax1.set_xlabel('RMSE (m.w.e.)', fontsize = 18)
plt.legend(fontsize = 18)


    
modelled_netmelt = np.transpose(np.sum(np.asarray(re_all), axis = 1))[0]
measured_netmelt = np.sum(np.asarray(mre_all), axis = 1)
diff = (modelled_netmelt-measured_netmelt)/455
plt.scatter(diff, rmse_l)
np.save('diff.npy',diff)
diff = np.load('diff.npy')

plt.close()


#############################Plot RMSE/RMSEP scatter ###########################

RMSE = np.loadtxt('rmse1.txt')
MAE = np.loadtxt('mae1.txt')
modelled = np.loadtxt('re_all1.txt')
measured = np.loadtxt('mre_all1.txt')

###normalize by days
length = np.loadtxt('duration1.txt')
norm_modelled = modelled/length
norm_measured = measured/length

time_length = length/np.nansum(length, axis = 0)

time_length_fac = (length/np.min(length))
time_total = np.sum(time_length_fac, axis = 0)

weighted_modelled = (norm_modelled*time_length_fac)
weighted_measured = (norm_measured*time_length_fac)


#norm_rmse = np.sqrt(np.sum(np.square(np.diff((norm_modelled, norm_measured), axis = 0)), axis = 2)/len(norm_modelled[0])) #OG!!
#norm_mae = np.sum(np.abs(np.diff((norm_modelled, norm_measured), axis = 0)), axis = 2)/len(norm_modelled[0]) #OG!!!
norm_rmse = np.sqrt(np.sum(np.square(np.diff((norm_modelled, norm_measured), axis = 0)))/len(norm_modelled))
norm_mae = np.sum(np.abs(np.diff((norm_modelled, norm_measured), axis = 0)))/len(norm_modelled)

diffs = []
abs_diffs = []
maxes = []
varis = []
stds = []
for i in range(0, len(modelled)):

    diff_fac = np.median((modelled[i] - measured[i])/measured[i])
    abs_diff_fac = np.median(np.abs((modelled[i] - measured[i]))/measured[i])   
    max_fac = np.max((modelled[i] - measured[i]))
    var_fac = np.var((modelled[i] - measured[i])/measured[i])
    std_fac = np.std((modelled[i] - measured[i])/measured[i])
    diffs.append(diff_fac)
    abs_diffs.append(abs_diff_fac)
    maxes.append(max_fac)
    varis.append(var_fac)
    stds.append(std_fac)

plt.scatter(MAE, RMSE)
plt.xlabel('MAE')
plt.ylabel('RMSE')
plt.xlim([0.006,0.016])
plt.ylim([0.006,0.02])
plt.grid()

plt.figure(figsize=(8,6))
#this is the plot from fig6 in Eriks paper:
plt.scatter(norm_mae, norm_rmse, c = 'k')
plt.title('Debris-Free Simulations')
plt.xlabel('MAE')
plt.ylabel('RMSE')
plt.grid()
plt.axhline(y=0.01,linewidth=4,xmin=0,xmax=0.4,color='red',linestyle='--')
plt.axhline(y=0.008,linewidth=4,xmin=0,xmax=0.4,color='red',linestyle='--')
plt.axvline(x=0.01,linewidth=4,ymin=0,ymax=0.17,color='red',linestyle='--')
plt.axhline(y=0.0099,linewidth=4,xmin=0,xmax=0.39,color='green',linestyle='--')
plt.axhline(y=0.008,linewidth=4,xmin=0,xmax=0.39,color='green',linestyle='--')
plt.axvline(x=0.0099,linewidth=4,ymin=0,ymax=0.16,color='green',linestyle='--')
plt.xlim(0.006, 0.016)
plt.ylim(0.008, 0.02)


plt.scatter(norm_mae, norm_rmse)
plt.xlabel('MAE')
plt.ylabel('RMSE')
plt.grid()

plt.figure(figsize=(8,6))
#this is also the plot from fig6:
fig = plt.scatter(diffs, abs_diffs, c = 'k')
plt.title('Debris-Free Simulations')
plt.xlabel('Median relative difference')
plt.ylabel('Median Absolute relative difference')
plt.grid()
plt.axhline(y=0.3,linewidth=4,xmin=0.333,xmax=0.555,color='red',linestyle='--')
plt.axhline(y=0.5,linewidth=4,xmin=0.333,xmax=0.554,color='red',linestyle='--')
plt.axvline(x=-0.2,linewidth=4,ymin=0,ymax=0.4,color='red',linestyle='--')
plt.axvline(x=0.2,linewidth=4,ymin=0,ymax=0.4,color='red',linestyle='--')
plt.axhline(y=0.3,linewidth=4,xmin=0.388,xmax=0.499,color='green',linestyle='--')
plt.axhline(y=0.44,linewidth=4,xmin=0.388,xmax=0.499,color='green',linestyle='--')
plt.axvline(x=-0.1,linewidth=4,ymin=0,ymax=0.3,color='green',linestyle='--')
plt.axvline(x=0.1,linewidth=4,ymin=0,ymax=0.3,color='green',linestyle='--')
plt.xlim(-0.75, 1.)
plt.ylim(0.3, 0.8)

plt.scatter(stds, RMSE)
plt.xlabel('Maximum melt difference')
plt.ylabel('RMSE')
plt.grid()


##get weighted median HAHA!
#diffs = []
#abs_diffs = []
#maxes = []
#varis = []
#stds = []
#for i in range(0, len(modelled)):
#
#    diff_fac = np.median(((weighted_modelled[i] - weighted_measured[i])/weighted_measured[i])*time_length_fac[i])
#    abs_diff_fac = np.median((np.abs((weighted_modelled[i] - weighted_measured[i]))/weighted_measured[i])*time_length_fac[i])   
#    max_fac = np.max((weighted_modelled[i] - weighted_measured[i]))
#    var_fac = np.var((weighted_modelled[i] - weighted_measured[i])/weighted_measured[i])
#    std_fac = np.std((weighted_modelled[i] - weighted_measured[i])/weighted_measured[i])
#    diffs.append(diff_fac)
#    abs_diffs.append(abs_diff_fac)
#    maxes.append(max_fac)
#    varis.append(var_fac)
#    stds.append(std_fac)
#    
#plt.scatter(diffs, abs_diffs)
#plt.xlabel('Median relative difference')
#plt.ylabel('Median Absolute relative difference')
#plt.grid()

###create lists of parameter values for finall runs

#index_list_rmse = np.where(np.logical_and(np.asarray(RMSE) < 1., np.asarray(MAE) < 0.8))
#index_list_pdiff = np.where(np.logical_and(np.asarray(diffs) > -0.2, np.asarray(diffs) < 0.2))
#index_list_abspdiff = np.where(np.asarray(abs_diffs) < 0.5)
#
#
#final_list = np.intersect1d(index_list_rmse, index_list_pdiff, index_list_abspdiff)


###final indices using updated metircs

# OG CONDITIONS FOR PASSING TUNING:
index_list_rmse = np.where(np.logical_and(np.asarray(norm_rmse[0]) < 0.01, np.asarray(norm_mae[0]) < 0.01))
index_list_pdiff = np.where(np.logical_and(np.asarray(diffs) > -0.2, np.asarray(diffs) < 0.2))
index_list_abspdiff = np.where(np.asarray(abs_diffs) < 0.5)

#index_list_rmse = np.where(np.logical_and(np.asarray(norm_rmse) < 0.01, np.asarray(norm_mae) < 0.01))
#index_list_pdiff = np.where(np.logical_and(np.asarray(diffs) > -0.2, np.asarray(diffs) < 0.2))
#index_list_abspdiff = np.where(np.asarray(abs_diffs) < 0.5)

# conditions which make 10 params
#index_list_rmse = np.where(np.logical_and(np.asarray(norm_rmse[0]) < 0.0096, np.asarray(norm_mae[0]) < 0.0096))
#index_list_pdiff = np.where(np.logical_and(np.asarray(diffs) > -0.1, np.asarray(diffs) < 0.1))
#index_list_abspdiff = np.where(np.asarray(abs_diffs) < 0.45)


final_list_first = np.intersect1d(index_list_rmse, index_list_pdiff)
final_list = np.intersect1d(final_list_first, index_list_abspdiff)

#covert to list for random sampling
finallist = final_list.tolist()
param_sample = random.sample(finallist,10)
#save in param file
aicel = np.loadtxt('aice1.txt')
asnowl = np.loadtxt('asnow1.txt')
MFl = np.loadtxt('MF1.txt')

np.savetxt('final_params_noBCacc.csv', [np.transpose(aicel[final_list]), np.transpose(asnowl[final_list]), np.transpose(MFl[final_list])])

#plot the distribution of the new params


###plots of parameter combinations used

aice = np.random.normal(0.000003396, 0.00000438)
asnow = np.random.normal(0.000001546,0.00000085)
MF = np.random.normal(0.0002707,0.0001632)


fig, axs = plt.subplots(1, 3, sharey=True, figsize = (5,2))
#aice
mu = 0.000003396
sigma = 0.00000438
x = np.linspace(0, mu + 3*sigma, 100)
axs[0].hist(I, bins = 12, color = 'gray', rwidth = 0.9)
axs[0].plot(x, (stats.norm.pdf(x, mu, sigma)/np.max(stats.norm.pdf(x, mu, sigma)))*190, linewidth = 3, color = 'k')
axs[0].set_ylabel('Frequency')
axs[0].set_xlim([0,mu + 3*sigma])
axs[0].axvline(mu, linestyle = '--', c = 'r')

#aice
mu = 0.000001546
sigma = 0.00000085
x = np.linspace(0, mu + 3*sigma, 100)
axs[1].hist(S, bins = 12,color = 'gray', rwidth = 0.9)
axs[1].plot(x, (stats.norm.pdf(x, mu, sigma)/np.max(stats.norm.pdf(x, mu, sigma)))*166, linewidth = 3, color = 'k')
axs[1].set_xlabel('Parameter value (m w.e./3h C [m2/W])')
axs[1].set_xlim([0,mu + 3*sigma])
axs[1].axvline(mu, linestyle = '--', c = 'r')

#aice
mu = 0.0002707
sigma = 0.0001632
x = np.linspace(0, mu + 3*sigma, 100)
axs[2].hist(Mf, bins = 12,color = 'gray', rwidth = 0.9)
axs[2].plot(x, (stats.norm.pdf(x, mu, sigma)/np.max(stats.norm.pdf(x, mu, sigma)))*203, linewidth = 3, color = 'k')
axs[2].set_xlim([0,mu + 3*sigma])
axs[2].axvline(mu, linestyle = '--', c = 'r')







