# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 08:40:57 2023

Plot snowline tuning scores 

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import os
import sys
import heapq
from scipy.stats import norm
import pylab
import random

# Import functions for model:
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import shiftedColorMap

years = np.arange(2012,2019+1)

path_to_obs_snowlines = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Snowlines/Rasterized_Observed_Snowlines'
OUTPUT_TUNING_RESULTS = 'D:\TuningOutputs\REF_MODEL\StageII_Tuning_Results'

all_snow_depth_dates = []
for file in os.listdir(path_to_obs_snowlines):
    if file.endswith('.npy'):
        date = pd.Timestamp(year=int(file[15:19]), month=int(file[20:22]), day=int(file[23:25]))
        all_snow_depth_dates.append(date)

# Load mass balance for each sim
# Get the tuning results from stage I:
# =============================================================================
params_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Tuning_Params_Initial10k.csv' # Path to csv file containing MF, aice, asnow.
model_name = 'REF_MODEL'

params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice = params[:,1]
asnow = params[:,2]
MF = params[:,3]
# 
mb0_1000 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_0-1000.csv',delimiter=',')
mb1001_2000 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_1001-2000.csv',delimiter=',')
mb2001_3000 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_2001-3000.csv',delimiter=',')
mb3001_4000 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_3001-4000.csv',delimiter=',')
mb4001_4999 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_4001-4999.csv',delimiter=',')
mb5000_5999 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_5000-5999.csv',delimiter=',')
mb6000_6999 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_6000-6999.csv',delimiter=',')
mb7000_7999 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_7000-7999.csv',delimiter=',')
mb8000_8999 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_8000-8999.csv',delimiter=',')
mb9000_9999 = np.loadtxt('D:/TuningOutputs/REF_MODEL/mb_results_sim_9000-9999.csv',delimiter=',')
mb = np.concatenate((mb0_1000,mb1001_2000,mb2001_3000,mb3001_4000,mb4001_4999,mb5000_5999,mb6000_6999,mb7000_7999,mb8000_8999,mb9000_9999),axis=0)

simID = np.arange(0,len(mb))  

# =============================================================================
# Remove sims where aice < asnow
# =============================================================================
mb_aice_geq_asnow = np.delete(mb,np.where(aice<asnow))
simID_aice_geq_asnow = np.delete(simID,np.where(aice<asnow))

target_mb = -0.46
tagret_mb_stddev = 0.17

# =============================================================================
# Get simIDs and param values for all sims (excluding ones with aice <asnow) where the 2007-18 MB is within -0.46 +- 3std.dev
# =============================================================================
simID_passing = simID_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (target_mb-(3*tagret_mb_stddev))) & (mb_aice_geq_asnow <= (target_mb+(3*tagret_mb_stddev))))]
mb_passing = mb_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (target_mb-(3*tagret_mb_stddev))) & (mb_aice_geq_asnow <= (target_mb+(3*tagret_mb_stddev))))]
aice_passing = aice[simID_passing]
asnow_passing = asnow[simID_passing]
MF_passing = MF[simID_passing]

sim_num = np.arange(0,len(mb_passing))
# =============================================================================
# CALCULATE THE FINAL SCORE FOR EACH IMAGE (SCORE X WEIGHT), FOR EACH SIM
# =============================================================================

all_sim_scores = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_scores[:] = np.nan

all_sim_weights = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_weights[:] = np.nan

all_sim_SA = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_SA[:] = np.nan

all_sim_SW = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_SW[:] = np.nan

all_sim_CA = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_CA[:] = np.nan

all_sim_NA = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_NA[:] = np.nan

all_sim_TR = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_TR[:] = np.nan

all_sim_weightedscores = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_weightedscores[:] = np.nan

all_sim_norm = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_norm[:] = np.nan

for sim in sim_num:
    print(sim)
    scores = np.loadtxt(os.path.join(OUTPUT_TUNING_RESULTS,'StageII_scores_sim' + str(sim) + '.txt'))
    all_sim_scores[sim,:] = scores[0]
    all_sim_weights[sim,:] = scores[1]
    all_sim_SA[sim,:] = scores[2]
    all_sim_SW[sim,:] = scores[3]
    all_sim_CA[sim,:] = scores[4]
    all_sim_NA[sim,:] = scores[5]
    all_sim_TR[sim,:] = scores[6]
    all_sim_weightedscores[sim,:] = (scores[0]*scores[1])
    all_sim_norm[sim,:] = scores[7]
    
# Implement the time-dependent averaging:
def get_final_timeaveraged_scores(score):
    
    time_averaged_score = np.array((np.arange(0,len(mb_passing))))
    
    for day in all_snow_depth_dates:
        #print(day)
        if day == pd.Timestamp(2016,7,4):
            pass
        elif day == pd.Timestamp(2017,7,24):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,7,24))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,7,25))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2017,7,25):
            pass
        
        elif day == pd.Timestamp(2017,8,6):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,8,6))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,8,8))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2017,8,8):
            pass
        
        elif day == pd.Timestamp(2017,9,13):
            pass
        
        elif day == pd.Timestamp(2017,9,20):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,9,20))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,9,23))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2017,9,23):
            pass  
        
        elif day == pd.Timestamp(2018,6,10):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,10))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,14))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,6,14):
            pass  
        
        elif day == pd.Timestamp(2018,6,25):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,25))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,27))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,6,27):
            pass  
        
        elif day == pd.Timestamp(2018,7,22):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,22))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,24))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,25))[0][0]])/3
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,7,24):
            pass  
        elif day == pd.Timestamp(2018,7,25):
            pass 
        
        elif day == pd.Timestamp(2018,7,30):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,30))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,8,1))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,8,1):
            pass 
        
        elif day == pd.Timestamp(2018,8,18):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,8,18))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,8,23))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,8,23):
            pass  
        
        else:
            time_averaged_score = np.column_stack((time_averaged_score,score[:,np.where(np.array(all_snow_depth_dates)==day)][:,0,0]))
        
    return time_averaged_score
        
time_avg_score = get_final_timeaveraged_scores(all_sim_weightedscores)[:,1:]
final_score = np.mean(time_avg_score,axis=1)
 
max_possible_score = get_final_timeaveraged_scores(np.ones((all_sim_scores.shape))*all_sim_weights)[:,1:]
maximum_final_score = np.max(np.mean(max_possible_score,axis=1))
       
#normalize the final scores:
final_scores_norm = np.array(final_score)/maximum_final_score

# Get final scores for each tributary:
time_avg_score = get_final_timeaveraged_scores(all_sim_SA)[:,1:]
final_score_SA = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_SW)[:,1:]
final_score_SW = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_CA)[:,1:]
final_score_CA = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_NA)[:,1:]
final_score_NA = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_TR)[:,1:]
final_score_TR = np.nanmean(time_avg_score,axis=1)

# =============================================================================
# Plot snowline score vs MB
# =============================================================================

plt.figure(figsize=(6,5.5))
plt.title('All sims within target MB $\pm$ 3$\sigma$\nTotal = ' + str(len(mb_passing)),fontsize=14)
plt.scatter(final_scores_norm,mb_passing,c='darkblue',s=20)
plt.xlabel('Snowline Score (normalized)',fontsize=14)
plt.ylabel('2007--2018 Mass Balance',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid()
plt.tight_layout()

# =============================================================================
# Get best snowline scores for sims that fit normal distribution on the mass balance
# =============================================================================

N = int(np.round(np.sqrt(len(mb_passing))))
mb_normaldist_bins =  np.linspace(-0.46-(3*0.17),-0.46+(3*0.17),N)
bin_centers = np.arange(-0.46-(3*0.17)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,-0.46+(3.1*0.17)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0]))

#simID_stageI = simID_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (-0.46-(3*0.17))) & (mb_aice_geq_asnow <= (-0.46+(3*0.17))))]
simID_stageI = np.arange(0,len(mb_passing))

simID_under_normaldist = []
mb_under_normaldist = []
snowlines_under_normaldist = []
for bin_edge in range(0,len(mb_normaldist_bins)-1):
    # Define the normal distribution:
    #x = np.arange(-0.46-(3*0.17)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,-0.46+(3*0.17)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0])) # aligns with the center of each bin
    p = norm.pdf(bin_centers, -0.46,0.17)
    p_scaled =  p/np.max(p)*np.histogram(mb_aice_geq_asnow,mb_normaldist_bins)[0][4]
    p_scaled =  p/np.max(p)*9.6
    
    # Define edges of the bin
    left_edge = np.round(mb_normaldist_bins[bin_edge],2) # bins are defined by [left_edge,right_edge), except the last bin where the right_edge is also included
    right_edge = np.round(mb_normaldist_bins[bin_edge+1],2)
    
    left_edge = mb_normaldist_bins[bin_edge] # bins are defined by [left_edge,right_edge), except the last bin where the right_edge is also included
    right_edge = mb_normaldist_bins[bin_edge+1]
    #print(bin_edge,left_edge,right_edge)
    
    bin_center = np.round((left_edge + right_edge)/2,2)
    print(bin_edge,bin_center)
    
    # Get number of sims in this bin allowed by the normal distribution
    N = int(np.round((p_scaled[np.where(np.round(bin_centers,2)==bin_center)[0][0]])))
    print(bin_center,p_scaled[np.where(np.round(bin_centers,2)==bin_center)[0][0]],N)
    
    if N == 0:
        pass
    else:
        
        # Get the sim IDs and snowline scores corresponding to the mass balance values within each bin 
        simID_bin = simID_stageI[np.where((mb_passing >= (left_edge)) & (mb_passing < (right_edge)))]
        snowline_scores_bin = final_scores_norm[np.where((mb_passing >= (left_edge)) & (mb_passing < (right_edge)))]
        
        # Sample N unique elements from the bin
        if N > len(simID_bin):
            N = len(simID_bin)
            
        # Instead of randomly sampling from inside the bin, get the N best snowline scores from inside the bin
        best_scores = (heapq.nlargest(N, snowline_scores_bin))
        # get Sim IDs corresponding to the best scores
        sims = simID_bin[np.where(snowline_scores_bin >= np.min(best_scores))]
        
        # save the sample to the normal distirbution lists
        simID_under_normaldist.extend(sims)
        mb_under_normaldist.extend(list(mb_passing[sims]))
        snowlines_under_normaldist.extend(list(final_scores_norm[sims]))
  
# Save final sims:
aice_final = aice_passing[simID_under_normaldist]
asnow_final = asnow_passing[simID_under_normaldist]
MF_final = MF_passing[simID_under_normaldist]
simID_final = simID_passing[simID_under_normaldist]

# =============================================================================
# Save final params for this model"
# =============================================================================

d = {'aice': list(aice_final), 'asnow': list(asnow_final), 'MF': list(MF_final)}
df = pd.DataFrame(data=d)
#df.to_csv('Tuning_Params_Final_' + model_name + '.csv')


# =============================================================================
# Value added analysis TEST 1
# Excluding snowline informations
# Sims are randomly selected from binned normal distrubtion (instead of selected by snowline score)
# =============================================================================

N = int(np.round(np.sqrt(len(mb_passing))))
mb_normaldist_bins =  np.linspace(-0.46-(3*0.17),-0.46+(3*0.17),N)
bin_centers = np.arange(-0.46-(3*0.17)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,-0.46+(3.1*0.17)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0]))

#simID_stageI = simID_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (-0.46-(3*0.17))) & (mb_aice_geq_asnow <= (-0.46+(3*0.17))))]
simID_stageI = np.arange(0,len(mb_passing))

simID_under_normaldist_VAAtest1 = []
mb_under_normaldist_VAAtest1 = []
snowlines_under_normaldist_VAAtest1 = []
for bin_edge in range(0,len(mb_normaldist_bins)-1):
    # Define the normal distribution:
    #x = np.arange(-0.46-(3*0.17)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,-0.46+(3*0.17)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0])) # aligns with the center of each bin
    p = norm.pdf(bin_centers, -0.46,0.17)
    p_scaled =  p/np.max(p)*np.histogram(mb_aice_geq_asnow,mb_normaldist_bins)[0][4]
    p_scaled =  p/np.max(p)*9.6
    
    # Define edges of the bin
    left_edge = np.round(mb_normaldist_bins[bin_edge],2) # bins are defined by [left_edge,right_edge), except the last bin where the right_edge is also included
    right_edge = np.round(mb_normaldist_bins[bin_edge+1],2)
    
    left_edge = mb_normaldist_bins[bin_edge] # bins are defined by [left_edge,right_edge), except the last bin where the right_edge is also included
    right_edge = mb_normaldist_bins[bin_edge+1]
    #print(bin_edge,left_edge,right_edge)
    
    bin_center = np.round((left_edge + right_edge)/2,2)
    print(bin_edge,bin_center)
    
    # Get number of sims in this bin allowed by the normal distribution
    N = int(np.round((p_scaled[np.where(np.round(bin_centers,2)==bin_center)[0][0]])))
    print(bin_center,p_scaled[np.where(np.round(bin_centers,2)==bin_center)[0][0]],N)
    
    if N == 0:
        pass
    else:
        
        # Get the sim IDs and snowline scores corresponding to the mass balance values within each bin 
        simID_bin = simID_stageI[np.where((mb_passing >= (left_edge)) & (mb_passing < (right_edge)))]
        snowline_scores_bin = final_scores_norm[np.where((mb_passing >= (left_edge)) & (mb_passing < (right_edge)))]
        
        # Sample N unique elements from the bin
        if N > len(simID_bin):
            N = len(simID_bin)
            
        # randomly select N sims from the bin (simID_bin)
        random.seed(10) # get the same result every time
        sims = random.sample(list(simID_bin),N)
        
        # save the snowline scores and mb from the sample
        simID_under_normaldist_VAAtest1.extend(sims)
        mb_under_normaldist_VAAtest1.extend(list(mb_passing[sims]))
        snowlines_under_normaldist_VAAtest1.extend(list(final_scores_norm[sims]))
  
# Save final sims:
aice_VAAtest1 = aice_passing[simID_under_normaldist_VAAtest1]
asnow_VAAtest1 = asnow_passing[simID_under_normaldist_VAAtest1]
MF_VAAtest1 = MF_passing[simID_under_normaldist_VAAtest1]
simID_VAAtest1 = simID_passing[simID_under_normaldist_VAAtest1]

# =============================================================================
# Save final params for this model"
# =============================================================================

d = {'aice': list(aice_VAAtest1), 'asnow': list(asnow_VAAtest1), 'MF': list(MF_VAAtest1)}
df = pd.DataFrame(data=d)
#df.to_csv('Tuning_Params_VAAtest1_' + model_name + '.csv')
    
# =============================================================================
# Plot the final parameter selection
# =============================================================================

fig = plt.figure(figsize=(9,4))
plt.subplot(1,2,1)
plt.title('a) Reference model',weight='bold',fontsize=14,loc='left')
plt.scatter(final_scores_norm,mb_passing,c='slategrey',s=20)
plt.scatter(snowlines_under_normaldist,mb_under_normaldist,c='navy',edgecolor='navy',s=20)
#plt.scatter(final_scores_norm[np.where(aice_passing>asnow_passing)],mb_passing[np.where(aice_passing>asnow_passing)],c='darkblue',label='a$_{ice}$ $\geq$ a$_{snow}$')
plt.xlabel('Snowline Score (normalized) (-)',fontsize=14)
plt.ylabel('2007-2018\nMass Balance (m w.e. a$^{-1}$)',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(0.8,0.9)
plt.grid()

plt.subplot(1,2,2)
plt.title('b) V.A.A. Test #1',weight='bold',fontsize=14,loc='left')
plt.scatter(final_scores_norm,mb_passing,c='slategrey',s=20,label='Simulations within $\pm$ 3$\sigma$ of $\dot{B}_{obs}$\n(Total = ' + str(len(mb_passing))+ ')')
plt.scatter(snowlines_under_normaldist_VAAtest1,mb_under_normaldist,c='navy',edgecolor='navy',s=20,label='Final simulations\n(Total = 100)')
#plt.scatter(final_scores_norm[np.where(aice_passing>asnow_passing)],mb_passing[np.where(aice_passing>asnow_passing)],c='darkblue',label='a$_{ice}$ $\geq$ a$_{snow}$')
plt.xlabel('Snowline Score (normalized) (-)',fontsize=14)
#plt.ylabel('2007-2018\nMass Balance (m w.e. a$^{-1}$)',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(0.8,0.9)
plt.grid()

fig.legend(fontsize=14,bbox_to_anchor=(0.97,0.95),ncol=2,loc='lower right')
fig.tight_layout()

#fig.savefig('VAAtest1_tuning_scatterplot.pdf',bbox_inches='tight')

# =============================================================================
# How many sims are in common between the ref set and VAA1 set?
# =============================================================================
# Convert arrays to sets and find the intersection
set1 = set(simID_under_normaldist)
set2 = set(simID_under_normaldist_VAAtest1)
common_elements = set1.intersection(set2)
# Get the count of common elements
num_common_elements = len(common_elements)

print(f"There are {num_common_elements} common integers between the two arrays.")
    







# Plot parameter selection for ref model and all 3 VAAA tests:

fig, (ax, ax2, cax) = plt.subplots(ncols=3,figsize=(10,3.6), gridspec_kw={"width_ratios":[1,1,1]})
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.grid(zorder=-20)
#im = plt.scatter(asnow,MF,c=mb,cmap='massbal',edgecolor='k',s=32,vmin=-0.63,vmax=-0.29
#im = plt.scatter(asnow,MF,c='lightgrey',s=32,zorder=10)
plt.scatter(asnow_final,MF_final,c='grey',s=32,zorder=10)
plt.scatter(asnow_VAAtest1,MF_VAAtest1,c='deeppink',s=32,zorder=10)
#plt.scatter(asnow_final_rouncedebris,MF_final_rouncedebris,c='maroon',s=32,zorder=10)
#plt.scatter(asnow_final_uncorrected_acc,MF_final_uncorrected_acc,c='cornflowerblue',s=32,zorder=10)
#plt.scatter(asnow_final_aws_bias,MF_final_aws_bias,c='indigo',s=32,zorder=10)
#plt.xlim(0-0.3e-6,np.max(asnow_final)+0.3e-6)
#plt.ylim(0-0.3e-4,np.max(MF_final)+0.3e-4)
plt.xlim(0,np.max(asnow_final)+0.3e-6)
plt.ylim(0,np.max(MF_final)+0.3e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)
plt.xlabel('a$_{snow}$',fontsize=14)
plt.ylabel('$MF$',fontsize=14)


plt.subplot(1,3,2)
#plt.scatter(aice,MF,c='lightgrey',s=32,zorder=10)
plt.scatter(aice_final,MF_final,c='grey',s=32,zorder=10,label='Reference model')
plt.scatter(aice_VAAtest1,MF_VAAtest1,c='deeppink',s=32,zorder=10,label='V.A.A test #1')
#plt.scatter(aice_final_rouncedebris,MF_final_rouncedebris,c='maroon',s=32,zorder=10,label='Rounce et al. (2021) debris model')
#plt.scatter(aice_final_uncorrected_acc,MF_final_uncorrected_acc,c='cornflowerblue',s=32,zorder=10,label='Model with uncorrected accumulation')
#plt.scatter(aice_final_aws_bias,MF_final_aws_bias,c='indigo',s=32,zorder=10,label='Model with ECCC precipitation bias correction')
plt.xlim(0,np.max(aice_final)+0.4e-6)
plt.ylim(0,np.max(MF_final)+0.3e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)
plt.xlabel('a$_{ice}$',fontsize=14)
plt.ylabel('$MF$',fontsize=14)
plt.grid()

plt.subplot(1,3,3)
#plt.scatter(aice,asnow,c='lightgrey',s=32,zorder=10)
plt.scatter(aice_final,asnow_final,c='grey',s=32,zorder=10)
plt.scatter(aice_VAAtest1,asnow_VAAtest1,c='deeppink',s=32,zorder=10)
#plt.scatter(aice_final_rouncedebris,asnow_final_rouncedebris,c='maroon',s=32,zorder=10)
#plt.scatter(aice_final_uncorrected_acc,asnow_final_uncorrected_acc,c='cornflowerblue',s=32,zorder=10)
#plt.scatter(aice_final_aws_bias,asnow_final_aws_bias,c='indigo',s=32,zorder=10)
plt.xlim(0,np.max(aice_final)+0.4e-6)
plt.ylim(0,np.max(asnow_final)+0.3e-6)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)
plt.xlabel('a$_{ice}$',fontsize=14)
plt.ylabel('a$_{snow}$',fontsize=14)
plt.grid()

fig.legend(fontsize=14,ncol=2,bbox_to_anchor=(1,0.98),loc='lower right')

plt.tight_layout()
#plt.savefig('All_Models_FinalParams_Scatter.pdf',bbox_inches='tight')



fig, axs = plt.subplots(2, 3, figsize=(9.8,6.7))

axs[0,0].text(0.2e-6,36,'a)',fontsize=14,weight='bold')
axs[0,1].text(0.1e-6,36,'b)',fontsize=14,weight='bold')
axs[0,2].text(0.2e-5,36,'c)',fontsize=14,weight='bold')
axs[1,0].text(0.2e-6,36,'d)',fontsize=14,weight='bold')
axs[1,1].text(0.1e-6,36,'e)',fontsize=14,weight='bold')
axs[1,2].text(0.2e-5,36,'f)',fontsize=14,weight='bold')

axs[0,0].grid(zorder=-20)
axs[0,0].hist(aice_final,color='grey',alpha=0.5,zorder=10,label='Reference model')
axs[0,0].hist(aice_final,color='grey',histtype='step',linewidth=3,zorder=10)
axs[0,0].hist(aice_VAAtest1,color='deeppink',alpha=0.5,zorder=10,label='V.A.A. test #1')
axs[0,0].hist(aice_VAAtest1,color='deeppink',histtype='step',linewidth=3,zorder=10)
#axs[0,0].hist(aice_final_rouncedebris,color='maroon',alpha=0.5,zorder=10,label='Rounce et al. (2021) debris model')
#axs[0,0].hist(aice_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
axs[0,0].set_ylabel('Frequency',fontsize=14)
axs[0,0].set_xlabel('a$_{ice}$',fontsize=14)
axs[0,0].set_xlim(0,np.max(aice_final)+0.4e-6)
axs[0,0].set_xticks(np.arange(0,4.01e-6,1e-6))
axs[0,0].tick_params(axis='both',labelsize=14)
axs[0,0].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
axs[0,0].set_ylim(0,35)
axs[0,0].set_yticks(np.arange(0,31,5))

axs[0,1].grid(zorder=-20)
axs[0,1].hist(asnow_final,color='grey',alpha=0.5,zorder=10)
axs[0,1].hist(asnow_final,color='grey',histtype='step',linewidth=3,zorder=10)
axs[0,1].hist(asnow_VAAtest1,color='deeppink',alpha=0.5,zorder=10)
axs[0,1].hist(asnow_VAAtest1,color='deeppink',histtype='step',linewidth=3,zorder=10)
#axs[0,1].hist(asnow_final_rouncedebris,color='maroon',alpha=0.5,zorder=10)
#axs[0,1].hist(asnow_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
axs[0,1].set_ylabel('Frequency',fontsize=14)
axs[0,1].set_xlabel('a$_{snow}$',fontsize=14)
axs[0,1].set_xlim(0,np.max(asnow_final)+0.3e-6)
axs[0,1].set_xticks(np.arange(0,2.01e-6,0.5e-6))
axs[0,1].tick_params(axis='both',labelsize=14)
axs[0,1].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
axs[0,1].set_ylim(0,35)
axs[0,1].set_yticks(np.arange(0,31,5))

axs[0,2].grid(zorder=-20)
axs[0,2].hist(MF_final,color='grey',alpha=0.5,zorder=10)
axs[0,2].hist(MF_final,color='grey',histtype='step',linewidth=3,zorder=10)
axs[0,2].hist(MF_VAAtest1,color='deeppink',alpha=0.5,zorder=10)
axs[0,2].hist(MF_VAAtest1,color='deeppink',histtype='step',linewidth=3,zorder=10)
#axs[0,2].hist(MF_final_rouncedebris,color='maroon',alpha=0.5,zorder=10)
#axs[0,2].hist(MF_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
axs[0,2].set_ylabel('Frequency',fontsize=14)
axs[0,2].set_xlabel('$MF$',fontsize=14)
axs[0,2].set_xlim(0,np.max(MF_final)+0.3e-4)
axs[0,2].set_xticks(np.arange(0,5.01e-4,1e-4))
axs[0,2].tick_params(axis='both',labelsize=14)
axs[0,2].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
axs[0,2].set_ylim(0,36)
axs[0,2].set_yticks(np.arange(0,31,5))

# =============================================================================

axs[1,0].grid(zorder=-20)
axs[0,0].hist(aice_final,color='grey',alpha=0.5,zorder=10,label='Reference model')
axs[0,0].hist(aice_final,color='grey',histtype='step',linewidth=3,zorder=10)
axs[0,0].hist(aice_VAAtest1,color='deeppink',alpha=0.5,zorder=10,label='V.A.A. test #1')
axs[0,0].hist(aice_VAAtest1,color='deeppink',histtype='step',linewidth=3,zorder=10)
#axs[0,0].hist(aice_final_rouncedebris,color='maroon',alpha=0.5,zorder=10,label='Rounce et al. (2021) debris model')
#axs[0,0].hist(aice_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
axs[0,0].set_ylabel('Frequency',fontsize=14)
axs[0,0].set_xlabel('a$_{ice}$',fontsize=14)
axs[0,0].set_xlim(0,np.max(aice_final)+0.4e-6)
axs[0,0].set_xticks(np.arange(0,4.01e-6,1e-6))
axs[0,0].tick_params(axis='both',labelsize=14)
axs[0,0].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
axs[0,0].set_ylim(0,35)
axs[0,0].set_yticks(np.arange(0,31,5))

axs[0,1].grid(zorder=-20)
axs[0,1].hist(asnow_final,color='grey',alpha=0.5,zorder=10)
axs[0,1].hist(asnow_final,color='grey',histtype='step',linewidth=3,zorder=10)
axs[0,1].hist(asnow_VAAtest1,color='deeppink',alpha=0.5,zorder=10)
axs[0,1].hist(asnow_VAAtest1,color='deeppink',histtype='step',linewidth=3,zorder=10)
#axs[0,1].hist(asnow_final_rouncedebris,color='maroon',alpha=0.5,zorder=10)
#axs[0,1].hist(asnow_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
axs[0,1].set_ylabel('Frequency',fontsize=14)
axs[0,1].set_xlabel('a$_{snow}$',fontsize=14)
axs[0,1].set_xlim(0,np.max(asnow_final)+0.3e-6)
axs[0,1].set_xticks(np.arange(0,2.01e-6,0.5e-6))
axs[0,1].tick_params(axis='both',labelsize=14)
axs[0,1].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
axs[0,1].set_ylim(0,35)
axs[0,1].set_yticks(np.arange(0,31,5))

axs[0,2].grid(zorder=-20)
axs[0,2].hist(MF_final,color='grey',alpha=0.5,zorder=10)
axs[0,2].hist(MF_final,color='grey',histtype='step',linewidth=3,zorder=10)
axs[0,2].hist(MF_VAAtest1,color='deeppink',alpha=0.5,zorder=10)
axs[0,2].hist(MF_VAAtest1,color='deeppink',histtype='step',linewidth=3,zorder=10)
#axs[0,2].hist(MF_final_rouncedebris,color='maroon',alpha=0.5,zorder=10)
#axs[0,2].hist(MF_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
axs[0,2].set_ylabel('Frequency',fontsize=14)
axs[0,2].set_xlabel('$MF$',fontsize=14)
axs[0,2].set_xlim(0,np.max(MF_final)+0.3e-4)
axs[0,2].set_xticks(np.arange(0,5.01e-4,1e-4))
axs[0,2].tick_params(axis='both',labelsize=14)
axs[0,2].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
axs[0,2].set_ylim(0,36)
axs[0,2].set_yticks(np.arange(0,31,5))
# 
# =============================================================================

#fig.legend(fontsize=14,ncol=2,bbox_to_anchor=(0.9,0.98),loc='lower right')

fig.tight_layout()

