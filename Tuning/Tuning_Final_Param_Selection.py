# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 08:40:57 2023

Plot snowline tuning scores 

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import heapq
from scipy.stats import norm

years = np.arange(2012,2019+1)

path_to_obs_snowlines = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Snowlines/Rasterized_Observed_Snowlines'
OUTPUT_TUNING_RESULTS = 'D:/TuningOutputs/StageII_Tuning_Results'

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
mb0_1000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_0-1000.csv',delimiter=',')
mb1001_2000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_1001-2000.csv',delimiter=',')
mb2001_3000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_2001-3000.csv',delimiter=',')
mb3001_4000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_3001-4000.csv',delimiter=',')
mb4001_4999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_4001-4999.csv',delimiter=',')
mb5000_5999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_5000-5999.csv',delimiter=',')
mb6000_6999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_6000-6999.csv',delimiter=',')
mb7000_7999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_7000-7999.csv',delimiter=',')
mb8000_8999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_8000-8999.csv',delimiter=',')
mb9000_9999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_9000-9999.csv',delimiter=',')
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
df.to_csv('Tuning_Params_Final_' + model_name + '.csv')
    
# =============================================================================
# Plot the final parameter selection
# =============================================================================

plt.figure(figsize=(12,5))
plt.grid(zorder=20)
#plt.ylim(0,170)
plt.title('Final param selection\nTotal = ' + str(len(mb_under_normaldist)),fontsize=14)
plt.hist(mb_passing,mb_normaldist_bins,rwidth=0.9,color='mediumblue',zorder=10)
plt.hist(mb_under_normaldist,mb_normaldist_bins,rwidth=0.9,color='orange',zorder=20,label='Selected from best snowline scores\nunder normal distribution')
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.xticks(bin_centers,np.round(bin_centers,2),fontsize=14,rotation=45)
x = np.arange(-0.98,0.11,0.02) # aligns with the center of each bin
p = norm.pdf(x, -0.46,0.17)
#plt.ylim(0,30)
#plt.plot(x, p/np.max(p)*np.histogram(mb_aice_geq_asnow,mb_normaldist_bins)[0][4], 'k', linewidth=4,zorder=30)
plt.plot(x, p/np.max(p)*9.6, 'k', linewidth=4,zorder=30)
plt.tight_layout()


plt.figure(figsize=(6,5.5))
plt.title('Final param selection\nTotal = ' + str(len(mb_under_normaldist)),fontsize=14)
plt.scatter(final_scores_norm,mb_passing,c='darkblue',s=20,label='sims within\nMB$_{mean}$ $\pm$ 3$\sigma$')
plt.scatter(snowlines_under_normaldist,mb_under_normaldist,c='red',edgecolor='red',s=20,label='Final params')
#plt.scatter(final_scores_norm[np.where(aice_passing>asnow_passing)],mb_passing[np.where(aice_passing>asnow_passing)],c='darkblue',label='a$_{ice}$ $\geq$ a$_{snow}$')
plt.xlabel('Snowline Score (normalized)',fontsize=14)
plt.ylabel('2007--2018 Mass Balance',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid()
plt.legend(fontsize=14,loc='upper left')
plt.tight_layout()


mb_passing_bins =  np.linspace(-0.46-(3*0.17),-0.46+(3*0.17),26)
bin_centers = np.arange(-0.46-(3*0.17)+(mb_passing_bins[1] - mb_passing_bins[0])/2,-0.46+(3.1*0.17)-(mb_passing_bins[1] - mb_passing_bins[0])/2,(mb_passing_bins[1] - mb_passing_bins[0]))
bin_centers = np.arange(-0.46-(3*0.17)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,-0.46+(3.1*0.17)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0]))


plt.figure(figsize=(7,3))
plt.title('Mass balance for final sims. Total = ' + str(len(mb_under_normaldist)) +'\n Mean = ' + str(np.round(np.mean(mb_under_normaldist),2)) + ' $\pm$ ' + str(np.round(np.std(mb_under_normaldist),2)) +  ' m w.e. yr$^{-1}$',fontsize=14)
plt.hist(mb_under_normaldist,mb_passing_bins,rwidth=0.9,color='mediumblue')
#plt.hist(mb_refmodel,mb_passing_bins,rwidth=0.9,color='red')
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
plt.xticks(bin_centers,np.round(bin_centers,2),fontsize=14,rotation=45)
plt.yticks(fontsize=14)
plt.grid()
plt.tight_layout()

snowline_passing_bins =  np.linspace(np.min(final_scores_norm),np.max(final_scores_norm),26)
#snowline_passing_bins =  np.linspace(np.min(snowlines_under_normaldist),np.max(snowlines_under_normaldist),24)
bin_centers = np.arange(np.min(final_scores_norm)+(snowline_passing_bins[1] - snowline_passing_bins[0])/2,np.max(final_scores_norm)+0.001-(snowline_passing_bins[1] - snowline_passing_bins[0])/2,(snowline_passing_bins[1] - snowline_passing_bins[0]))


plt.figure(figsize=(7,3))
plt.title('Snowline scores for final sims. Total = ' + str(len(snowlines_under_normaldist)) +'\n Mean = ' + str(np.round(np.mean(snowlines_under_normaldist),2)) + ' $\pm$ ' + str(np.round(np.std(snowlines_under_normaldist),3)),fontsize=14)
plt.hist(snowlines_under_normaldist,snowline_passing_bins,rwidth=0.9,color='mediumblue')
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('Snowline score',fontsize=14)
#plt.xticks(np.arange(0.8,0.9,0.01),np.round(np.arange(0.8,0.9,0.01),3),fontsize=14,rotation=45)
plt.xticks(bin_centers,np.round(bin_centers,3),fontsize=14,rotation=45)
plt.yticks(fontsize=14)
plt.xlim(np.min(final_scores_norm),np.max(final_scores_norm))
plt.grid()
plt.tight_layout()

# =============================================================================
# Compare MB from final sims with model run with daily outputs (mb should all be the same!)
# =============================================================================
mb_refmodel = np.loadtxt('D:/Model Runs/REF_MODEL/mb_results_REF_MODEL.csv')
plt.plot(mb_under_normaldist-mb_refmodel)
plt.ylabel('MB from initial runs minus\n MB from final runs',fontsize=14)
plt.xlabel('Sim',fontsize=14)

# =============================================================================
# Plot param distirbutions for top 200% of snowline scores
# =============================================================================
aice_bins = np.linspace(0,2.1e-5,26)
asnow_bins = np.linspace(0,4.7e-6,26)
MF_bins = np.linspace(0,8.9e-4,26)

fig, (ax) = plt.subplots(ncols=3,figsize=(12,5),gridspec_kw={"width_ratios":[1,1,1]},sharey=True)
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.title('a$_{ice}$',fontsize=14)
plt.hist(aice_passing,aice_bins,rwidth=0.9,color='lightgrey',label='Sims within\nMB$_{mean}$ $\pm$ 3$\sigma$')
plt.hist(aice_final,aice_bins,rwidth=0.9,color='mediumblue',label='Final Parameters')
plt.ylim(0,260)
plt.ylabel('Frequency',fontsize=14)
plt.vlines(x=0.000003396,ymin=0,ymax=1200,linestyle='--')
plt.xticks([0,0.000003396],['0','3.4x10$^{-6}$'])
x = np.linspace(np.min(aice), np.max(aice), 100)
p = norm.pdf(x, 0.000003396, 0.00000438)
plt.plot(x, p/np.max(p)*np.max(np.histogram(aice_passing,aice_bins)[0]), 'k', linewidth=2,label='Normal distribution')
plt.margins(x=0)

plt.subplot(1,3,2)
plt.title('a$_{snow}$',fontsize=14)
plt.hist(asnow_passing,asnow_bins,rwidth=0.9,color='lightgrey')
plt.hist(asnow_final,asnow_bins,rwidth=0.9,color='mediumblue')
plt.ylim(0,260)
plt.vlines(x=0.000001546,ymin=0,ymax=1200,linestyle='--')
plt.xticks([0,0.000001546],['0','1.5x10$^{-6}$'])
x = np.linspace(np.min(asnow), np.max(asnow), 100)
p = norm.pdf(x, 0.000001546,0.00000085)
plt.plot(x, p/np.max(p)*np.max(np.histogram(asnow_passing,asnow_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

plt.subplot(1,3,3)
plt.title('$MF$',fontsize=14)
plt.hist(MF_passing,MF_bins,rwidth=0.9,color='lightgrey')
plt.hist(MF_final,MF_bins,rwidth=0.9,color='mediumblue')
plt.ylim(0,260)
plt.vlines(x=0.0002707,ymin=0,ymax=1200,linestyle='--')
plt.xticks([0,0.0002707],['0','2.7x10$^{-4}$'])
x = np.linspace(np.min(MF), np.max(MF), 100)
p = norm.pdf(x, 0.0002707,0.0001632)
plt.plot(x, p/np.max(p)*np.max(np.histogram(MF_passing,MF_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

#handles, labels = ax.get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
fig.legend(fontsize=14,bbox_to_anchor=(0.71,1.1), ncol=3, borderaxespad=0.19)
plt.tight_layout()


# =============================================================================
# =============================================================================
# # Compare final sims from ref model to those from ROUNCE_DEBRIS Model:
# =============================================================================
# =============================================================================
rouncedebris = np.loadtxt('D:/TuningOutputs/ROUNCE_DEBRIS/MB_Final_ROUNCE_DEBRIS.csv',skiprows=1,delimiter=',') 
simID_final_rouncedebris = rouncedebris[:,1]
mb_under_normaldist_rouncedebris = rouncedebris[:,2]
snowlines_under_normaldist_rouncedebris = rouncedebris[:,3]

rouncedebris_params = np.loadtxt('D:/TuningOutputs/ROUNCE_DEBRIS/Params_Final_ROUNCE_DEBRIS.csv',skiprows=1,delimiter=',') 
aice_final_rouncedebris = rouncedebris_params[:,1]
asnow_final_rouncedebris = rouncedebris_params[:,2]
MF_final_rouncedebris = rouncedebris_params[:,3]

# How many sims are in this range for both the ref model and the rounce model?
# Convert arrays to sets and find the intersection
set1 = set(simID_final)
set2 = set(simID_final_rouncedebris)
common_elements = set1.intersection(set2)

# Get the count of common elements
num_common_elements = len(common_elements)

print(f"There are {num_common_elements} common integers between the two arrays.")
    
plt.figure(figsize=(6,6))
plt.title('Final param selection\nTotal = ' + str(len(mb_under_normaldist)),fontsize=14)
#plt.scatter(final_scores_norm,mb_passing,c='darkblue',s=20,label='sims within\nMB$_{mean}$ $\pm$ 3$\sigma$')
plt.scatter(snowlines_under_normaldist,mb_under_normaldist,c='red',edgecolor='red',s=20,label='REF_MODEL')
plt.scatter(snowlines_under_normaldist_rouncedebris,mb_under_normaldist_rouncedebris,c='mediumblue',edgecolor='mediumblue',s=20,label='ROUNCE_DEBRIS')
#plt.scatter(final_scores_norm[np.where(aice_passing>asnow_passing)],mb_passing[np.where(aice_passing>asnow_passing)],c='darkblue',label='a$_{ice}$ $\geq$ a$_{snow}$')
plt.xlabel('Snowline Score (normalized)',fontsize=14)
plt.ylabel('2007--2018 Mass Balance',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid()
plt.xticks(np.arange(0.85,0.9,0.01))
plt.legend(fontsize=14,loc='upper left')
plt.tight_layout()

aice_bins = np.linspace(0,2.1e-5,20)
asnow_bins = np.linspace(0,4.7e-6,20)
MF_bins = np.linspace(0,8.9e-4,20)

fig, (ax) = plt.subplots(ncols=3,figsize=(12,5),gridspec_kw={"width_ratios":[1,1,1]},sharey=True)
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.title('a$_{ice}$',fontsize=14)
plt.hist(aice_final,aice_bins,rwidth=1,color='red',zorder=5,alpha=0.5,label='REF_MODEL')
plt.hist(aice_final_rouncedebris,aice_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5,label='ROUNCE_DEBRIS')
plt.hist(aice_final,aice_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=3,zorder=10)
plt.hist(aice_final_rouncedebris,aice_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=3,zorder=10)
plt.ylim(0,75)
plt.ylabel('Frequency',fontsize=14)
plt.vlines(x=np.mean(aice_final),ymin=0,ymax=1200,linestyle='--',color='red',linewidth=2)
plt.vlines(x=np.mean(aice_final_rouncedebris),ymin=0,ymax=1200,linestyle='--',color='mediumblue',linewidth=2)
plt.vlines(x=np.mean(aice_final_rouncedebris),ymin=-10,ymax=-9,linestyle='--',color='k',linewidth=2,label='Mean param value')
#plt.xticks([0,0.000003396],['0','3.4x10$^{-6}$'])
plt.xticks(np.arange(0,2.2e-5,0.3e-5))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
x = np.linspace(np.min(aice), np.max(aice), 100)
p = norm.pdf(x, 0.000003396, 0.00000438)
plt.plot(x, p/np.max(p)*40, 'k', linewidth=2,label='Normal distribution')
plt.margins(x=0)

plt.subplot(1,3,2)
plt.title('a$_{snow}$',fontsize=14)
plt.hist(asnow_final,asnow_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(asnow_final_rouncedebris,asnow_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
plt.hist(asnow_final,asnow_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=3,zorder=10)
plt.hist(asnow_final_rouncedebris,asnow_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=3,zorder=10)
plt.ylim(0,75)
plt.vlines(x=np.mean(asnow_final),ymin=0,ymax=1200,linestyle='--',color='red',linewidth=2)
plt.vlines(x=np.mean(asnow_final_rouncedebris),ymin=0,ymax=1200,linestyle='--',color='mediumblue',linewidth=2)
plt.xticks(np.arange(0,4.7e-6,0.6e-6))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
x = np.linspace(np.min(asnow), np.max(asnow), 100)
p = norm.pdf(x, 0.000001546,0.00000085)
plt.plot(x, p/np.max(p)*40, 'k', linewidth=2)
plt.margins(x=0)

plt.subplot(1,3,3)
plt.title('$MF$',fontsize=14)
plt.hist(MF_final,MF_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(MF_final_rouncedebris,MF_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
plt.hist(MF_final,MF_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=3,zorder=10)
plt.hist(MF_final_rouncedebris,MF_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=3,zorder=10)
plt.ylim(0,75)
plt.vlines(x=np.mean(MF_final),ymin=0,ymax=1200,linestyle='--',color='red',linewidth=2)
plt.vlines(x=np.mean(MF_final_rouncedebris),ymin=0,ymax=1200,linestyle='--',color='mediumblue',linewidth=2)
plt.xticks(np.arange(0,8.9e-4,1.2e-4))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
x = np.linspace(np.min(MF), np.max(MF), 100)
p = norm.pdf(x, 0.0002707,0.0001632)
plt.plot(x, p/np.max(p)*40, 'k', linewidth=2)
plt.margins(x=0)

fig.legend(fontsize=14,bbox_to_anchor=(0.95,1.1), ncol=4, borderaxespad=0.19)
plt.tight_layout()
