# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 15:37:13 2024

Generate combinations of the melt model parameters a_ice, a_snow, MF
for tuning the mass balance model with an Environment Canada AWS derived
accumulation bias correction. The range of values for each parameter will be
limited to simulations that produced modelled mass balances within +- 4 std. dev.
of the target MB, to speed up the rate at which new parameter combinations pass the 
tuning process. 

@author: katierobinson
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm


# =============================================================================
# Load the initial 20k runs with the AWS bias correction model, determine
# the simulations that produced modelled mass balances within +- 4 std. dev
# =============================================================================
params_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Tuning_Params_Initial10k.csv' # Path to csv file containing MF, aice, asnow.
model_name = 'AWS_BIAS'

params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice1 = params[:,1]
asnow1 = params[:,2]
MF1 = params[:,3]
# 
mb1 = np.loadtxt('D:/TuningOutputs/AWS_BIAS/StageI/mb_results_0-9999.csv',delimiter=',')

simID1 = np.arange(0,len(mb1))  


################################################

params_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Tuning_Params_Extra10k.csv' # Path to csv file containing MF, aice, asnow.

params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice2 = params[:,1]
asnow2 = params[:,2]
MF2 = params[:,3]
# 
mb2 = np.loadtxt('D:/TuningOutputs/AWS_BIAS/ExtraRuns/mb_results_10k-20k.csv',delimiter=',')

simID2 = np.arange(0,len(mb2)) 

aice = np.concatenate((aice1,aice2))
asnow = np.concatenate((asnow1,asnow2))
MF = np.concatenate((MF1,MF2))
mb = np.concatenate((mb1,mb2))
simID =  np.arange(0,len(mb)) 

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

simID_passing_4sigma = simID_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (target_mb-(4*tagret_mb_stddev))) & (mb_aice_geq_asnow <= (target_mb+(4*tagret_mb_stddev))))]
mb_passing_4sigma = mb_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (target_mb-(4*tagret_mb_stddev))) & (mb_aice_geq_asnow <= (target_mb+(4*tagret_mb_stddev))))]
aice_passing_4sigma = aice[simID_passing_4sigma]
asnow_passing_4sigma = asnow[simID_passing_4sigma]
MF_passing_4sigma = MF[simID_passing_4sigma]

# =============================================================================
# Plot histogram of params that were tested
# =============================================================================
N = int(np.sqrt(len(mb_passing)))
#aice_bins = np.arange(0,2.1e-5,0.2e-5)
#asnow_bins = np.arange(0,4.7e-6,0.4e-6)
#MF_bins = np.arange(0,8.9e-4,0.8e-4)

aice_bins = np.linspace(0,2.1e-5,N)
asnow_bins = np.linspace(0,4.7e-6,N)
MF_bins = np.linspace(0,8.9e-4,N)

aicestd = 0.00000438
asnowstd = 0.00000085
MFstd = 0.0001632

fig, (ax) = plt.subplots(ncols=3,figsize=(10,3.5),gridspec_kw={"width_ratios":[1,1,1]},sharey=True)
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.title('a$_{ice}$',fontsize=14)
#plt.hist(aice,aice_bins,rwidth=1,color='lightgrey',label='All simulations')
#plt.hist(aice_passing_allsims,aice_bins,rwidth=1,color='red',zorder=5,alpha=0.5,label='$\dot{b}_{mod}$ = $\dot{b}_{obs}$ $\pm$ 3$\sigma$')
plt.hist(aice_passing,aice_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5,label='$\dot{b}_{mod}$ = $\dot{b}_{obs}$ $\pm$ 3$\sigma$\nand a$_{ice}$ $\geq$ a$_{snow}$',align='mid')
#plt.hist(aice_passing_allsims,aice_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(aice_passing,aice_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10,align='mid')
plt.ylim(0,200)
plt.ylabel('Frequency',fontsize=14)
plt.vlines(x=0.000003396,ymin=0,ymax=2500,linestyle='--')
plt.vlines(x= max(aice_passing_4sigma),ymin=0,ymax=2500,linestyle='--',color='cyan',linewidth=3,label='Maximum param from\n$\dot{b}_{mod}$ = $\dot{b}_{obs}$ $\pm$ 4$\sigma$\nand a$_{ice}$ $\geq$ a$_{snow}$')
#plt.xticks([0,0.000003396,0.000003396+(1*aicestd),0.000003396+(2*aicestd),0.000003396+(3*aicestd)],['0','3.4x10$^{-6}$','7.8x10$^{-6}$','1.2x10$^{-5}$','1.6x10$^{-5}$'],rotation=45)
plt.xticks([0,0.000003396,0.000003396+(3*aicestd)],['0','3.4x10$^{-6}$','1.6x10$^{-5}$'],rotation=0)
x = np.linspace(np.min(aice), np.max(aice), 100)
p = norm.pdf(x, 0.000003396, 0.00000438)
plt.plot(x, p/np.max(p)*np.max(np.histogram(aice_passing,aice_bins)[0]), 'k', linewidth=2,label='Normal distribution')
plt.margins(x=0)


plt.subplot(1,3,2)
plt.title('a$_{snow}$',fontsize=14)
#plt.hist(asnow,asnow_bins,rwidth=1,color='lightgrey')
#plt.hist(asnow_passing_allsims,asnow_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(asnow_passing,asnow_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
#plt.hist(asnow_passing_allsims,asnow_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(asnow_passing,asnow_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10)
plt.ylim(0,100)
plt.vlines(x=0.000001546,ymin=0,ymax=2500,linestyle='--')
plt.vlines(x= max(asnow_passing_4sigma),ymin=0,ymax=2500,linestyle='--',color='cyan',linewidth=3)
plt.xticks([0,0.000001546,0.000001546+(3*asnowstd)],['0','1.5x10$^{-6}$','4.1x10$^{-6}$'])
x = np.linspace(np.min(asnow), np.max(asnow), 100)
p = norm.pdf(x, 0.000001546,0.00000085)
plt.plot(x, p/np.max(p)*np.max(np.histogram(asnow_passing,asnow_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

plt.subplot(1,3,3)
plt.title('$MF$',fontsize=14)
#plt.hist(MF,MF_bins,rwidth=1,color='lightgrey')
#plt.hist(MF_passing_allsims,MF_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(MF_passing,MF_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
#plt.hist(MF_passing_allsims,MF_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(MF_passing,MF_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10)
plt.ylim(0,100)
plt.vlines(x=0.0002707,ymin=0,ymax=2500,linestyle='--')
plt.vlines(x= max(MF_passing_4sigma),ymin=0,ymax=2500,linestyle='--',color='cyan',linewidth=3)
plt.xticks([0,0.0002707,0.0002707+(3*MFstd)],['0','2.7x10$^{-4}$','7.6x10$^{-4}$'])
x = np.linspace(np.min(MF), np.max(MF), 100)
p = norm.pdf(x, 0.0002707,0.0001632)
plt.plot(x, p/np.max(p)*np.max(np.histogram(MF_passing,MF_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

#handles, labels = ax.get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
fig.legend(fontsize=14,bbox_to_anchor=(1.2,0.8), ncol=1, borderaxespad=0.19)
plt.tight_layout()

# Copied from Param_Generator.py

#np.random.seed(6) # This is the seed for the Tuning_Params_Extra10k.csv file
#np.random.seed(28) # This is the seed for the Tuning_Params_20k-30k.csv file
np.random.seed(4) # This is the seed for the Tuning_Params_20k-30k.csv file

# N is the number of parameters to be generated.
N = 10000

aice_list = []
asnow_list = []
MF_list = []
for i in range(0,N):

    aice = np.random.normal(0.000003396, 0.00000438) #OG distribution
    while (aice < 0) or (aice > max(aice_passing_4sigma)):
        aice = np.random.normal(0.000003396, 0.00000438) #OG
    aice_list.append(aice)
    
    asnow = np.random.normal(0.000001546,0.00000085) #OG
    while (asnow < 0) or (asnow > max(asnow_passing_4sigma)):
        asnow = np.random.normal(0.000001546,0.00000085) #OG
    asnow_list.append(asnow)
        
    MF = np.random.normal(0.0002707,0.0001632) #OG
    while (MF < 0) or (MF > max(MF_passing_4sigma)):
        MF = np.random.normal(0.0002707,0.0001632) #OG
    MF_list.append(MF)

d = {'aice': aice_list, 'asnow': asnow_list, 'MF': MF_list}
df = pd.DataFrame(data=d)

#df.to_csv('Tuning_Params_20k-30k_Truncated_for_AWSBIAS.csv')

# =============================================================================
# Plot param distirbutions with the truncated 10k params in grey
# =============================================================================

aice = np.concatenate((aice1,aice2))
asnow = np.concatenate((asnow1,asnow2))
MF = np.concatenate((MF1,MF2))

N = int(np.sqrt(len(aice_list)))
#aice_bins = np.arange(0,2.1e-5,0.2e-5)
#asnow_bins = np.arange(0,4.7e-6,0.4e-6)
#MF_bins = np.arange(0,8.9e-4,0.8e-4)

aice_bins = np.linspace(0,2.1e-5,N)
asnow_bins = np.linspace(0,4.7e-6,N)
MF_bins = np.linspace(0,8.9e-4,N)

aicestd = 0.00000438
asnowstd = 0.00000085
MFstd = 0.0001632

fig, (ax) = plt.subplots(ncols=3,figsize=(10,3.5),gridspec_kw={"width_ratios":[1,1,1]},sharey=True)
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.title('a$_{ice}$',fontsize=14)
plt.hist(aice_list,aice_bins,rwidth=1,color='lightgrey',label='All simulations')
#plt.hist(aice_passing_allsims,aice_bins,rwidth=1,color='red',zorder=5,alpha=0.5,label='$\dot{b}_{mod}$ = $\dot{b}_{obs}$ $\pm$ 3$\sigma$')
plt.hist(aice_passing,aice_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5,label='$\dot{b}_{mod}$ = $\dot{b}_{obs}$ $\pm$ 3$\sigma$\nand a$_{ice}$ $\geq$ a$_{snow}$',align='mid')
#plt.hist(aice_passing_allsims,aice_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(aice_passing,aice_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10,align='mid')
plt.ylim(0,350)
plt.ylabel('Frequency',fontsize=14)
plt.vlines(x=0.000003396,ymin=0,ymax=2500,linestyle='--')
plt.vlines(x= max(aice_passing_4sigma),ymin=0,ymax=2500,linestyle='--',color='cyan',linewidth=3,label='Maximum param from\n$\dot{b}_{mod}$ = $\dot{b}_{obs}$ $\pm$ 4$\sigma$\nand a$_{ice}$ $\geq$ a$_{snow}$')
#plt.xticks([0,0.000003396,0.000003396+(1*aicestd),0.000003396+(2*aicestd),0.000003396+(3*aicestd)],['0','3.4x10$^{-6}$','7.8x10$^{-6}$','1.2x10$^{-5}$','1.6x10$^{-5}$'],rotation=45)
plt.xticks([0,0.000003396,0.000003396+(3*aicestd)],['0','3.4x10$^{-6}$','1.6x10$^{-5}$'],rotation=0)
x = np.linspace(np.min(aice), np.max(aice), 100)
p = norm.pdf(x, 0.000003396, 0.00000438)
plt.plot(x, p/np.max(p)*np.max(np.histogram(aice_list,aice_bins)[0]), 'k', linewidth=2,label='Normal distribution')
plt.margins(x=0)

plt.subplot(1,3,2)
plt.title('a$_{snow}$',fontsize=14)
plt.hist(asnow_list,asnow_bins,rwidth=1,color='lightgrey')
#plt.hist(asnow_passing_allsims,asnow_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(asnow_passing,asnow_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
#plt.hist(asnow_passing_allsims,asnow_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(asnow_passing,asnow_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10)
plt.ylim(0,450)
plt.vlines(x=0.000001546,ymin=0,ymax=2500,linestyle='--')
plt.vlines(x= max(asnow_passing_4sigma),ymin=0,ymax=2500,linestyle='--',color='cyan',linewidth=3)
plt.xticks([0,0.000001546,0.000001546+(3*asnowstd)],['0','1.5x10$^{-6}$','4.1x10$^{-6}$'])
x = np.linspace(np.min(asnow), np.max(asnow), 100)
p = norm.pdf(x, 0.000001546,0.00000085)
plt.plot(x, p/np.max(p)*np.max(np.histogram(asnow_list,asnow_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)
plt.subplot(1,3,3)
plt.title('$MF$',fontsize=14)
plt.hist(MF_list,MF_bins,rwidth=1,color='lightgrey')
#plt.hist(MF_passing_allsims,MF_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(MF_passing,MF_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
#plt.hist(MF_passing_allsims,MF_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(MF_passing,MF_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10)
plt.ylim(0,350)
plt.vlines(x=0.0002707,ymin=0,ymax=2500,linestyle='--')
plt.vlines(x= max(MF_passing_4sigma),ymin=0,ymax=2500,linestyle='--',color='cyan',linewidth=3)
plt.xticks([0,0.0002707,0.0002707+(3*MFstd)],['0','2.7x10$^{-4}$','7.6x10$^{-4}$'])
x = np.linspace(np.min(MF), np.max(MF), 100)
p = norm.pdf(x, 0.0002707,0.0001632)
plt.plot(x, p/np.max(p)*np.max(np.histogram(MF_list,MF_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

#handles, labels = ax.get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
fig.legend(fontsize=14,bbox_to_anchor=(1.2,0.8), ncol=1, borderaxespad=0.19)
plt.tight_layout()
