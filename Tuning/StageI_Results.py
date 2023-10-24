# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 11:33:55 2023

Plot the results from Stage 1 of tuning

@author: katierobinson
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import pandas as pd
from scipy.stats import norm
import sys
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\RunModel')
from Model_functions_ver4 import shiftedColorMap

# Get the tuning results:
# =============================================================================
params_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Tuning_Params_Initial10k.csv' # Path to csv file containing MF, aice, asnow.
params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice = params[:,1][0:5000]
asnow = params[:,2][0:5000]
MF = params[:,3][0:5000]
# 
mb0_1000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_0-1000.csv',delimiter=',')
mb1001_2000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_1001-2000.csv',delimiter=',')
mb2001_3000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_2001-3000.csv',delimiter=',')
mb3001_4000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_3001-4000.csv',delimiter=',')
mb4001_4999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_4001-4999.csv',delimiter=',')
mb = np.concatenate((mb0_1000,mb1001_2000,mb2001_3000,mb3001_4000,mb4001_4999),axis=0)
#      

mb_passing = mb[np.where((mb >= (-0.46-0.17)) & (mb <= (-0.46+0.17)))]
aice_passing = aice[np.where((mb >= (-0.46-0.17)) & (mb <= (-0.46+0.17)))]
asnow_passing = asnow[np.where((mb >= (-0.46-0.17)) & (mb <= (-0.46+0.17)))]
MF_passing = MF[np.where((mb >= (-0.46-0.17)) & (mb <= (-0.46+0.17)))]

# =============================================================================
# Save params that pass tuning stage I in new csv:
# =============================================================================
d = {'aice': aice_passing, 'asnow': asnow_passing, 'MF': MF_passing}
df = pd.DataFrame(data=d)
#df.to_csv('Params_Passed_StageI.csv')

# =============================================================================
# Plot each param against the other params
# =============================================================================
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.916666666,stop=1,name='massbal')

fig, (ax, ax2, cax) = plt.subplots(ncols=3,figsize=(9,9), 
                  gridspec_kw={"width_ratios":[1,1, 0.05]})
fig.subplots_adjust(wspace=0.3)
#plt.suptitle('Tuning parameters',fontsize=14,y=1.02)

plt.subplot(3,3,1)
plt.text(0.45,0.5,'$MF$',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)

plt.subplot(3,3,2)
#im = plt.scatter(asnow,MF,c=mb,cmap='massbal',edgecolor='k',s=40,vmin=-0.63,vmax=-0.29
im = plt.scatter(asnow,MF,c='lightgrey',s=40)
plt.scatter(asnow_passing,MF_passing,c=mb_passing,cmap='RdYlBu',edgecolor='k',s=40,vmin=-0.63,vmax=-0.29)
plt.xlim(np.min(asnow)-0.3e-6,np.max(asnow)+0.3e-6)
plt.ylim(np.min(MF)-0.3e-4,np.max(MF)+0.3e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)

plt.subplot(3,3,3)
plt.scatter(aice,MF,c='lightgrey',s=40)
im = plt.scatter(aice_passing,MF_passing,c=mb_passing,cmap='RdYlBu',edgecolor='k',s=40,vmin=-0.63,vmax=-0.29)
plt.xlim(np.min(aice)-0.4e-6,np.max(aice)+0.4e-6)
plt.ylim(np.min(MF)-0.3e-4,np.max(MF)+0.3e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)

plt.subplot(3,3,4)
plt.scatter(MF,asnow,c='lightgrey',s=40)
plt.scatter(MF_passing,asnow_passing,c=mb_passing,cmap='RdYlBu',edgecolor='k',s=40,vmin=-0.63,vmax=-0.29)
plt.ylim(np.min(asnow)-0.3e-6,np.max(asnow)+0.3e-6)
plt.xlim(np.min(MF)-0.3e-4,np.max(MF)+0.3e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)

plt.subplot(3,3,5)
plt.text(0.35,0.5,'a$_{snow}$',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)

plt.subplot(3,3,6)
plt.scatter(aice,asnow,c='lightgrey',s=40)
plt.scatter(aice_passing,asnow_passing,c=mb_passing,cmap='RdYlBu',edgecolor='k',s=40,vmin=-0.63,vmax=-0.29)
plt.xlim(np.min(aice)-0.4e-6,np.max(aice)+0.4e-6)
plt.ylim(np.min(asnow)-0.3e-6,np.max(asnow)+0.3e-6)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)

plt.subplot(3,3,7)
plt.scatter(MF,aice,c='lightgrey',s=40)
plt.scatter(MF_passing,aice_passing,c=mb_passing,cmap='RdYlBu',edgecolor='k',s=40,vmin=-0.63,vmax=-0.29)
plt.ylim(np.min(aice)-0.4e-6,np.max(aice)+0.4e-6)
plt.xlim(np.min(MF)-0.3e-4,np.max(MF)+0.3e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)

plt.subplot(3,3,8)
plt.scatter(asnow,aice,c='lightgrey',s=40)
plt.scatter(asnow_passing,aice_passing,c=mb_passing,cmap='RdYlBu',edgecolor='k',s=40,vmin=-0.63,vmax=-0.29)
plt.ylim(np.min(aice)-0.4e-6,np.max(aice)+0.4e-6)
plt.xlim(np.min(asnow)-0.3e-6,np.max(asnow)+0.3e-6)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.xticks(fontsize=14);plt.yticks(fontsize=14)

plt.subplot(3,3,9)
plt.text(0.4,0.5,'a$_{ice}$',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)

cb_ax = fig.add_axes([1.03, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
cbar.ax.set_ylabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=23)
cbar.ax.tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('TuningParams_Stage1.png',bbox_inches='tight')


# =============================================================================
# Plot histogram of params that were tested
# =============================================================================
aice_bins = np.arange(0,2.1e-5,0.2e-5)
asnow_bins = np.arange(0,4.7e-6,0.4e-6)
MF_bins = np.arange(0,8.9e-4,0.8e-4)

fig, (ax) = plt.subplots(ncols=3,figsize=(12,5),gridspec_kw={"width_ratios":[1,1,1]},sharey=True)
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.title('a$_{ice}$',fontsize=14)
plt.hist(aice,aice_bins,rwidth=0.9,color='lightgrey',label='All parameters')
plt.hist(aice_passing,aice_bins,rwidth=0.9,color='mediumblue',label='Parameters within\ntarget MB')
plt.ylim(0,1200)
plt.ylabel('Frequency',fontsize=14)
plt.vlines(x=0.000003396,ymin=0,ymax=1200,linestyle='--')
plt.xticks([0,0.000003396],['0','3.4x10$^{-6}$'])
x = np.linspace(np.min(aice), np.max(aice), 100)
p = norm.pdf(x, 0.000003396, 0.00000438)
plt.plot(x, p/np.max(p)*np.max(np.histogram(aice,aice_bins)[0]), 'k', linewidth=2,label='Normal distribution')
plt.margins(x=0)

plt.subplot(1,3,2)
plt.title('a$_{snow}$',fontsize=14)
plt.hist(asnow,asnow_bins,rwidth=0.9,color='lightgrey')
plt.hist(asnow_passing,asnow_bins,rwidth=0.9,color='mediumblue')
plt.ylim(0,1200)
plt.vlines(x=0.000001546,ymin=0,ymax=1200,linestyle='--')
plt.xticks([0,0.000001546],['0','1.5x10$^{-6}$'])
x = np.linspace(np.min(asnow), np.max(asnow), 100)
p = norm.pdf(x, 0.000001546,0.00000085)
plt.plot(x, p/np.max(p)*np.max(np.histogram(asnow,asnow_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

plt.subplot(1,3,3)
plt.title('$MF$',fontsize=14)
plt.hist(MF,MF_bins,rwidth=0.9,color='lightgrey')
plt.hist(MF_passing,MF_bins,rwidth=0.9,color='mediumblue')
plt.ylim(0,1200)
plt.vlines(x=0.0002707,ymin=0,ymax=1200,linestyle='--')
plt.xticks([0,0.0002707],['0','2.7x10$^{-4}$'])
x = np.linspace(np.min(MF), np.max(MF), 100)
p = norm.pdf(x, 0.0002707,0.0001632)
plt.plot(x, p/np.max(p)*np.max(np.histogram(MF,MF_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

#handles, labels = ax.get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
fig.legend(fontsize=14,bbox_to_anchor=(0.71,1.1), ncol=3, borderaxespad=0.19)
plt.tight_layout()
#plt.savefig('TuningParams_Stage1_histograms.png',bbox_inches='tight')

# =============================================================================
# Plot distribution (histogram) of 2007-2018 mass balances from sims that passed tuning
# =============================================================================
mb_passing_bins =  np.arange(-0.46-0.17,-0.28,0.02)

plt.figure(figsize=(12,5))
plt.title('Runs that pass tuning Stage I\nTotal = ' + str(len(mb_passing)),fontsize=14)
plt.hist(mb_passing,mb_passing_bins,rwidth=0.9,color='mediumblue')
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
plt.xticks(mb_passing_bins,np.round(mb_passing_bins,2),fontsize=14,rotation=45)
plt.yticks(fontsize=14)
plt.grid()
plt.tight_layout()
#plt.savefig('TuningParams_Stage1_MBdistribution.png',bbox_inches='tight')

# =============================================================================
# Select mass balance results from normal distribution

# x[np.where((p/np.max(p)*11)>=1)] # 11 is the number of sims in the bin centered on -0.46
# approx range of mb_vals inside normal distribution is -1 to 0.1
# =============================================================================
mb_normaldist_bins =  np.arange(-0.99,0.12,0.02)
mb_normaldist = mb[np.where((mb >= (-0.99)) & (mb <= (0.11)))]

plt.figure(figsize=(12,5))
plt.grid(zorder=20)
plt.title('Stage I',fontsize=14)
plt.hist(mb_normaldist,mb_normaldist_bins,rwidth=0.9,color='mediumblue',zorder=10)
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
plt.xticks(np.arange(-0.98,0.12,0.04),np.round(np.arange(-0.98,0.12,0.04),2),fontsize=14,rotation=45)
plt.yticks(np.arange(0,10,1),fontsize=14)
#x = np.linspace(np.min(mb_passing), np.max(mb_passing), 100)
x = np.linspace(-1,0.12,100)
p = norm.pdf(x, -0.46,0.17)
plt.plot(x, p/np.max(p)*np.histogram(mb_passing,mb_passing_bins)[0][8], 'k', linewidth=2,zorder=12)
plt.tight_layout()

