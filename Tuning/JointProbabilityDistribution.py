# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 10:41:52 2022

develop a joint probability distribution on a_snow, a_ice, MF

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.stats.stats import pearsonr 
from scipy.stats import multivariate_normal


#start by plotting asnow vs aice vs MF: colour code based on resulting NET MB
# make one plot for debris case and one for debris-free

#load param combos from debris and debris free cases
debris_params = 'F:\\Mass Balance Model\\Kaskawulsh-Mass-Balance\\RunModel\\final_params_deb.csv'
params = np.loadtxt(debris_params) #namelist!!
aice_deb = params[0,:]
asnow_deb = params[1,:]
MF_deb = params[2,:]

nondebris_params = 'F:\\Mass Balance Model\\Kaskawulsh-Mass-Balance\\RunModel\\final_params.csv'
params = np.loadtxt(nondebris_params) #namelist!!
aice_nodeb = params[0,:]
asnow_nodeb = params[1,:]
MF_nodeb = params[2,:]

MFmean, MFstd = 2.707e-4, 1.632e-4
aicemean, aicestd = 3.396e-6, 2.65e-6
asnowmean, asnowstd = 1.546e-6, 0.85e-6
#aice = np.random.normal(0.000003396, 0.00000438)

#get the net mass balance for each of the 12 / 25 runs 
deb_mbs_npy = np.load('F:/Mass Balance Model/OUTPUTS/Kaskonly_Baseline/Plots/allrunMBs_Kaskonly_Baseline.npy')
deb_mbs = []
for i in range(0,len(deb_mbs_npy)):
    deb_mbs.append(np.nanmean(deb_mbs_npy[i]))
    
nondeb_mbs_npy = np.load('F:/Mass Balance Model/OUTPUTS/Plots/Baseline_NoDebris_Allparams/allrunMBs_Baseline_NoDebris_Allparams.npy')
nondeb_mbs = []
for i in range(0,len(nondeb_mbs_npy)):
    nondeb_mbs.append(np.nanmean(nondeb_mbs_npy[i]))


###############################################################################
###############################################################################
###############################################################################

#recreate the 3x3 plots from Surjanovic, 2016 (https://summit.sfu.ca/item/16642)
#debris case first
#plt.figure(figsize=(9,9))
fig, (ax, ax2, cax) = plt.subplots(ncols=3,figsize=(9,9), 
                  gridspec_kw={"width_ratios":[1,1, 0.05]})
fig.subplots_adjust(wspace=0.3)
plt.suptitle('Model parameters - debris case',fontsize=14,y=1.02)
plt.subplot(3,3,1)
plt.text(0.45,0.5,'MF',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,2)
im = plt.scatter(asnow_deb,MF_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60) #,levels = np.linspace(-0.63,-0.29, 20))
plt.xlim(np.min(asnow_deb)-0.1e-6,np.max(asnow_deb)+0.1e-6)
plt.ylim(np.min(MF_deb)-0.1e-4,np.max(MF_deb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,3)
plt.scatter(aice_deb,MF_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(aice_deb)-0.1e-6,np.max(aice_deb)+0.1e-6)
plt.ylim(np.min(MF_deb)-0.1e-4,np.max(MF_deb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
#plt.tick_params(bottom=False, top=True, left=False, right=True)
#plt.tick_params(labelbottom=False, labeltop=True, labelleft=False, labelright=True)
plt.subplot(3,3,4)
plt.scatter(MF_deb,asnow_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(asnow_deb)-0.1e-6,np.max(asnow_deb)+0.1e-6)
plt.xlim(np.min(MF_deb)-0.1e-4,np.max(MF_deb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,5)
plt.text(0.35,0.5,'asnow',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,6)
plt.scatter(aice_deb,asnow_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(aice_deb)-0.1e-6,np.max(aice_deb)+0.1e-6)
plt.ylim(np.min(asnow_deb)-0.1e-6,np.max(asnow_deb)+0.1e-6)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,7)
plt.scatter(MF_deb,aice_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(aice_deb)-0.1e-6,np.max(aice_deb)+0.1e-6)
plt.xlim(np.min(MF_deb)-0.1e-4,np.max(MF_deb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,8)
plt.scatter(asnow_deb,aice_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(aice_deb)-0.1e-6,np.max(aice_deb)+0.1e-6)
plt.xlim(np.min(asnow_deb)-0.1e-6,np.max(asnow_deb)+0.1e-6)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,9)
plt.text(0.4,0.5,'aice',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
cb_ax = fig.add_axes([1.03, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
cbar.ax.set_ylabel('Net Mass Balance (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=23)
plt.tight_layout()
#plt.savefig('Meltmodelparams_debcase.png',bbox_inches = 'tight')

###############################################################################
###############################################################################
###############################################################################

fig, (ax, ax2, cax) = plt.subplots(ncols=3,figsize=(9,9), 
                  gridspec_kw={"width_ratios":[1,1, 0.05]})
fig.subplots_adjust(wspace=0.3)
plt.suptitle('Model parameters - non-debris case',fontsize=14,y=1.02)
plt.subplot(3,3,1)
plt.text(0.45,0.5,'MF',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,2)
im = plt.scatter(asnow_nodeb,MF_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(asnow_nodeb)-0.1e-6,np.max(asnow_nodeb)+0.1e-6)
plt.ylim(np.min(MF_nodeb)-0.1e-4,np.max(MF_nodeb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,3)
plt.scatter(aice_nodeb,MF_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(aice_nodeb)-0.1e-6,np.max(aice_nodeb)+0.1e-6)
plt.ylim(np.min(MF_nodeb)-0.1e-4,np.max(MF_nodeb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
#plt.tick_params(bottom=False, top=True, left=False, right=True)
#plt.tick_params(labelbottom=False, labeltop=True, labelleft=False, labelright=True)
plt.subplot(3,3,4)
plt.scatter(MF_nodeb,asnow_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(asnow_nodeb)-0.1e-6,np.max(asnow_nodeb)+0.1e-6)
plt.xlim(np.min(MF_nodeb)-0.1e-4,np.max(MF_nodeb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,5)
plt.text(0.35,0.5,'asnow',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,6)
plt.scatter(aice_nodeb,asnow_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(aice_nodeb)-0.1e-6,np.max(aice_nodeb)+0.1e-6)
plt.ylim(np.min(asnow_nodeb)-0.1e-6,np.max(asnow_nodeb)+0.1e-6)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,7)
plt.scatter(MF_nodeb,aice_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(aice_nodeb)-0.1e-6,np.max(aice_nodeb)+0.1e-6)
plt.xlim(np.min(MF_nodeb)-0.1e-4,np.max(MF_nodeb)+0.1e-4)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,8)
plt.scatter(asnow_nodeb,aice_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(aice_nodeb)-0.1e-6,np.max(aice_nodeb)+0.1e-6)
plt.xlim(np.min(asnow_nodeb)-0.1e-6,np.max(asnow_nodeb)+0.1e-6)
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,9)
plt.text(0.4,0.5,'aice',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
cb_ax = fig.add_axes([1.03, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
cbar.ax.set_ylabel('Net Mass Balance (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=23)
plt.tight_layout()
#plt.savefig('Meltmodelparams_nodebcase.png',bbox_inches = 'tight')


#######################################################################
# calculate the correlation matrix:
#deb
correlation_snowMF = pearsonr(asnow_deb,MF_deb)[0]
correlation_iceMF = pearsonr(aice_deb,MF_deb)[0]
correlation_snowice = pearsonr(aice_deb,asnow_deb)[0]

#no deb
correlation_snowMF = pearsonr(asnow_nodeb,MF_nodeb)[0]
correlation_iceMF = pearsonr(aice_nodeb,MF_nodeb)[0]
correlation_snowice = pearsonr(aice_nodeb,asnow_nodeb)[0]

#multivariate guassian distribution:
# https://numpy.org/doc/stable/reference/random/generated/numpy.random.multivariate_normal.html

#get covariance matrix: m has 1 row per variable and one column per value (ie. 3 x 25 for non deb case)
m_deb = np.array([MF_deb,aice_deb,asnow_deb])
cov_deb = np.cov(m_deb)

means_deb = np.array([np.mean(MF_deb),np.mean(aice_deb),np.mean(asnow_deb)])
a_deb = np.random.multivariate_normal(means_deb, cov_deb, 1000).T
MF_deb_jpd = a_deb[0]
aice_deb_jpd = a_deb[1]
asnow_deb_jpd = a_deb[2]




fig, (ax, ax2, cax) = plt.subplots(ncols=3,figsize=(9,9), 
                  gridspec_kw={"width_ratios":[1,1, 0.05]})
fig.subplots_adjust(wspace=0.3)
plt.suptitle('Multivariate normal distribution - debris case',fontsize=14,y=1.02)
plt.subplot(3,3,1)
plt.text(0.45,0.5,'MF',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,2)
plt.scatter(asnow_deb_jpd, MF_deb_jpd,zorder=-1,alpha=0.3)
im = plt.scatter(asnow_deb,MF_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(asnow_deb_jpd),np.max(asnow_deb_jpd))
plt.ylim(np.min(MF_deb_jpd),np.max(MF_deb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,3)
plt.scatter(aice_deb_jpd, MF_deb_jpd,zorder=-1,alpha=0.3)
plt.scatter(aice_deb,MF_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(aice_deb_jpd),np.max(aice_deb_jpd))
plt.ylim(np.min(MF_deb_jpd),np.max(MF_deb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
#plt.tick_params(bottom=False, top=True, left=False, right=True)
#plt.tick_params(labelbottom=False, labeltop=True, labelleft=False, labelright=True)
plt.subplot(3,3,4)
plt.scatter(MF_deb_jpd, asnow_deb_jpd,zorder=-1,alpha=0.3)
plt.scatter(MF_deb,asnow_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(asnow_deb_jpd),np.max(asnow_deb_jpd))
plt.xlim(np.min(MF_deb_jpd),np.max(MF_deb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,5)
plt.text(0.35,0.5,'asnow',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,6)
plt.scatter(aice_deb_jpd, asnow_deb_jpd,zorder=-1,alpha=0.3)
plt.scatter(aice_deb,asnow_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(asnow_deb_jpd),np.max(asnow_deb_jpd))
plt.xlim(np.min(aice_deb_jpd),np.max(aice_deb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,7)
plt.scatter(MF_deb_jpd, aice_deb_jpd,zorder=-1,alpha=0.3)
plt.scatter(MF_deb,aice_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(aice_deb_jpd),np.max(aice_deb_jpd))
plt.xlim(np.min(MF_deb_jpd),np.max(MF_deb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,8)
plt.scatter(asnow_deb_jpd, aice_deb_jpd,zorder=-1,alpha=0.3)
plt.scatter(asnow_deb,aice_deb,c=deb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(asnow_deb_jpd),np.max(asnow_deb_jpd))
plt.ylim(np.min(aice_deb_jpd),np.max(aice_deb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,9)
plt.text(0.4,0.5,'aice',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
cb_ax = fig.add_axes([1.03, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
cbar.ax.set_ylabel('Net Mass Balance (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=23)
plt.tight_layout()
#plt.savefig('Meltmodelparams_debcase_multivariatenorm.png',bbox_inches = 'tight')

#########################################################################
m_nodeb = np.array([MF_nodeb,aice_nodeb,asnow_nodeb])
cov_nodeb = np.cov(m_nodeb)

means_nodeb = np.array([np.mean(MF_nodeb),np.mean(aice_nodeb),np.mean(asnow_nodeb)])
a_nodeb = np.random.multivariate_normal(means_nodeb, cov_nodeb, 1000).T
MF_nodeb_jpd = a_nodeb[0]
aice_nodeb_jpd = a_nodeb[1]
asnow_nodeb_jpd = a_nodeb[2]

fig, (ax, ax2, cax) = plt.subplots(ncols=3,figsize=(9,9), 
                  gridspec_kw={"width_ratios":[1,1, 0.05]})
fig.subplots_adjust(wspace=0.3)
plt.suptitle('Multivariate normal distribution - non-debris case',fontsize=14,y=1.02)
plt.subplot(3,3,1)
plt.text(0.45,0.5,'MF',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,2)
plt.scatter(asnow_nodeb_jpd, MF_nodeb_jpd,zorder=-1,alpha=0.3)
im = plt.scatter(asnow_nodeb,MF_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(asnow_nodeb_jpd),np.max(asnow_nodeb_jpd))
plt.ylim(np.min(MF_nodeb_jpd),np.max(MF_nodeb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,3)
plt.scatter(aice_nodeb_jpd, MF_nodeb_jpd,zorder=-1,alpha=0.3)
plt.scatter(aice_nodeb,MF_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(aice_nodeb_jpd),np.max(aice_nodeb_jpd))
plt.ylim(np.min(MF_nodeb_jpd),np.max(MF_nodeb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
#plt.tick_params(bottom=False, top=True, left=False, right=True)
#plt.tick_params(labelbottom=False, labeltop=True, labelleft=False, labelright=True)
plt.subplot(3,3,4)
plt.scatter(MF_nodeb_jpd, asnow_nodeb_jpd,zorder=-1,alpha=0.3)
plt.scatter(MF_nodeb,asnow_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(asnow_nodeb_jpd),np.max(asnow_nodeb_jpd))
plt.xlim(np.min(MF_nodeb_jpd),np.max(MF_nodeb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,5)
plt.text(0.35,0.5,'asnow',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.subplot(3,3,6)
plt.scatter(aice_nodeb_jpd, asnow_nodeb_jpd,zorder=-1,alpha=0.3)
plt.scatter(aice_nodeb,asnow_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(asnow_nodeb_jpd),np.max(asnow_nodeb_jpd))
plt.xlim(np.min(aice_nodeb_jpd),np.max(aice_nodeb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,7)
plt.scatter(MF_nodeb_jpd, aice_nodeb_jpd,zorder=-1,alpha=0.3)
plt.scatter(MF_nodeb,aice_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.ylim(np.min(aice_nodeb_jpd),np.max(aice_nodeb_jpd))
plt.xlim(np.min(MF_nodeb_jpd),np.max(MF_nodeb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,8)
plt.scatter(asnow_nodeb_jpd, aice_nodeb_jpd,zorder=-1,alpha=0.3)
plt.scatter(asnow_nodeb,aice_nodeb,c=nondeb_mbs,cmap='RdYlBu',vmin=-0.63,vmax=-0.29,edgecolor='k',s=60)
plt.xlim(np.min(asnow_nodeb_jpd),np.max(asnow_nodeb_jpd))
plt.ylim(np.min(aice_nodeb_jpd),np.max(aice_nodeb_jpd))
pylab.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
plt.subplot(3,3,9)
plt.text(0.4,0.5,'aice',fontsize=15)
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
cb_ax = fig.add_axes([1.03, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
cbar.ax.set_ylabel('Net Mass Balance (m w.e. $a^{-1}$)', rotation=270,fontsize=14,labelpad=23)
plt.tight_layout()
#plt.savefig('Meltmodelparams_nodebcase_multivariatenorm.png',bbox_inches = 'tight')

#plot the original 1000 param 1-d distributions + the new multivariate distributions on top 
#goal is to show that we can test less parameter combos 

fig = plt.figure(figsize=(12,5))
plt.subplot(1,3,1)
mu, sigma = aicemean, aicestd # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)
count, bins, ignored = plt.hist(s, 30, density=False, color='k')
count2, bins2, ignored = plt.hist(aice_deb_jpd, 30, density=False)
pylab.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
plt.xlim(0)
plt.xlabel('a$_{ice}$',fontsize=16)
plt.ylabel('Frequency',fontsize=14)
plt.subplot(1,3,2)
mu, sigma = asnowmean, asnowstd # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)
count, bins, ignored = plt.hist(s, 30, density=False, color='k')
count2, bins2, ignored = plt.hist(asnow_deb_jpd, 30, density=False)
pylab.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
plt.xlim(0)
plt.xlabel('a$_{snow}$',fontsize=16)
plt.ylabel('Frequency',fontsize=14)
plt.subplot(1,3,3)
mu, sigma = MFmean, MFstd # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)
count, bins, ignored = plt.hist(s, 30, density=False, color='k')
count2, bins2, ignored = plt.hist(MF_nodeb_jpd, 30, density=False)
pylab.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
plt.xlim(0)
plt.xlabel('MF',fontsize=16)
plt.ylabel('Frequency',fontsize=14)
plt.tight_layout()
fig.legend(['Full parameter distribution','Joint probability distribution (debris case)'],fontsize=14,bbox_to_anchor=(0.85, 1.05),ncol=2)
#plt.savefig('Param_distributions',bbox_inches = 'tight')
