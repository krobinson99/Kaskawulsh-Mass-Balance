# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:10:32 2023

Value added analysis

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import stats

# =============================================================================
# Load the full population of parameters used in the tuning
# =============================================================================
params_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/Tuning/Tuning_Params_Initial10k.csv' # Path to csv file containing MF, aice, asnow.
params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice = params[:,1][0:7000]
asnow = params[:,2][0:7000]
MF = params[:,3][0:7000]
# 
mb0_1000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_0-1000.csv',delimiter=',')
mb1001_2000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_1001-2000.csv',delimiter=',')
mb2001_3000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_2001-3000.csv',delimiter=',')
mb3001_4000 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_3001-4000.csv',delimiter=',')
mb4001_4999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_4001-4999.csv',delimiter=',')
mb5000_5999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_5000-5999.csv',delimiter=',')
mb6000_6999 = np.loadtxt('D:/TuningOutputs/Tuning_Sept23/mb_results_sim_6000-6999.csv',delimiter=',')

mb = np.concatenate((mb0_1000,mb1001_2000,mb2001_3000,mb3001_4000,mb4001_4999,mb5000_5999,mb6000_6999),axis=0)

simID = np.arange(0,len(mb))     

# =============================================================================
# Calculate the stability of the distributions on aice, asnow, MF
# =============================================================================

def Kolmogorov_Smirnov_test(sample_size,full_population):
    
    sample = random.sample(list(full_population),k=sample_size)
    ks_statistic, p_value = stats.ks_2samp(full_population,np.array(sample))
    
    return ks_statistic, p_value, np.array(sample)

def KS_test(stat):
    '''
    stat = which statistic to calculate from the KS test:
        0 = KS statistic
        1 = p-value
    '''
    sample_sizes = np.arange(1000,len(mb)+1,1000)
    results = []
    for sample_size in sample_sizes:
        print('sample size =',sample_size)
        
        # repeat the KS test 100 times
        aice_pvals = []
        asnow_pvals = []
        MF_pvals = []
        mb_pvals = []
        for i in range(1,101):
            aice_pvals.append(Kolmogorov_Smirnov_test(sample_size,aice)[stat])
            asnow_pvals.append(Kolmogorov_Smirnov_test(sample_size,asnow)[stat])
            MF_pvals.append(Kolmogorov_Smirnov_test(sample_size,MF)[stat])
            mb_pvals.append(Kolmogorov_Smirnov_test(sample_size,mb)[stat])
        
        results.extend([np.mean(aice_pvals),np.mean(asnow_pvals),np.mean(MF_pvals),np.mean(mb_pvals)])
    
    KS_results = np.reshape(np.array(results),(len(sample_sizes),4))
    return KS_results

KS_pval_results = KS_test(1)
KS_stat_results = KS_test(0)

def plot_KS_pval(i):
    plt.bar(np.arange(1000,len(mb)+1,1000),KS_pval_results[:,i],align='center',color='navy',width=500)
    plt.xlabel('Sample size',fontsize=14)
    plt.ylabel('p-value',fontsize=14)
    plt.ylim(0.5,1.02)
    plt.grid()

plt.figure(figsize=(7,6))
titles = ['a$_{ice}$','a$_{snow}$','$MF$','Mass Balance']
for i in range(0,4):
    plt.subplot(2,2,i+1)
    plt.title(titles[i],fontsize=14)
    plot_KS_pval(i)
plt.tight_layout()

def plot_KS_scores(i):
    plt.bar(np.arange(1000,len(mb)+1,1000),KS_stat_results[:,i],align='center',color='navy',width=500)
    plt.xlabel('Sample size',fontsize=14)
    plt.ylabel('KS statistic',fontsize=14)
    plt.ylim(0,0.026)
    plt.grid()

plt.figure(figsize=(7,6))
titles = ['a$_{ice}$','a$_{snow}$','$MF$','Mass Balance']
for i in range(0,4):
    plt.subplot(2,2,i+1)
    plt.title(titles[i],fontsize=14)
    plot_KS_scores(i)
plt.tight_layout()


# =============================================================================
# Calculate the mean and std dev of samples of aice, asnow, MF as they
# approach the total population size
# =============================================================================

def params_mean_std(sample_size,full_population):
    sample = random.sample(list(full_population),k=sample_size)
    return np.mean(sample),np.std(sample)

sample_sizes = np.arange(100,len(mb)+1,100)
mean_results, std_results = [], []
for sample_size in sample_sizes:
    print('sample size =',sample_size)
    
    # repeat the KS test 100 times
    aice_means, aice_stds = [], []
    asnow_means, asnow_stds = [], []
    MF_means, MF_stds = [], []
    mb_means, mb_stds = [], []
    for i in range(1,101):
        aice_means.append(params_mean_std(sample_size,aice)[0])
        asnow_means.append(params_mean_std(sample_size,asnow)[0])
        MF_means.append(params_mean_std(sample_size,MF)[0])
        mb_means.append(params_mean_std(sample_size,mb)[0])
        
        aice_stds.append(params_mean_std(sample_size,aice)[1])
        asnow_stds.append(params_mean_std(sample_size,asnow)[1])
        MF_stds.append(params_mean_std(sample_size,MF)[1])
        mb_stds.append(params_mean_std(sample_size,mb)[1])
    
    mean_results.extend([np.mean(aice_means),np.mean(asnow_means),np.mean(MF_means),np.mean(mb_means)])
    std_results.extend([np.mean(aice_stds),np.mean(asnow_stds),np.mean(MF_stds),np.mean(mb_stds)])

mean_results = np.reshape(np.array(mean_results),(len(sample_sizes),4))
std_results = np.reshape(np.array(std_results),(len(sample_sizes),4))

def plot_sample_means(i):
    plt.plot(sample_sizes,mean_results[:,i],color='navy')
    plt.hlines()
    plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    plt.xlabel('Sample size',fontsize=14)
    plt.ylabel('Mean ' + titles[i],fontsize=14)
    plt.grid()

plt.figure(figsize=(7,6))
titles = ['a$_{ice}$','a$_{snow}$','$MF$','Mass Balance']
for i in range(0,4):
    plt.subplot(2,2,i+1)
    plt.title(titles[i],fontsize=14)
    plot_sample_means(i)
plt.tight_layout()

def plot_sample_mean_diffs(i):
    plt.plot(sample_sizes[1:],np.abs(np.diff(mean_results[:,i])),color='navy')
    plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    plt.xlabel('Sample size',fontsize=14)
    plt.ylabel('Absolute difference',fontsize=14)
    plt.grid()

plt.figure(figsize=(7,6))
titles = ['a$_{ice}$','a$_{snow}$','$MF$','Mass Balance']
for i in range(0,4):
    plt.subplot(2,2,i+1)
    plt.title(titles[i],fontsize=14)
    plot_sample_mean_diffs(i)
plt.tight_layout()

def plot_sample_stds(i):
    plt.plot(sample_sizes,std_results[:,i],color='navy')
    plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    plt.xlabel('Sample size',fontsize=14)
    plt.ylabel('std. dev.',fontsize=14)
    plt.grid()

plt.figure(figsize=(7,6))
titles = ['a$_{ice}$','a$_{snow}$','$MF$','Mass Balance']
for i in range(0,4):
    plt.subplot(2,2,i+1)
    plt.title(titles[i],fontsize=14)
    plot_sample_stds(i)
plt.tight_layout()

def plot_sample_std_diffs(i):
    plt.plot(sample_sizes[1:],np.abs(np.diff(std_results[:,i])),color='navy')
    plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    plt.xlabel('Sample size',fontsize=14)
    plt.ylabel('Absolute difference',fontsize=14)
    plt.grid()

plt.figure(figsize=(7,6))
titles = ['a$_{ice}$','a$_{snow}$','$MF$','Mass Balance']
for i in range(0,4):
    plt.subplot(2,2,i+1)
    plt.title(titles[i],fontsize=14)
    plot_sample_std_diffs(i)
plt.tight_layout()


