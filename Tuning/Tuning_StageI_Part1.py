# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 10:06:05 2023

Calculate 2007--2018 mass balance from the tuning outputs
DEMs used to calculate the 2007--2018 mass balance (Young et al. 2021) were from
(primarily) Sept 3 2007 and Oct 1 2018. 

We calculate the mass balance for this period and divide by the number of years (~11.09)
for each simulation.

The 2007-2018 modelled mass balance is saved in a txt file.

@author: katierobinson
"""

import numpy as np
from netCDF4 import Dataset
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd

INPUT_PATH = 'D:/TuningOutputs/Tuning_Sept23'   # Path to tuning outputs

sim = 1
# Get sim number from job array
#sim = int(sys.argv[1])
#print('this is sim',sim)
#print(type(sim))

KRH_tributaries = np.loadtxt('F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Tributaries.txt')

print('sim',sim)
sys.stdout.flush()
sim_mbvals = []
running_mb = np.zeros(KRH_tributaries.shape)
for year in range(2007,2018+1):
    print(year)
    sys.stdout.flush()
    dataset = Dataset(os.path.join(INPUT_PATH,'Netbalance_KRH' + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
    mb = dataset.variables['Net balance'][:]
    sys.stdout.flush()
    dataset.close()
    
    if year == 2007:
        dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq=str(3)+'H')
        DEMdate = np.where(dates == pd.Timestamp(str(year)+'-09-03T00'))[0][0]
        endofyear_mb = np.sum(mb[DEMdate:],axis=0)
        
    elif year == 2018:
        dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq=str(3)+'H')
        DEMdate = np.where(dates == pd.Timestamp(str(year)+'-10-01T00'))[0][0]
        endofyear_mb = np.sum(mb[:DEMdate],axis=0)
        
    else:
        endofyear_mb = np.sum(mb,axis=0)
        
    running_mb += endofyear_mb
    
num_years = len(pd.date_range(start= str(2007) + '-09-03 00:00:00',end= str(2018) + '-10-01 00:00:00',freq='D'))/365
annual_mb = running_mb/num_years
        
KW_mb = np.mean(annual_mb[np.isfinite(KRH_tributaries)])
 
np.savetxt(os.path.join(INPUT_PATH,'sim' + str(sim) + '_mb.txt'),[KW_mb])

