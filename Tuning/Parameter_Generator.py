# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:45:10 2023

Generate N random combinations of a_snow, a_ice, MF, and save them to a csv

@author: katierobinson
"""

import numpy as np
import pandas as pd

# N is the number of parameters to be generated.
N = 10000

aice_list = []
asnow_list = []
MF_list = []
for i in range(0,N):

    aice = np.random.normal(0.000003396, 0.00000438) #OG distribution
    while aice < 0:
        aice = np.random.normal(0.000003396, 0.00000438) #OG
    aice_list.append(aice)
    
    asnow = np.random.normal(0.000001546,0.00000085) #OG
    while asnow < 0:
        asnow = np.random.normal(0.000001546,0.00000085) #OG
    asnow_list.append(asnow)
        
    MF = np.random.normal(0.0002707,0.0001632) #OG
    while MF < 0:
        MF = np.random.normal(0.0002707,0.0001632) #OG
    MF_list.append(MF)

d = {'aice': aice_list, 'asnow': asnow_list, 'MF': MF_list}
df = pd.DataFrame(data=d)

#df.to_csv('Tuning_Params_10k.csv')
