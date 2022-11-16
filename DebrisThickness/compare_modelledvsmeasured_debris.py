# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 12:10:12 2022

Compare what the melt factors would have been if we used the measured debris thicknesses
in each gridcell compared to the modelled debris thicknesses 

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
from generate_meltfactors import dirtyicemelt_obs,peakmelt_obs,peakthickness_ref,transition_ref
sys.path.insert(1,'F:\Mass Balance Model\Kaskawulsh-Mass-Balance\FieldWork2022')
from debrissurvey import dr_estimate_avg as MODELLED_debris
from debrissurvey import debthickness_obs_avg as MEASURED_debris
from debrissurvey import reachedice_avg as reached_ice


#params
peakM = 6.62306028549649 # from hdopt_params.csv (melt_mwea_2cm)
b0 = 11.0349260206858 #from hdopt_params.csv (b0)
k = 1.98717418666925 # from hdopt_params.csv
cleaniceM = 2.9200183038347
peakthickness = 0.02 #m

#transect A = [:18]
#transect B = [18:38]
#transect C = [38:]

#determine the melt factor assosciated with each MEASURED and MODELLED debris thickness

def meltfactor(dh,cleaniceM=dirtyicemelt_obs,peakM=peakmelt_obs,peakM_thickness=peakthickness_ref,transition_thickness=transition_ref,b0 = 11.0349260206858,k = 1.98717418666925):
    '''
    #input a value for debris thickness (dh; cm) and return the melt factor based on the other values
    #function should take in:
        # debris val,
        # cleanicemelt and peak melt value
        # transition thickness and thickness at peak melt 
        # Rounce et al. model params for the kaskawulsh
    #function should return:
        # 1 val for melt factor. should be 1 for debris_h = 0
    '''
    # convert h from cm to m:
    h = dh/100                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    if h == 0:
        meltfactor = 1
    elif h <= peakM_thickness:
        Melt = (((peakM - cleaniceM)/peakM_thickness)*h) + cleaniceM
        meltfactor = Melt/cleaniceM
    elif (h > peakM_thickness) and (h <= transition_thickness):
        yintercept = peakM - (((cleaniceM-peakM)/(transition_thickness-peakM_thickness))*peakM_thickness)
        Melt = (((cleaniceM-peakM)/(transition_thickness-peakM_thickness))*h) + yintercept
        meltfactor = Melt/cleaniceM
    elif h > transition_thickness:
        #determine how much this part of the curve needs to be shifted over to meet the linear parts
        h_curvestart = (b0-cleaniceM)/(cleaniceM*k*b0)
        horizontal_shift = transition_thickness - h_curvestart
        Melt = b0/(1+(k*b0*(h-horizontal_shift))) #calculate the melt on the original curve(take away the shift) so that you're actually getting the melt at the desired thickness
        meltfactor = Melt/cleaniceM
    else:
        print('debris thickness not accounted for in if statements')

    return meltfactor 
    
h = 0.006
MF = meltfactor(h)
print(MF)

MODELLED_MF = []
MEASURED_MF = []
for i in range(0,len(MODELLED_debris[:49])):
    MODELLED_MF.append(meltfactor(MODELLED_debris[i]))
    MEASURED_MF.append(meltfactor(MEASURED_debris[i]))
    
plt.figure(figsize=(6,6))
plt.scatter(MEASURED_MF[:18],MODELLED_MF[:18],color='mediumblue',label='Transect A',s=60,edgecolor='k')
plt.scatter(MEASURED_MF[18:38],MODELLED_MF[18:38],color='darkorange',label='Transect B',s=60,edgecolor='k')
plt.scatter(MEASURED_MF[38:],MODELLED_MF[38:],color='crimson',label='Transect C',s=60,edgecolor='k')
plt.scatter(np.array(MEASURED_MF[:18])[reached_ice[:18]==0],np.array(MODELLED_MF[:18])[reached_ice[:18]==0],color='white',s=60,edgecolor='mediumblue')
plt.scatter(np.array(MEASURED_MF[18:38])[reached_ice[18:38]==0],np.array(MODELLED_MF[18:38])[reached_ice[18:38]==0],color='white',s=60,edgecolor='darkorange')
plt.scatter(np.array(MEASURED_MF[38:])[reached_ice[38:49]==0],np.array(MODELLED_MF[38:])[reached_ice[38:49]==0],color='white',s=60,edgecolor='crimson')
plt.scatter(MEASURED_MF[0],MODELLED_MF[0],color='white',edgecolor='k',s=60,zorder=-1,label='Measurement did not\nreach ice surface') #fake point used to add empty circle to legend, zorder hides it behind a blue point
plt.legend(fontsize=12)
plt.xlabel('Melt factor of measured debris',fontsize=14)
plt.ylabel('Melt factor of modelled debris',fontsize=14)
plt.xlim(0,1.1)
plt.ylim(0,1.1)
plt.plot([0,1.1],[0,1.1],color='grey',linestyle='--',zorder=-2)
#plt.savefig('modelledvsmeasured_MF.png',bbox_inches = 'tight')

